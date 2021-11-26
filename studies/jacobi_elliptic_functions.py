import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy
import scipy
import scipy.special

# EllipK
def rf_duplicate_v2(x, y, z, lam = 0, Q = None):
    Am = (x + y + z)/3
    A = Am + lam
    X = (Am - x) / A
    Y = (Am - y) / A
    Z = (Am - z) / A
    E2 = X*Y + X*Z + Y*Z
    E3 = X*Y*Z
    if not Q:
        r = numpy.finfo(float).eps / 4
        Q = (3*r)**(-1/8) * max(abs(A - x), abs(A - y), abs(A - z))
    if Q < A:
        return 1/math.sqrt(A) * (E2*E2*E3/16 - 5*E2*E2*E2/208 + 3*E3*E3/104 - 3*E2*E3/44 + E2*E2/24 + E3/14 - E2/10 + 1)

    l1 = math.sqrt(x*y + (x+y)*lam + lam*lam)
    l2 = math.sqrt(x*z + (x+z)*lam + lam*lam) 
    l3 = math.sqrt(z*y + (z+y)*lam + lam*lam)
    lam_inc = l1 + l2 + l3
    return rf_duplicate_v2(x/4, y/4, z/4, (lam + lam_inc)/4, Q/4)


def rf_duplicate(x, y, z, Q = None):
    A = (x + y + z)/3
    X = 1 - (x/A)
    Y = 1 - (y/A)
    Z = -X - Y
    E2 = X*Y + X*Z + Y*Z
    E3 = X*Y*Z
    if not Q:
        r = numpy.finfo(float).eps / 4
        Q = (3*r)**(-1/8) * max(abs(A - x), abs(A - y), abs(A - z))
    if Q < A:
        return 1/math.sqrt(A) * (E2*E2*E3/16 - 5*E2*E2*E2/208 + 3*E3*E3/104 - 3*E2*E3/44 + E2*E2/24 + E3/14 - E2/10 + 1)

    lam = math.sqrt(x*y) + math.sqrt(y*z) + math.sqrt(z*x)
    return rf_duplicate((x + lam)/4, (y+lam)/4, (z+lam)/4, Q/4)

def ellipk(k):
    return rf_duplicate_v2(0, 1-k*k, 1)

# Maclaurin
def sncndn_maclaurin(x, k):
    if k < 0.5:
        k2 = k*k
        s = math.sin(x)
        c = math.cos(x)
        sn = s - k2/4*(x - s*c)*c
        cn = c + k2/4*(x-s*c)*s
        dn = 1 - k2/2*s*s
        return (sn, cn, dn)
    else:
        kp2 = 1 - k*k
        sh = math.sinh(x)
        ch = math.cosh(x)
        th = math.tanh(x)
        sn = th - kp2/4*(x/ch/ch - th)
        cn = 1/ch + kp2/4*(x/ch - sh)*th
        dn = 1/ch + kp2/4*(x/ch + sh)*th
        return (sn, cn, dn)

# AGM
def agm(a, b, c, x, n=0, epsilon=0):
    if epsilon == 0:
        epsilon = c * 1e-17
    an = 0.5*(a + b)
    bn = math.sqrt(a*b)
    cn = 0.5*(a-b)

    if (abs(cn - c) <= epsilon):
        return 2**n*an*x
    phi = agm(an, bn, cn, x, n+1, epsilon)
    phim1 = 0.5*(phi + math.asin(cn/an*math.sin(phi)))
    return phim1

def am_agm(x, k):
    kprime = math.sqrt(1 - k*k)
    return agm(1, kprime, math.sqrt(kprime), x)

def sncndn_agm(x, k):
    am = am_agm(x, k)
    sn = math.sin(am)
    cn = math.cos(am)
    dn = math.sqrt(1 - k*k*sn*sn)
    return (sn, cn, dn)

# Landen
def sncndn_maclaurin_inv(x, k):
    if k < 0.5:
        k2 = k*k
        s = math.sin(x)
        c = math.cos(x)
        sn = s - k2/4*(x - s*c)*c
        cn = c + k2/4*(x-s*c)*s
        dnm = k2/2*s*s
        return (sn, cn, dnm)

def sncndn_landen(x, k, initial = True):
    if k < 1e-3:
        return sncndn_maclaurin(x, k) if initial else sncndn_maclaurin_inv(x, k)

    kprime = math.sqrt(1-k*k)
    k1 = (1-kprime)/(1+kprime)
    x1 = x / (1 + k1)

    (sn1, cn1, dn1m) = sncndn_landen(x1, k1, False)
    sn = sn1*(1+k1) / (1 + k1*sn1*sn1)
    cn = cn1*(1 - dn1m) / (1 + k1*sn1*sn1)
    dnm = (4*dn1m - 2*dn1m*dn1m) / (k1 + (2*dn1m - dn1m*dn1m))
    return (sn, cn, 1 - dnm) if initial else (sn, cn, dnm)

def reduce_range(x):
    K = scipy.special.ellipk(k*k)
    xprime = math.remainder(x, 2*K)
    n = 0.5*(x - xprime)/K
    if abs(math.floor(math.remainder(n, 2) + 0.1)) > 0.5:
        return (xprime, -1, -1, 1)
    return (xprime, 1, 1, 1)

# Combined
def sncndn(x, k):
    # fix k
    if k > 1:
        (sn,cn,dn) = sncndn(x*k, 1.0/k)
        return (sn/k, dn, cn)
    if (k < 0):
        return sncndn(x, -k)
    # reduce x range
    x, sign_sn, sign_cn, sign_dn = reduce_range(x)
    #sign_sn, sign_cn, sign_dn = 1, 1, 1
    if (k > 0.999):
        (sn,cn,dn) = sncndn_maclaurin(x, k)
    if (k < 0.001):
        (sn,cn,dn) = sncndn_maclaurin(x, k)
    (sn,cn,dn) = sncndn_agm(x, k)
    return (sign_sn*sn, sign_cn*cn, sign_dn*dn)

# Complex argument
def sn(z, k):
    if k > 1:
        return sn(z*k, 1.0/k)/k
    if k < 0:
        return sn(z, -k)

    (snr, cnr, dnr) = sncndn(z.real, k)
    (sni, cni, dni) = sncndn(z.imag, math.sqrt(1 - k*k))

    d = cni*cni + k*k*sni*sni*snr*snr
    r = snr*dni
    i = sni*cni*cnr*dnr
    return r/d + 1j*i/d

# Inverse functions
def arcsn(u, k):
    s = u
    c = math.sqrt(1 - u**2)
    return s*rf_duplicate(c**2, 1 - k**2 * s**2, 1)

def arccn(u, k):
    s = math.sqrt(1 - u**2)
    c = u
    return s*rf_duplicate(c**2, 1 - k**2 * s**2, 1)

def arcdn(u, k):
    s = math.sqrt((1 - u*u) / (k*k))
    return arcsn(min(1, s), k)

# Utility
def ulp_error(approx, truth):
    epsilon = numpy.finfo(float).eps
    error = approx - truth
    ulp = approx*epsilon
    return math.log2(max(1, abs(float(error / ulp))))

mp.mp.prec = 256

# Test single values
k = 0.4
x = 3.6
print("x = ", x, "k = ", k)
truth_v = (mp.ellipfun('sn', x, k=k), mp.ellipfun('cn', x, k=k), mp.ellipfun('dn', x, k=k))
agm_v = sncndn_agm(x, k)
landen_v = sncndn_landen(x, k)
combine_v = sncndn(x, k)
print("combine: ", combine_v)
print("landen:  ", landen_v)
print("agm:     ", agm_v)
print("mpmath:  ", truth_v[0], truth_v[1], truth_v[2])
print("inverse: ", (arcsn(combine_v[0], k), arccn(combine_v[1], k), arcdn(combine_v[2], k)))

#exit()

# Calculate errors
max_period = 77.632480664198080934069420189 # Periodicity at 1-epsilon
#ks = [0, 0.00009, 0.2, 0.8, 0.999, 0.99991, 0.999999]
#xs = numpy.linspace(-max_period/2, max_period/2, 1000)

ks = [k]
xs = numpy.linspace(-5, 5, 10000)

if True:
    xs_plot = xs
    ys_plot = [sncndn(x, k) for x in xs_plot]
    plt.figure()
    plt.plot(xs_plot, [y[0] for y in ys_plot])
    plt.plot(xs_plot, [y[1] for y in ys_plot])
    plt.plot(xs_plot, [y[2] for y in ys_plot])
    plt.legend(["sn(u,k)", "cn(u,k)", "dn(u,k)"])

ulp_sn = []
ulp_cn = []
ulp_dn = []
for k in ks:
    print(ellipk(k))
    approx = [sncndn(x, k) for x in xs]
    truth = [(mp.ellipfun('sn', x, k=mp.mpmathify(k)), mp.ellipfun('cn', x, k=mp.mpmathify(k)), mp.ellipfun('dn', x, k=mp.mpmathify(k))) for x in xs]
    ulp_sn.append([ulp_error(v[0][0], v[1][0]) for v in zip(approx, truth)])
    ulp_cn.append([ulp_error(v[0][1], v[1][1]) for v in zip(approx, truth)])
    ulp_dn.append([ulp_error(v[0][2], v[1][2]) for v in zip(approx, truth)])

legend = ["k=" + str(k) for k in ks]

plt.figure()
plt.title("ULP error for sn(x,k)")
for p in zip(ks, ulp_sn):
    plt.plot(xs, numpy.abs(p[1]))
plt.legend(legend)

plt.figure()
plt.title("ULP error for cn(x,k)")
for p in zip(ks, ulp_cn):
    plt.plot(xs, numpy.abs(p[1]))
plt.legend(legend)

plt.figure()
plt.title("ULP error for dn(x,k)")
for p in zip(ks, ulp_dn):
    plt.plot(xs, numpy.abs(p[1]))
plt.legend(legend)

plt.show()