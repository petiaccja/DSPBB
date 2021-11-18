import numpy
import math
import cmath
import matplotlib.pyplot as plt

def reciprocal(tau):
    if tau.imag == 0 and tau.real == 0:
        return 1j
    return -1/tau

def shift(tau):
    (fracp, intp) = math.modf(tau.real)
    if fracp < -0.5:
        fracp = fracp + 1
    if fracp > 0.5:
        fracp = fracp - 1
    return complex(fracp, tau.imag)

def transform_tau_3_count(tau, max_iters=5):
    iter = 0
    if tau.imag < 0:
        return (-1, 0)
    while (tau.imag < 0.5 and iter < max_iters):
        tau = reciprocal(shift(tau))
        iter +=1
    return (tau, iter)

def transform_tau_3(tau, max_iters=5):
    (tau, iter) = transform_tau_3_count(tau, max_iters)
    return tau

# plot transform
resolution = 10001
real_span = 2
imag_span = 1
real_line = numpy.linspace(-real_span/2, real_span/2, resolution)
imag_line = numpy.linspace(imag_span/resolution, imag_span, resolution - 1)
plane = numpy.reshape(real_line, (1,resolution)) + 1j*numpy.reshape(imag_line, (resolution-1,1))

plane_transformed = numpy.zeros_like(plane)
worst_iter = 1
worst_tau = 0
for x in range(plane.shape[0]):
    for y in range(plane.shape[1]):
        (tau, iter) = transform_tau_3_count(plane[x,y], 15)
        if iter > worst_iter:
            worst_tau = plane[x,y]
            worst_iter = iter
        plane_transformed[x,y] = tau

print(worst_tau)
print(worst_iter)

plt.figure()
plt.imshow(numpy.clip(numpy.imag(plane_transformed), -2, 2), extent=(-real_span/2, real_span/2, 0, imag_span), origin='lower')
plt.show()