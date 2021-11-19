#!/usr/bin/env wolframscript

tau = -1.55 + 0.4*I
q = Exp[I*Pi*tau]
invTau = -1/tau
invQ = Exp[I*Pi*invTau]
z = 0.7 + 0.3*I

Print[EllipticTheta[1, z, invQ]]
Print[-I*Sqrt[tau/I]*Exp[I*tau*z^2/Pi]*EllipticTheta[1, tau*z, q]]

Print[Series[1, z, q]]