#!/usr/bin/env wolframscript

Print[JacobiCN[-9.97, SetPrecision[0.8, 100]^2]]

Print[CarlsonRF[1+2*I, 0.6+3*I, 1.55]]
Print[JacobiAmplitude[3.6, 0.4^2]]


(* Carlson symmetric form arguments for elliptic integral 1 / (S1*S2*S3) *)
Carlsonify[S1_, S2_, S3_, y_, x_] := Module[
    {
        X1 = S1[x],
        X2 = S2[x],
        X3 = S3[x],
        X4 = 1,

        Y1 = S1[y],
        Y2 = S2[y],
        Y3 = S3[y],
        Y4 = 1
    },
    U12 = (X1*X2*Y3*Y4 + Y1*Y2*X3*X4) / (x - y);
    U13 = (X1*X3*Y2*Y4 + Y1*Y3*X2*X4) / (x - y);
    U23 = (X2*X3*Y1*Y4 + Y2*Y3*X1*X4) / (x - y);
    {U12^2, U13^2, U23^2}
]

(* Inverse functions *)
s1SN[u_] := Sqrt[u]
s2SN[u_] := Sqrt[1 - u]
s3SN[u_] := Sqrt[1 - k*k*u]
carlsonArcSN = Carlsonify[s1SN, s2SN, s3SN, 0, x^2]

s1CN[u_] := Sqrt[u]
s2CN[u_] := Sqrt[1 - u]
s3CN[u_] := Sqrt[(1 - k*k) + k*k*u]
carlsonArcCN = Carlsonify[s1CN, s2CN, s3CN, x^2, 1]

s1DN[u_] := Sqrt[u]
s2DN[u_] := Sqrt[1 - u]
s3DN[u_] := Sqrt[u - (1 - k*k)]
carlsonArcDN = Carlsonify[s1DN, s2DN, s3DN, x^2, 1]

Print[carlsonArcSN]
Print[carlsonArcCN]
Print[carlsonArcDN]

(* Print[Integrate[1/Sqrt[(1-t^2)*(t^2-0.3^2)], {t, 0.5, 1}]] *)

snv = JacobiSN[0.7, 0.7^2]
cnv = JacobiCN[0.7, 0.7^2]
dnv = JacobiDN[0.7, 0.7^2]
srule = {x->snv, k->0.7}
crule = {x->cnv, k->0.7}
drule = {x->dnv, k->0.7}
Print[CarlsonRF[carlsonArcSN[[1]]/.srule, carlsonArcSN[[2]]/.srule, carlsonArcSN[[3]]/.srule]]
Print[CarlsonRF[carlsonArcCN[[1]]/.crule, carlsonArcCN[[2]]/.crule, carlsonArcCN[[3]]/.crule]]
Print[CarlsonRF[carlsonArcDN[[1]]/.drule, carlsonArcDN[[2]]/.drule, carlsonArcDN[[3]]/.drule]]