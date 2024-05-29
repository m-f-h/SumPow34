# SumPow34
In this GitHub repository we collect data and PARI/GP and Python code related to our paper
<!--M.F.Hasler &amp; G.Melfi (2024):--> 
"On sums of distinct powers of 3 and 4", to appear in "Combinatorics and Number Theory".

More precisely, the computational part of this research work consisted of the study of the function k(c) defined in Sect. 2, Def. 1, eq. (3), 
and the calculation of the numerical value of the exponent &gamma; in the asymptotic bound for the counting function,

> P<sub>{3,4}</sub>(x)  &Gt; x<sup>&gamma;</sup> ,

given in the main result, Theorem 1. As detailed especially towards the end of the proof of this theorem,
this requires the numerical integration of the function -log(g(t)) = -log(k(exp(t))) over t = 0 .. log(4/3),

Details and Python code concerning this numerical integration can be found in the folder [integral](integral). 

We can compute the function k(c), as well as its derivative, either exactly (in rational points the function has rational values) or numerically (floating point approximation) in given points.

As we show in the paper, the function k(c) has the property that it is piecewise affine, k(c) = a + b*c, or hyperbolic, k(c) = p + q/c.

We can determine the intervals J[i] = [c[i-1], c[i]] on which k(c) is of a given form, with fixed parameters (a, b) or (p, q).
These intervals have borders which are rational numbers,
except at points where the nature of the function changes from (piecewise) affine to hyperbolic or conversely.
Those points are of the form  c[i] = (m[i] +- sqrt(d[i]))/n[i] = r[i] +- sqrt(s[i]), where m[i], n[i], d[i] are integers and r[i], s[i] are rationals.

The folder [intervals](intervals) lists data we have obtained concerning these intervals, and the PARI/GP programs used to compute it.
(This data is provided on informative grounds, but it is in the current approach not used for the calculation of the asymptotic bound.)
