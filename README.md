# SumPow34
In this GitHub repository we collect data and PARI/GP and Python code related to our paper
<!--M.F.Hasler &amp; G.Melfi (2024):--> 
"On sums of distinct powers of 3 and 4", to appear in "Combinatorics and Number Theory".

More precisely, the computational part of this research work consisted of the study of the function k(c) defined in Sect. 2, Def. 1, eq. (3), and the numerical integration of log(g(t)) = log(k(exp(t))) over t = 0 .. log(4/3). 

Details and Python code concerning this numerical integration can be found in the folder [integral](integral). 

We can compute the function k(c), as well as its derivative, either exactly (in rational points the function has rational values) or numerically (floating point approximation) in given points.

As we show in the paper, the function k(c) has the property that it is piecewise affine, k(c) = a + b*c, or hyperbolic, k(c) = p + q/c.

We can determine the intervals J[i] = [c[i-1], c[i]] on which k(c) is of a given form, with fixed parameters (a, b) or (p, q).
These intervals have borders which are rational numbers, except at points where the nature of the function changes from (piecewise) affine to hyperbolic or conversely: 
These points are of the form  c[i] = (m[i] +- sqrt(d[i]))/n[i] = r[i] +- sqrt(s[i]), where m[i], n[i], d[i] are integers and/or r[i], s[i] are rationals.

The folder [intervals](intervals) lists the data we have collected concerning these intervals, and the PARI/GP programs used to compute this data.
