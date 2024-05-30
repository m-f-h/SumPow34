# integral

In this subfolder we collect the Python programs used to compute
the exponent  &gamma;  in the asymptotic bound for the counting function

>  P<sub>{3,4}</sub>(x)  &Gt;  x<sup>&gamma;</sup> 

which is the result presented in Theorem 1 of the paper.

The main program [main.py](main.py) uses functions defined elsewhere in order to compute this exponent,
which requires performing the numerical integration of the function -log g(t) = -log k(e^t) 
over  t &isin; [0, log(4/3)].

The function `k(c)` as well as the points where it has its discontinuity and other "irregularities" are defined in [function_k.py](function_k.py).

The function `Romberg` used for the numerical integration is defined in [integrals.py](integrals.py).
