""" main.py - (c) 2023 - 2024 by M. F. Hasler & G. Melfi

The purpose of this program is to compute the lower bound, 
given in Theorem 1 of the paper, on the counting function  
P{3,4}(x)  of the set  Sum(Pow({3,4}, 1)).

More precisely, we prove in the paper that  P{3,4}(x)  >>  x^gamma
with  gamma = 1 - (tau/2) (1/log 3 - 1/log 4),
where  tau = 1/log(4/3) integral_{0 < t < log(4/3)} -log(g(t)) dt.

In this "main" program we evaluate this integral of 
-log g(t) = -log k(c = e^t)  over  0 < t < log(4/3),
and the resulting constants  tau  and  gamma.

The function  g(t) = k(e^t)  is continuous except for a jump at the point c_o = 81/64.
However, there are some other points where the function has (inverse/downward) peaks,
cf. Figure 1 in the paper.
In order to minimize the impact of these "irregularities" on the precision of the numerical 
integration, we split up the integration over subintervals delimited by these points.
    
To keep things simple yet efficient, we use Romberg extrapolation,
with numerical evaluation of the function at points t[k] = a + k*(b-a)/2^N
for k = 0, ..., N increasingly larger until the required precision is reached.
"""
from function_k import k_min, points_of_interest
from integrals import Romberg

# In the final application, we will provide a more complete list
# of points where to split up the integral in several parts. 
# By default, we just list the only point of discontinuity.
# Also, the default precision 'eps' and 'Nmax' are rather low,
# for testing purpose.

c_o = 81/64  # point of discontinuity of the function k

def compute_exponent(k: callable, points: list = [1, c_o, 4/3], 
                     eps: float = 1e-5, Nmax: int = 9, verbose: int = 1):
    """Compute integral of  -log k(c = exp(t))  over  t in log(convex_hull(points)),
    using Romberg extrapolation on each of the subintervals delimited 
    by the given points (assumed to be given in increasing order).
    """
    if 'exp' not in vars(): from math import exp,log
    f = lambda t: -log(k(exp(t)))  # the integrand
    J = sum( Romberg(f, log(a), log(b), Nmax=Nmax, eps=eps, verbose=verbose)
             for a,b in zip(points, points[1:]))
    print("Integral =", J)
    # This and the following refer to definitions in the paper:
    print("tau =", t := J / log(4/3))  
    print("(1/ln 3 - 1/ln 4)*tau =", c := t/log(3) - t/log(4) )
    # There's a factor 1/2 if we use "supercycles" (which is what we do)
    print("Exponent 1 - ... =", gamma := 1 - (c/2 if supercycle else c) )
    return gamma

if __name__ == '__main__':
    # use up to 2^12 = 4096 function evaluations *on each subinterval*
    compute_exponent(k_min, sorted(points_of_interest.values()),
                     eps = 1e-7, Nmax = 12, verbose = 1) # fct. k_min

""" Output:
Integral over (0, 0.10947): After 9 iterations, error = 2.3444070431505315e-08.
Integral over (0.10947, 0.11778): After 9 iterations, error = -7.958040960779657e-08.
Integral over (0.11778, 0.11956): After 6 iterations, error = 5.55499209645041e-08.
Integral over (0.11956, 0.12045): After 2 iterations, error = 2.7679838492079226e-08.
Integral over (0.12045, 0.18232): After 11 iterations, error = -2.5973919757713015e-07.
Integral over (0.18232, 0.18345): After 3 iterations, error = 3.8671855372661695e-08.
Integral over (0.18345, 0.18345): After 1 iterations, error = 2.5463386827488014e-10.
Integral over (0.18345, 0.23312): After 9 iterations, error = -1.115865004563954e-10.
Integral over (0.23312, 0.23547): After 6 iterations, error = 1.5191516109305803e-08.
Integral over (0.23547, 0.23862): After 7 iterations, error = -8.999005272208871e-08.
Integral over (0.23862, 0.28766): After 7 iterations, error = 6.153258948551898e-08.
Integral = 0.06771071456280571
tau = 0.23536647238995012
(1/ln 3 - 1/ln 4)*tau = 0.04445877454785996
Exponent 1 - ... = 0.97777061272607
"""
