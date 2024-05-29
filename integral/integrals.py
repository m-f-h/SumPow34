""" integrals.py
    (c) 2023-2024 by M.F.Hasler & G.Melfi

Here we provide the functions used for numerical integration:

- The function `Romberg()` implements Romberg extrapolation,
  using the trapedoidal rule and Richardson extrapolation,
  with evaluation of the function at points t[k] = a + k*(b-a)/2^N
  for k = 0, ..., 2^N with increasingly larger N.

- The function `Riemann_sum()` can be useful to integrate a function
  for which one knows the exact value in given points (that are a priori
  not equidistant as the above) and its form between these points.
  
  This is in theory the case for the function k(c) we consider, 
  of which we can determine the exact (piecewise) form
  (piecewise affine or hyperbolic) and values within given subintervals,
  cf. data given in the folder [intervals](../intervals).
  However, in some regions of the integration domain there are 
  tens of thousands of extremely small subintervals (width ~ 10^-7 or smaller)
  on which k(c) has a given form.
  
  Also, we must integrate not just the function `k(c)`, but the more complicated
  function -log(k(exp(t))) in points  t = log(c),  where c are the (mostly rational) 
  points between which the function is affine or hyperbolic with given coefficients.
  
  In view of these complications it is in the case at hand not interesting to use
  this "exact" approach.
"""

def Romberg( f, a, b, Nmax=10, eps=1e-6, verbose=0 ):
    """Compute integral of the function f over the interval [a,b] 
    using Romberg extrapolation (with 2^N subintervals, N = 0, 1, ...),
    until the value changes by less than eps, or stop after at most Nmax iterations.
    """
    h = b - a	            # initial step: length of the interval
    S = [(f(a)+f(b))/2*h]	# starting value (Trapezoidal rule over [a,b])
    if verbose: print(end=f"Integral over ({a:.5g}, {b:.5g}): ")
    for n in range(1,Nmax):	# compute approx.integral for 2^n subintervals + extrapolation
        R = S			  # store previous value of the row S
        h /= 2			# step size divided by 2 : twice the number of subintervals
        if verbose > 1:
            print(f"Doing 2^{n-1} = {2**(n-1)} additional function evaluations.")
        # use the previously computed sum / 2 + sum over the new = odd-indexed points
        S = [R[0]/2 + h*sum(f(a+h*k) for k in range(1,2**n,2))]

        # Now use Richardson extrapolation; column m is estimator of order 2m
        for m in range(1,n+1):
            S.append((4**m*S[-1]-R[m-1])/(4**m-1))
        if verbose > 1:
            print(f"With a total of 2^{n}+1={2**n+1} evaluations, R =")
            print(S)	    	# Display the entire row/diagonal n° n
            print("Error =", S[-1] - R[-1])   # difference w.r.t. previous result
        if abs(R[-1] - S[-1]) < eps: break    # required precision is reached
    if verbose: 
        print(f"After {n} iterations, error is less than {R[-1] - S[-1]}.")
    return S[-1] 		# the last (thus most precise) approximation
  
""" Alternate integration method:
    For a given list of values (x[k], f(x[k])), compute the Riemann sum, using on each 
    subinterval [x[k-1], x[k]]  the average of the function values at the endpoints.
    This yields the exact integral if the function is affine on each of the subintervals.
    
    (Note: the implementation below must be modified in order to yield the exact 
    result for integration of the function  f(t) = -log(k(exp(t)))  if only the 
    function `k(c)` but not `f(t)` is piecewise affine in the given points.)
"""
    
def riemann_sum(f, x0 = None, f0 = None):
    """Compute Riemann sum (sum of (mean value = average of left & right border)
    times width of subintervals) for the map f = {x: f(x) ; x ∈ { x1, x2, ..., xN }}.
    (x0,f0) can be used to add a border value [possibly useful if the first point
    listed in f is not the same as the left border of the integration interval].
    """
    if x0 is not None: # if x0 is given but not f0, use f0 = f(min{x_i})
        f[x0] = f[min(f)] if f0 is None else f0
    x = sorted(f)
    return sum( (b-a)*(f[b]+f[a]) for a,b in zip(x,x[1:]) )/2
