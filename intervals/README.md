# intervals
In this directory we put .gp files (PARI/GP source files) with the data for the intervals we have computed with PARI/GP.

We are interested in (maximal) intervals  J[i] = [c[i-1], c[i]]  on which the function k(c) has a continuous derivative and a given form,

- either affine,  k(c)  =  k[i-1] + slope[i] * (c - c[i-1])
                        =   k[i] +  slope[i] * (c - c[i]),   <br/>
    where the last parenthesis is negative for c[i-1] < c < c[i]
- or hyperbolic,  k(c)  =  p[i] + q[i] / c
 
with fixed  slope  or  p, q  values.  (Obviously,  k[i]  =  k( c[i] )  =  p[i] + q[i]/c[i]  if hyperbolic.
So, the p-value is redundant when c, k(c) and q are known.)

Each of the files named `"intervals xxx - yyy.gp"` essentially contains one long list (formerly called `c_list_xxx`, now renamed to `interval_list_xxx`)
whose items are 3-element lists which code the intervals in the range from xxx to yyy.

Each interval is given as a list, say `interval_list[i]`, with three elements, of the form:

  - either:   [ c[i] /* = approx \*/, k[i], [ slope[i+1] ]     ], /* dist = c[i] - c[i-1] */
  - or:       [ c[i] /* = approx \*/, k[i], [ p[i+1], q[i+1] ] ], /* dist = c[i] - c[i-1] */

Here, parts within /* ... \*/, are comments, ignored by PARI/GP, but useful for the human reader:
- on the one hand, "approx" is the decimal approximation of the exact value of c_i, given as fraction a/b or possibly a/b +- 1/b*sqrt(d).
- on the other hand, 'dist' gives the width/size/length of the *preceding* (!) interval, also as decimal approximation in scientific notation,
  m.mmm e-EE (denoting the value m.mmm * 10^-EE,  where  1 <= m.mmm < 10  and  EE is an integer, here typically between 4 and 9).

The nature of the function  k  over the interval  J[i]  is thus obviously manifest from the length of the list which is the third and last element,
which has one single value for intervals on which  k(c)  is affine, and two values for intervals on which  k(c)  is hyperbolic.

There may be other comments in these files:
- Usually there is a header with the file name and summary of the data : how many data items (i.e., "intervals"),
  at which indices/c-values occur the transitions from affine to hyperbolic (and conversely), ...
- There may also be single-line comments of the form `\\ ...` preceding some of the data lines, with details about the computation.
  (For example, when the algorithm suggested to back up one or two steps and re-do the scan with a smaller
  step size, in order to avoid to miss a subinterval. 
  [This applies only to output produced by the older function `adaptive_scan(h, c0)`, formerly `pq_list()`; 
  the more recent function `intervals(a, b)` does not work the same way.])
- Sometimes, there is a multi-line comment /* ... */ between two data lines,
  when there is something noteworthy to say at that point.
  (Affine / hyperbolic transition, ...)
