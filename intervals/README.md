# intervals
In this directory we put .gp files (PARI/GP source files) with the data for the intervals we have computed with PARI/GP.

We are interested in (maximal) intervals  J[i] = [c[i-1], c[i]]  on which the function k(c) has a given form,

- either affine,  k(c)  =  k[i-1] + slope[i] * (c - c[i-1])  =  k[i] + slope[i] * (c - c[i])  [where the last parenthesis <= 0] 
- or hyperbolic,  k(c)  =  p[i] + q[i] / c
 
with fixed  slope  or  p, q  values.  (Obviously,  k[i] = k( c[i] ) = p[i] + q[i]/c[i].)

Each of the files named `"intervals xxx - yyy.gp"` essentially contains one long list (usually called `c_list_xxx`)
whose items are 3-element lists which code the intervals in the range from xxx to yyy.

Each interval is given as a list, say `c_list[i]`, with three elements, of the form:

  - either:   [ c_i /* = approx \*/, k_i, [ slope[i+1] ]     ], /* dist = (c_i - c_{i-1}) */
  - or:       [ c_i /* = approx \*/, k_i, [ p[i+1], q[i+1] ] ], /* dist = (c_i - c_{i-1}) */

Here, parts within /* ... \*/, are comments, ignored by PARI/GP, but useful for the human reader:
- on the one hand, "approx" is the decimal approximation of the exact value of c_i, given as fraction a/b or possibly a/b +- 1/b*sqrt(d).
- on the other hand, 'dist' lists the width/size/length of the *preceding* (!) interval, also as decimal approximation in scientific notation,
  m.mmm e-EE (denoting the value m.mmm * 10^-EE,  where  1 <= m.mmm < 10  and  EE is an integer, here typically between 4 and 9).

The nature of the function  k  over the interval J[i]  is thus obviously manifest from the length of the list which is the third and last element,
which has one single value for intervals on which  k(c)  is affine, and two values for intervals on which  k(c)  is hyperbolic.

There may be other comments in these files:
- Usually there is a header with the file name and summary of the data : how many data items (i.e., "intervals"),
  at which indices/c-values occur the transitions from affine to hyperbolic (and conversely), ...
- There may also be single-line comments of the form "\\ ..." preceding some of the data lines, with information about "incidents" during the computation.
  (For example, when heuristics suggested to back up one step and re-do the scan with a smaller step size, in order to avoid to miss the border of an interval.)
