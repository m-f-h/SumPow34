""" function_k.py
    (c) 2023 - 2024 by M.F.Hasler & G.Melfi

Here we provide the definition of the function k used for the numerical integration,
and some data concerning this function, in particular the dictionary `points_of_interest`.

Specifically:
- Sc(c) returns the list of all (up to) 2^18 sums of distinct powers of 3 and of 4, the latter being multiplied by c.
- k_full(c) runs over all the sums a,b in S(c) where a and b are "nearest neighbors",
  and check for holes <=> a+d < b, and if such a hole is detected, 
  the size of the gap is added to a total "offset", which gives the value of b - Ac(b).
- k_min is a more efficient function to compute just k(c) (for numerical integration and plotting, for example),

  k(c) = min { (x - offset at x)/x ; x ∈ Sc(c) } = 1 - max { (offset at x)/x ; x ∈ Sc(c) }
"""

# "points of interest"
# This lists remarkable c-values where the function has special
# properties (minima, maxima, ...) or where transitions occur.
# ("hand made" from studies detailed elsewhere.)

if 'sqrt' not in vars(): sqrt = lambda x: x**.5
points_of_interest = {
        # some of the values we know exactly have been replaced by slightly rounded
        # up or down values in view of the use of this list in the numerical integration
        'start':    1,  #
        #---        #1.11565 # = 22313/20000 EXACT : 0.7992 : max
        'max1':     # = 1.1156916664123537, # k = 0.799218587
        #'c_157':   # = 1.1156963701360634, # k = 0.79921970974441
                    (sqrt(1317648910247977)-35239219)/950272, 
        'min1':     1.125,  # = 9/8 (EXACT) : k = 56107/73728 = 0.760999891493... = min
        'min1_end': 1.1264,  # 0.78929, # where the *average* slope becomes << 1
        'max2':     1.127125, #=9017/8000 : k = 0.79109779      (much zigzag here)
        #'min2a':   1.1277, # k = 0.79035
        #'max2b':   1.12794, # k = 0.791058
        #           1.128 = 141/125: k = 3654443/4620288 = 0.790955672... = max
        #'min2b':   1.12804,       # k = 0.7908079
        #'max2c':   1.128075,      # k = 0.79093694466
        'min2':     1.2,    # (EXACT) 0.76323... = min
        'max_ll':   1.20135498, # : 0.7675    # upper bound of LL, but rather a min of k !
        #'c_ll':    3**9/4**7, #(EXACT) = 1.20135498046875 : k = 533311/663552 ~ 0.8037215
        'min_ls':     1.201355, # k=0.8037215,    # lower bound of LS
        'min3_start': 1.26253,  # k=0.81, slope starts decreasing
        'min3':     1.2655,     # k=0.7940172
        'min3_end': 1.2695, # 0.813 , slope becomes flat again
        'end':      1.333333, # end
}

supercycle=1 # set this to 1 to use super-cycles, else normal (short or long) cycles

c_ll = 3**9/4**7 # = 1.20135498046875 (exactly).
if vars().get('verbose'):
    print("Upper limit of c-values for long+long supercycles: 3^9/4^7 =", c_ll)

def Sc(c, limit = None, append_3powR = True):
    """Return the sorted list of all (up to) 2^18 sums of distinct powers of 3
    (up to 3^R depending on the c-value) and distinct powers of 4,
    the latter being multiplied by c. If limit is given, use only terms <= limit.
    """
    if globals().get('supercycle'):
        R,L = (10,8) if c < c_ll else (9,7)
    else:
        R,L = (5,4) if c < 81/64 else (4,3)
    S = {0}
    if limit:
            if'log'not in vars(): from math import log
            if limit < 3**R: append_3powR = False
            elif limit > 3**R: limit = 3**R
            R = min(R, int(log(limit+0.9,3))+1) # we go only up to power R-1
            L = min(L, int(log(limit/c+0.5,4))+1) # we go only up to power L-1
            news = lambda x: {s+x for s in S if s <= limit - x}
    else:   news = lambda x: {s+x for s in S}
    for x in(  3**k for k in range(R)): S |= news(x)
    for x in(c*4**k for k in range(L)): S |= news(x)
    S = sorted(S)
    if append_3powR: S . append(3**R)
    if globals().get('verbose'): print([c, R, len(S)]) # some debugging info: R and #sums
    return S

def k_full(c, epsilon = 1e-8):
    """Return (k_min, x_min, k2, x2, g_max, x_max, x_max2) for a (super-)cycle 
    with the given c-value. The epsilon parameter is used to recognize 
    'gaps of equal size' even if there are small rounding errors. 
    NOTE: make sure that h >> epsilon when computing slope as (k(c+h)-k(c))/h !
    """
    a = o = t = G = o2 = t2 = x2 = xG2 = 0 ; d = 1/2 + c/3
    for b in Sc(c, append_3powR = False): # initial 0 could be ignored but does no harm
        if b > a+d:
            g = b - (a+d) ; o += g ; o2 += g
            if o/b > t: t,x = o/b, b
            if g >= G:
                rg = round(g/epsilon-1)*epsilon # slightly smaller to "capture" g == G
                if rg == G: xG2 = b                 # in case g has negative error
                else:
                    xG = b ; xG2 = o2 = 0 ; G = rg
            if o2 and o2/(b-xG) > t2: t2,x2 = o2/(b-xG), b
        a = b
    # it turns out that the x_max have tiny rounding errors (xxx.00000000012) which we round off
    return 1-t, x, 1-t2, x2, G+epsilon, round(xG,7), round(xG2,7), xG2-xG

# fonction simplifiée / optimisée qui calcule juste le k_min
def k_min(c):
    ''' # using conservative bounds
        #(i.e., yield the higher limit "in case of doubt")
    limit = (243 if c < 1.115696    else 19683 if c <= c_ll
        else 243 if c < 1.2433963   else    81 if c <= 1.2452
        else 243 if c < 1.2468      else    81 if c <= 1.262
        else 5200 if  c < 1.270     else    81 )
    '''
    limit = ((19683 if 1.115696 < c <= c_ll else 243) if c < 1.2433964
        else   5200 if 1.262    < c < 1.270
        else    243 if 1.2452   < c < 1.2468   else 81 )
    a = o = t = 0 ; d = 1/2+c/3
    for b in Sc(c, limit = limit, append_3powR = False):
        if b > a+d:
            o += b - a - d
            if o/b > t:
                t = o/b
            if limit and b >= limit: break
        a = b
    return 1-t
