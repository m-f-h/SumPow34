/***
    s34.gp - (c) 2023 - 2024 by M.F.Hasler

This file comprises several functions used for a detailed study
of the function
				k(c) = inf { A_c(x) / x ; 0 < x < max(A_c) },

where
				A_c(x) = integral_{t = 0..x} chi[A_c](t) dt

is the measure of the intersection of A_c with [0,x], here
written as integral over the characteristic function chi[A_c].

The set A_c is defined as

			A_c = Union_{t \in S_c} [t, t + 1/2 + c/3]

and finally,
			S_c = { alpha + c*beta ; alpha in S[3,R], beta in S[4,L] }

where
			S[3, R] = { sum_{ r in J } 3^r ;  J subset {0, ..., R-1} }

and analogous for S[4, L], and the limits (R, L) are (10,8) or (9,7)
or (5,4) or (4,3) depending on whether we use "super-cycles" or simple 
cycles, and whether  c < 3^9/4^7  or  c > 3^9/4^7  (for supercycles,
resp.   c < 81/64  or  c > 81/64  for simple cycles).

* The PARI/GP function  Sc(c,...)  defined below computes this 
  set  S_c  that may contain up to 2^(10 + 8 = 18) elements.

* The function  k_full(c)  returns  (k(c), x_min, count),
  where  x_min  is the x-value where the minimum was reached,
  and 'count' is the index of the corresponding element in  S_c.

* The function  k_normal(c)  returns  just the first component, k(c).

* The function  k_smart(c)  checks whether our database  ("c_list")
  contains information about the nature of  k  at  c, and the parameters
  (a,b)  or  (p,q)  such that  k(c) = a + b*c  or  k(c) = p + q/c.
  If so, this simple formula is used to compute  k(c), 
  otherwise the standard function k_normal(c).

  We have indeed proved in the paper that this function k(c) is piecewise 
  * either affine, i.e., of the form  k° + m*c 
  * or hyperbolic, i.e., of the form  p + q/c
  where  k°, m  resp.  p, q  are *rational* constants.

* The function "intervals(a, b)" (recursively) computes 
  the list of intervals on which  k  is of given form, 
  and the corresponding parameters  (k°, m)  or  (p, q).
  This computation is recursive and highly nontrivial.
  The current implementation may not always succeed, 
  and e.g. "get stuck" at some point.
  
  If it succeeds, the function returns a list  "c_list"  with items
  of the form  c_list[i] = [ c[i], k(c[i]),  (slope m)  or  p, q ]
  where c[1] = a  and  c[#c_list] = b.
 
* The function adaptive_scan(h, start)  (formerly pq_list(h, start))
  starts at  c = start,  with an initial step size 'stepsize',
  to look for intervals on which k(c) is of a given form.
  This function will rather print out its results in the 
  format detailed in the data files  "intervals xxx - yyy.gp".

* Several miscellaneous functions for custom symbolic calculations
  and other operations (intersection, extrapolation, ...) 
  related to (piecewise) affine or hyperbolic functions. 
  (See "arithmetics" and "curves" further below.)
  
***/

/********** obscure technical stuff **********
* Please ignore this initial part which is here for
* system-specific technical, non-mathematical reasons. 
*/
s34path = "C:/Users/MHasler/Dropbox/Recherche/Maths/NumberTheory/Melfi-Erdos"
\\ Append this to 'path', the list of directories where
\\ gp looks for files to include, unless already done.
if( #s34path > # PATH = default(path), default(path, s34path ";" PATH));
r34()={ read(); }\\ "reset" that may be required after s.th. went wrong
stats(reset=0)={
	if( reset>0, 	r_min = 99; cnt = 0; TOL = 1024.>>default(realbitprecision); print1("resetting "));
	iferr(printf("stats: r_min = %d, cnt = %d, TOL = %.3g. ", r_min, cnt, TOL), E, reset|| stats(1));
	if( reset == 0, print1("now "); stats(1));}
stats(1);
verbose=1;
/********** end of obscure technical stuff **********/

\\ here is the function k = k_normal (by default)
k = k_normal(c) = k_full(c)[1];

{addhelp(k_prime, "k_prime(c, k=10^-9):
	Compute  (k(c+h)-k(c))/h  which yields the exact slope if k(c)
	is affine on  [c, c+h], or else an approximate derivative.");
}
k_prime(c,h=10^-9) = (k(c+h)-k(c))/h;

{addhelp(k_full, "k_full(c):
	Compute  (k(c), x_min, count), where k(c) = min { A_c(x) / x ; 0 < x < max A_c }
	and  x_min  is the x-value where this minimum is reached,
	and 'count' counts the number of "gaps" or holes found in in A_c.
	Stops the search when  x >= current_limit  (a global variable).");
}
k_full(c)={ my( a=0, gap=0, sum_gaps = 0, max_ratio = 0, d = 1/2+c/3, cnt = 0);
    foreach ( Sc(c, x_limit(c)), b,
            gap = b - a - d;
            if ( gap > 0, cnt++;
                sum_gaps += gap;
                if ( sum_gaps / b > max_ratio,
                    x_min = b ; max_ratio = sum_gaps / b );
                if ( b >= current_limit, break );
            );
            a = b
    );
    [1-max_ratio, x_min, cnt];
}

\\ here is a smart version which uses the available list of intervals
\\ The list must be in the format [ J[i], i=1..N ] where J[i] = [x, kx, s, q]
\\ where x is the start of the range where this J[] is valid,
\\ kx = value of k(x), or 0 if this entry marks the end of a valid J[i-1]
\\ s = slope, if the function is linear on J[i], i.e., k(c) = kx + s*(c-x)
\\ q = q-value if the function is hyperbolic on J[i], i.e., k(c) = kx + q*(1/c - 1/x)
k_smart(c, list = 0)={ 
    my( ce=eval(c));
    if( list && (c = setsearch(list, [ce*1.,0,0,0], 1)) \\ index where [c,0,0,0] would be inserted
            && c <= #list && (c = list[c])[2],\\ our function is never zero. if needed this could be changed to [2..4]
            c[2] + if(c[4], c[4]*(1/ce-1/c[1]), c[3]*(ce-c[1]))
    ,   \\ fallback to normal k
        k_normal(ce)
    );}

supercycle=1;
c_ll = 3^9/4^7; \\ = 1.20135498046875 : upper limit of c-values for long+long supercycles
process_time = getabstime;

VERBOSE_Sc = 1<<8;
Sc(c, limit = 0, append_3powR = 0, timer = 0)={
    /* Return the sorted list of all (up to) 2^18 sums of distinct powers of 3
     * (up to 3^R depending on the c-value) and distinct powers of 4, the latter
     * being multiplied by c. If limit is given, use only terms <= limit.
     */
    my ( verbose = bitand(verbose, VERBOSE_Sc ), z=imag(c), c0=real(c));
    if ( timer, timer = process_time());
	\\ z && print("re="c0", im="z);
    [R,L] = if ( supercycle,    if ( c0 < c_ll, [10,8] , [9,7])
                ,               if ( c0 < 81/64, [5,4] , [4,3]));
    S = [0];
    if ( limit,
            if ( limit < 3^R, append_3powR = 0);
            R = min(R, logint(limit   ,3)+1); \\# we go only up to power R-1
            L = min(L, logint(max(limit\c0,1),4)+1); \\# we go only up to power L-1
            news = x-> [ s+x | s <- S, s+x <= limit ];
    ,/*else:*/
            news = x-> [ s+x | s <- S ]
        );
    for ( k=0, R-1, S = concat( S, news(3^k))); \\print(S);
    timer &&
        print1("done with 3^k. elapsed time = ",
               -timer + timer=process_time(), ". ");
    if ( z,
        if ( limit,
            news =  x-> vecsort([ s+x | s <- S , real(s+x) <= limit ], key, 8));
        key = x -> [real(x),imag(x)]
        );
    for ( k = 0, L-1, S = concat(S, news(c*4^k)));
    if ( timer,
        print1("done adding c*4^k. elapsed time = ",
               -timer + timer=process_time(), ". "));
    S = if ( z, vecsort(S, key, 8), vecsort(S,, 8));
    if ( timer,
        print1("done sorting. elapsed time =",
               -timer + timer=process_time(), ". \n"));
    \\#globals()['last_sorted'] = c0 - 1;
     
    if ( append_3powR, S = concat( S, 3^R));
    if ( verbose && limit && S[#S] != limit, print("limit = "limit", S[-1] = "S[#S]));
    current_limit = if ( limit , limit , S[#S] );
    verbose && print([c, R, #S]); \\# some debugging info: R and #sums
    S;
}

\\ cx=1.1156963701360625

{addhelp( x_limit, "x_limit(c = 0, z = 0):
    Return a "conservative" (i.e., in case of doubt, the larger one)
    limit adequate near c = 1 + z. (Either can be specified,
    but if c < 0 is given, we assume that z is meant with this.
    So if you really want c < 0, give a negative z < -1.)")
;}

x_limit(c = 0, z=0)={
    if ( !c,   c = 1 + z    \\# since c isn't, z *must* be defined
    ,   c < 1, c = 1 + c);  \\# given as 1st element, but it's z in reality
    if (c < 1.115696370136062, 
			243,\\ also in a tiny range 1.115 738 .. 1.115 755
		c < c_ll, \\ c_ll = 3^9/4^7 = 1.20135498046875
			19683, \\ceil(c<<14),	\\ c*4^7 <= 19683
		\\c < 1.2433963, 		\\ (3*sqrt(15889) - 219)/128
		\\	c*4^3, \\ ~ 80
		1.269806 <= c, 81 \\= x_min  for c > (2080353 - 3*√429430429945)/90112 ~ 1.2698059495
    ,   1.26253  <  c, 5300 \\ x_min = c*4^6  for c >= 4131/3272 ~ 1.26253  
    \\,   1.2452 < c && c < 1.2468,    243	????????????
    ,/*else*/ 243 \\ in most of the remaining region, c*4^3 <= 81 should suffice but never mind.
	)
;}



/*********** arithmetics with  x + y*sqrt_xxx_ *********
 *
 * Since PARI/GP would evaluate expressions like  1 + 2*sqrt(3)
 * (which occur as coordinates of the intersection of an affine
 *  function  a + b*x  and a hyperbolic function  p + q/x)
 * at once as a floating point number ("t_REAL"),
 * we store these as *polynomial* ("t_POL")  1 + 2*sqrt_3  
 * where sqrt_3 or sqrt_3_ is just an undefined variable name.
 *
 */

\\ inverse(a+b*sqrt(c)) = (a-b*sqrt(c))/(a^2-b^2*c)
inv(x)={if( poldegree(x)>0
			, 	my(v=variable(x)); subst(x,v,-v)/(polcoeff(x,0)^2-argument(v)*polcoeff(x,1)^2)
		,	1/x);}

argument(v)=eval(strsplit(Str(v), "_")[2]);

mul(a,b)=	if( poldegree(a *= b)>1, simplifySqrt(a), a);


addhelp(simplifySqrt, {"simplifySqrt(p):
	simplify polynomial in a variable of the form s_x_ by replacing s_x_^2 with x.");
}
simplifySqrt(p, v=variable(p))={ if( poldegree(p)>1, substpol(p, v^2, argument(v)), p);}	

/*
addhelp(sqrt2pol, "sqrt2pol(s): convert string involving 'Sqrt(xxx)' to polynomial in Sqrt_xxx_.");
sqrt2pol(s)={ (Sqrt(x)=eval(Str("Sqrt_"x"_"))); eval(Str(s));}
*/

\\return >0 if a>b, 0 if a=b, <0 if a<b
cmpf(a,b)= evalf(a - b);

\\evaluate after replacing variables using var2fct (=> t_REAL if it had sqrt_xxx_)
evalf(x)=eval(if( poldegree(x), var2value(x), x));

addhelp(var2fct, "var2fct(p): convert any pattern of the form 'f_x_' to string with 'f(x)'.");
var2fct(p, f=2)={ concat([ if( c=="_", ["(", ")"][ f=3-f ], c) | c<-Vec(Str(p)) ]);}

addhelp(var2value, {"var2value(x): replace variables in x by their value listed in SYMBOLS;
	if not yet listed there, determine the value using var2fct and store it.");
}
{var2value(x)=substvec(x, x=variables(x), [iferr(x = mapget(SYMBOLS, v), E,
	mapput(SYMBOLS, v, x = eval(var2fct(v))); x) | v<-x]);
}
SYMBOLS=Map();

val(x)=strprintf("%s = %.8f",x,evalf(x)*1.);


/*
{c_old=0;step=.5e-6;forstep(c=1.269 78, 1.269 81, step,
	c_old&& print(c=bestappr(c,1\/step),"  "interval_check(c_old,c)*1.);c_old=c)}
*/

/**************** main program ************/

{addhelp(interval_check, "interval_check(a, b, {ka=k(a)}, {kb=k(b)}, {m=(a+b)/2}, {km=k(ma)}): 
	Check whether the function k() appears to be mixed, linear or hyperbolic on the interval [a, b]
	by considering its value just at the middle point m = (a+b)/2.
    If linear, return [m, k(m), s] such that k(t) = k(a) + s*(t - a) on [a,b]
	if hyperbolic, return [m, k(m), 0, q] such that k(t) = k(a) + q/t - q/a on [a,b]
	otherwise return just [m, k(m)].");}

{interval_check(a, b, ka=k(a), kb=k(b),	m=(a+b)/2, km=iferr(k(m), E, k(m=bestappr(evalf(m)))))=
	\\ new style: a=[a,ka], b=[b,kb], m=[m,km] or omitted => we compute m=(a+b)/2 and km=k(m) 
	\\ if(	type(b) != "t_LIST",	[a, b, km] = [[a,m], [b,if(km,km,k(b))], k(m=(a+b)/2)]	);
	
    \\ m *must* be (a+b)/2 ! (simplifies many computations)
	\\ In case of an arbitrary "mid"-point m, we'd have to use:
    \\ if( (b-m)*(km-ka) == (m-a)*(km-ka), 1, (b-m)*(m*km-a*ka)==(m-a)*(b*kb-m*km), 2, 0)
	\\	b km  - b ka - m km + m ka = m km - m ka - a km + a ka
	\\  =
	if(	ka+kb==2*km, 		[m, km, (kb-ka)/(b-a)]
	,	a*ka+b*kb==2*m*km,	[m, km, 0, (ka-kb)*a*b/(b-a)]
	,	!poldegree(a+b),	[m,km]
	,	TOL>abs(ka+kb-2*km)
		, 	[m, km, if(	poldegree(a), (kb-km)/(b-m), (km-ka)/(m-a))]
	,	TOL>abs(evalf(a*ka+b*kb-2*m*km))
		,	[m, km, 0, 	\\simplifySqrt((ka-kb)*a*b*inv(b-a))
						\\bestappr((ka-kb)*evalf(a*b/(b-a)),1/TOL) 
						if( poldegree(a), m*b*(km-kb)/(b-m), m*a*(km-ka)/(a-m)) ]
\\	,	abs(	ka+kb - km*2)<TOL,  print("TOL-lin: ",[a,b]);		[m, km, (kb-ka)/(b-a)]
\\	,	abs(a*ka+b*kb - m*km*2)<TOL,  print("TOL-hyp: ",[a,b]);	[m, km, 0, a*b*(ka-kb)/(b-a)]
	,	[m, bestappr(km,1/TOL)]
	);}

\\ In order to be able to use "symbolic" (sqrt_xxx) expressions:
kbf(c)=iferr(k_normal(c), E, /*bestappr*/(k_normal(evalf(c))));

{addhelp(intervals, 
	"intervals(a, b, {r=30}, {ka = k(a)}, {kb = k(b)}, {m = (a+b)/2}, {km = k(m)}): 
	Make the list [J1, J2, ..., Jn] of intervals in [a,b] on which k(x) has 'constant slope'.
    'r' is the maximum recursion depth: if r <= 0 ==> error.
    Otherwise, do recursive checks to find the subintervals on which k() appears to be 
	linear or hyperbolic. The "intervals" are J[i] = [a[i], k[i], m[i], q[i]], where:
    a[i] = left border of the i-th interval
    k[i] = value of the function k() at a[i]
    m[i] = the slope, if k(x) = k[i] + m[i]*(x - a[i]) on the interval [a[i], b[i] = a[i+1]], else 0
    q[i] = the   q,   if k(x) = k[i] + q[i]*(1/x - 1/a[i])  -//-
	The last interval Jn in the list should be [b,kb] (or [b,kb,0,0]) meaning we don't know what comes after.
");}

{intervals(	a, b, r=30, ka=kbf(a), kb=kbf(b), m=(a+b)/2, \\ if a or b is irrational
			km=iferr(k(m), E, k(m=bestappr(evalf(m)))))=
	verbose>2&& printf("Intervals(%s): ",[val(a),val(b)]);
	cnt++ && r < r_min && if( 0 >= r_min=r, error("Deep recursion at a = "val(a)", b = "val(b)));

	\\ if(	type(b)!="t_LIST", [a,b,m]=[[a,m],[b,if(km,km,k(b)],[m=(a+b)/2,k(m)]);

    \\ First check whether we have a given type on either half,
	\\ computing two more values, L[1..2] and R[1..2].
    my( L = interval_check(a,m, ka,km), R = interval_check(m,b, km,kb));
	if( is_same(L,R), L[1]=a; L[2]=ka; return([L, [b,kb]]));

	\\ otherwise, if one side isn't OK, replace half of it with intervals(), then join the two
	join( 	if( 2 < #L,
					L[1]=a; L[2]=ka; [L, [m,km]],
				intervals(a, L[1], r-1, ka, L[2])
			)
		,	if( 2 < #R,
					R[1]=m; R[2]=km; [R, [b,kb]],
				intervals(R[1], b, r-1, R[2], kb);
			)
		,	r
	);
}

join(JL, JR, r)={
	verbose>1 && print("Joining "JL" and "JR);
	my( L, R);
	\\ compute intersection	of rightmost of left, and leftmost of right, if possible
	while( (L=evalf(JL[#JL][1])) != JR[1][1],\\ actually <, but... cmpf(JL[#JL][1], JR[1][1])<0,
		\\ to do : double check that we have > !
		verbose>5&& printf("Intersecting JL[%d]=%s and JR[-%d]=%s.\n", #JL-1, JL[#JL-1], #JR, JR[1]);
		my(x, ex, kx, JM);
		if( JL[#JL-1][3..-1] == JR[1][3..-1]
			,	\\if(	is_same( JL[#JL-1], JR[1] ) \\ better make the "refine" check.
				\\	,	verbose>5&& print("They are the same!");
				\\		return(concat(JL[^-1], JR[^1])));
				print("Same slope: can't intersect, must refine!");
		,	cmpf( JR[1][1], x = intersect( JL[#JL-1], JR[1] ))<0
			,	verbose>5&& print("Intersection "val(x)" to the right of JR[1]: refine.")
		,	cmpf(x, JL[#JL][1])<0
			,	verbose>5&& print("Intersection "val(x)" to the left of JL[-1]: refine.")
		,	x == JL[#JL][1]
			,	verbose>5&& print("*** Unexpected: Intersection "val(x)" = JL[-1].");
				JM = intervals(x, JR[1][1], r-1, JL[#JL][2], JR[1][2]);
		,	x == JR[1][1]
			,	verbose>5&& print("No interior intersection point, but x = JR[1].");
				JM = intervals(JL[#JL][1], x, r-1, JL[#JL][2], JR[1][2]);
		,	kx = iferr(k(x),E,k(evalf(x)));
			verbose&& print("found interior intersection point x = "val(x));
			JM = concat(intervals(JL[#JL][1], x, r-1, JL[#JL][2], kx)[^-1],
						intervals(x, JR[1][1], r-1, kx, JR[1][2]));
		);
		\\ do we need further subdivisions?
		#JM || JM = intervals( JL[#JL][1], JR[1][1], r-1, JL[#JL][2], JR[1][2] );
		if(#JL < #JR, JL = join(JL, JM), JR = join(JM, JR));
    );
	JL[#JL][1] == JR[1][1] || error("JL must end at the same point where JR starts.");
	while( #JL > 2 && is_same(JL[#JL-2], JL[#JL-1])
	, 	verbose&& print("Merging  JL[-2] !");
		JL = JL[^#JL-1] 
	);
	while( #JR > 2 && is_same(JR[1], JR[2])
	, 	verbose&& print("Merging  JR[2] !");
		JR = JR[^2]
	);
	while( #JL > 1 && is_same(JL[#JL-1], JR[1])
	, 	verbose>5&& print("Merging end of JL with start of JR.");
		JR = JR[^1]; JL[#JL] = JR[1]
	);
	verbose>2&& printf("OK, result: %s & %s (= %.6f, %.6f)\n",JL,JR,evalf(JL),evalf(JR));
    concat(JL[^-1], JR);
}


{addhelp(curves, "*** Curves and intervals ***
	Curves and intervals are coded as lists of the form J = [x0, y0, m, {q}, ...]
	which represents the function/curve y(x) = y0 + m(x - x0) + q/x - q/x0 + ...
	If the list is of length 2 only, it means we don't know of the functions form to the right of x0.
	If the list is of length 3, we have an affine function.
	If the list is of length 4 with nonzero q = J[4], we have a hyperbolic function.
	We can compute (extra-/interpolate) the function at any value  x  using  J.y(x).
	We can intersect curves using: intersect( J1, J2 ).
	We can check whether two curves are the same function using is_same( J1, J2 )
	");}
	
{addhelp(is_same, "is_same(curve1, curve2): 
	check whether J1 and J2 correspond to the same function.");}
is_same( J, K )={
	if(	#J>2 && #K>2 && J[3] && TOL>abs(J[3] - K[3])
		,	TOL>abs(evalf(K.y(J[1]) - J[2])) || print("yes but no.")
	,	#J>3 && #K>3 &&	J[4] && TOL>abs(J[4] - K[4])
		,	TOL>abs(evalf(K.y(J[1]) - J[2])) || print("Yes, but No.")
	);}

\\ Extrapolate a curve  C = [x0, y0, m, q] <=> y(x) = y0 + m*(x - x0) + q*(1/x - 1/x0)
C.y = if(#C>3 && C[4], x -> C[4]/x - C[4]/C[1] + C[2], x -> C[2] + C[3]*(x-C[1]));

{addhelp( intersect, "intersect(J, K): compute intersection point of curves J & K given in the form:
    J = [x0, y0, m, q] <=> y(x) = y0 + m(x - x0) + q/x - q/x0.
	(In current implementation, only one of m or q may be nonzero.)
Example: 
	gp > intersect([0,0,1,0],[4,0,0,1])
	%173 = \"-1/8 +1/8*sqrt(65)\"
	gp > [4,0,0,1].y(eval(%))
	%174 = 0.88278221853731870654582665378797139139
	gp > [0,0,1,0].y(eval(%173))
	%175 = 0.88278221853731870654582665378797139139
");}

intersect( J, K, D = if(#J==#K || J=Vec(J,#K=Vec(K,max(#J,#K))), J-K) )={
    if( #D>3 && J[4] && K[4]					\\  y(x) = y0 + q0*(1/x - 1/x0) = y1 + q1*(1/x - 1/x1)
		,	D[4] / ( J[4]/J[1] - K[4]/K[1] - D[2] )	\\	=> y0 - y1 + q1/x1 - q0/x0 = (q1 - q0)/x
	,	#D<4 || !D[4]	\\  both q's are 0:	y(x) = y0 + m0*(x - x0) = y1 + m1*(x - x1)
		,	( J[3]*J[1] - K[3]*K[1] - D[2] ) / D[3]	\\ => y0 - y1 + m1*x1 - m0*x0 = x*(m1 - m0)
    ,   \\ exactly one has q != 0: mixed: swap so that J is the affine curve
        #J>3 && J[4] && [J,K,D] = [K,J,-D];
        \\ m(x-x0)+y0 = y1 + q1/x - q/x1  <=>  m*x² - x*(m*x0 - y0 + y1 - q/x1) = q
        \\ <=>  x² - x (x0+(y1-q/x1-y0)/m) = q/m :  -b = x0+(y1-q/x1-y0)/m ; -c = q/m.
        my( b = (J[1] - ( D[2]+K[4]/K[1] )/J[3])/2  \\ This is -b,
        ,    s = iferr( make_sqrt( b^2 + K[4]/J[3], 0 ),  		\\ sqrt(b²/4 - c)
						E, make_sqrt( bestappr( b^2 + K[4]/J[3], 0 )))
		);
        \\ The intersection points are x_{1/2} = -b/2 +- sqrt(b²/4 - c).
		\\ We are interested in the smallest positive one.
        \\ iferr( if( b>s, b-s, b+s ), E, Str( b, if( b > eval(s)," -", " +" ), s))
		if( b>0 && K[4]/J[3] < 0, b-s, b+s) 
	);}

{addhelp( make_sqrt, "make_sqrt(x, {str=1}):
	Return an expression the form 'f*sqrt(s)' equal to sqrt(x), where 'f' is a rational and 
	's' is a square-free integer, either as str, or as polynomial in a variable 'sqrt_a_'.
	If |f| = 1, then 'f*' is omitted (for f = 1) or replaced by '-' (for f = -1).");}

make_sqrt(x, str=1)={
    my(f = factor(x), s = factorback([t[1] | t <- f~, t[2]%2]));\\ factors to an odd power
    f[,2] \= 2; f = factorback(f);
	if(	s==1,	f
	,	str	\\ distinguish special case f = 1
		,	Str(if(abs(f) != 1, Str(f"*"), f==-1,"-", ""),"sqrt("s")")
	,	s>0 || ((f*=I) && s*=-1)
		,	f*eval(Str("sqrt_"s"_"))
	);}


{addhelp( polroot, "polroot(P, {v}):
    Find the root(s) of the polynomial P in the variable v. 
    For degree > 2 (or < 1) we return the result of PARI/GP's polroots(), otherwise the exact solution(s):
    For degree 1, P = ax + b, =>  -b/a
    For degree 2, P = ax² + bx + c => [-b/2a +- sqrt(b²-4ac)/2a] (list of 2 values, unless sqrt = 0).
    Irrational roots are returned as strings \"A +-B*sqrt(C)\" where 
    coefficients A, B are rational and C is an integer, cf ?make_sqrt.")}

polroot(P, v=[])={ P=Vec(if(v, Pol(P,v), P));
    if( #P==2, -P[2]/P[1],
        #P==3, P /= P[1]; my( mb2 = -P[2]/2, sqrt = make_sqrt( mb2^2 - P[3] ));
            if( !sqrt,  mb2, 
                iferr([mb2-r, mb2+r], E,\\ error <=> make_sqrt has returned a string
                    [Str(mb2, sign, sqrt) | sign <- [" -", " +"]])),
        polroots(Pol(P)) \\ numerical
        \\error("Degree = ",#P-1," not yet implemented.");
    );}

