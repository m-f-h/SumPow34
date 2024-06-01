/* 	intervals 1.201355 - 1.243396
*	  (c) Nov. 2023 - 31.5.2024 by M.F.H.
*	
*	  This list actually starts at:
*		      c_LL = 3^9/4^7 = 19683/16384 = 1.20135498046875 (exactly)
*	  where k(c) = 533311/663552 = 0.8037214867862654320987654320987654321,
*
*   and ends at (notation of the paper, Table 1):
*		      c(5) = (3*sqrt(15889)-219)/128 ~ 1.243396329970176203817682739865128173
*   where k(c) = (sqrt(15889) + 343)/576 = 0.814325654071864903986495017767787272
*
*   The table below was produced from c_list_1201355 given towards the end of the file, using:
*/
\\(02:12) gp > make_intervals([x[1] | x <- c_list_1201355]) \\ minor editing by hand below:
c_list_1201355={[/***  List of "intervals" from c = 19683/16384 = 1.20135498046875 to (3*sqrt(15889)-219)/128 = 1.24339632997 ***/
[/*1.20135498 =*/  19683/16384,       533311/663552,           [446/243]],      \\dist = 0.e-19
[/*1.20165746 =*/      435/362,           7861/9774,           [614/729]],      \\dist = 0.0003024781
[/*1.20312500 =*/        77/64,         18791/23328,           [998/729]],      \\dist = 0.001467541
[/*1.20466321 =*/      465/386,          8417/10422,             [68/81]],      \\dist = 0.001538212
[/*1.20491803 =*/      147/122,            887/1098,           [124/729]],      \\dist = 0.0002548204
[/*1.20652174 =*/       111/92,           3011/3726,          [-112/243]],      \\dist = 0.001603706
[/*1.20779221 =*/        93/77,           1439/1782,           [196/243]],      \\dist = 0.001270469
[/*1.20905172 =*/      561/464,           7597/9396,           [124/729]],      \\dist = 0.001259516
[/*1.21601942 =*/      501/412,         13511/16686,         [-1112/729]],      \\dist = 0.006967693
[/*1.21727749 =*/      465/382,         24995/30942,         [-2258/729]],      \\dist = 0.001258069
[/*1.21739130 =*/        28/23,         27077/33534,           [226/729]],      \\dist = 0.0001138174
[/*1.21875000 =*/        39/32,             349/432,           [406/243]],      \\dist = 0.001358696
[/*1.22020725 =*/      471/386,           2815/3474,            [20/243]],      \\dist = 0.001457254
[/*1.22246696 =*/      555/454,          9935/12258,          [-848/729]],      \\dist = 0.002259706
[/*1.22368421 =*/        93/76,           7471/9234,           [976/729]],      \\dist = 0.001217250
[/*1.22489083 =*/      561/458,         10025/12366,            [20/243]],      \\dist = 0.001206619
[/*1.22802198 =*/      447/364,           3985/4914,          [-304/729]],      \\dist = 0.003131148
[/*1.22950820 =*/        75/61,         24023/29646,           [428/729]],      \\dist = 0.001486219
[/*1.23000000 =*/      123/100,           3283/4050,            [76/243]],      \\dist = 0.0004918033
[/*1.23076923 =*/        16/13,           5123/6318,           [206/243]],      \\dist = 0.0007692308
[/*1.23097826 =*/      453/368,         36263/44712,           [250/729]],      \\dist = 0.0002090301
[/*1.23214286 =*/        69/56,           5521/6804,          [-142/729]],      \\dist = 0.001164596
[/*1.23529412 =*/        21/17,           2233/2754,           [266/729]],      \\dist = 0.003151261
[/*1.23670213 =*/      465/376,         12355/15228,          [-538/243]],      \\dist = 0.001408010
[/*1.23809524 =*/        26/21,          8249/10206,           [806/243]],      \\dist = 0.001393110
[/*1.23883929 =*/      555/448,         44129/54432,          [1970/729]],      \\dist = 0.0007440476
[/*1.23947368 =*/      471/380,           3751/4617,            [70/729]],      \\dist = 0.0006343985
[/*1.24000000 =*/        31/25,           5923/7290,           [970/729]],      \\dist = 0.0005263158
[/*1.24038462 =*/      129/104,         10273/12636,           [254/243]],      \\dist = 0.0003846154
[/*1.24115044 =*/      561/452,         22346/27459,           [310/729]],      \\dist = 0.0007658271
[/*1.24218750 =*/      159/128,             469/576,              [2/27]]       \\dist = 0.001037058
/*
 * last point (not in this list but next point in subsequent list)
 * is c(5) in the notation of our paper, Table 1:
 *
 * 1.24339632997 =* /   c5 = (3*sqrt(15889)-219)/128,
 * 					    k5 = (sqrt(15889) + 343)/576 /* = 0.81432565407186490398649501776778727208 * /
 * (k(c(5)) can be computed as: k5 = 469/576 + (c5 - 159/128) * 2/27 	* /		, []],
 */
]} /*** end of clist_1201355 ***/
/*
*	Below how the above list can be re-computed using the intervals() function,
* which returns the c_list is in the format [ [c[i], k(c[i]), slope[i+1] ; i=1..#c_list ]
*
(00:49) gp > intervals(bestappr(3^9/4^7),bestappr(1.24339632997)) \\ could also be 1.2433963299701762038176827398651281732
found interior intersection point x = 435/362 = 1.20165746
found interior intersection point x = 77/64 = 1.20312500
found interior intersection point x = 1053/874 = 1.20480549
found interior intersection point x = 465/386 = 1.20466321
Same slope: can't intersect, must refine!
found interior intersection point x = 147/122 = 1.20491803
Same slope: can't intersect, must refine!
found interior intersection point x = 561/464 = 1.20905172
Same slope: can't intersect, must refine!
found interior intersection point x = 93/77 = 1.20779221
found interior intersection point x = 111/92 = 1.20652174
found interior intersection point x = 293/236 = 1.24152542
found interior intersection point x = 561/452 = 1.24115044
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
found interior intersection point x = 159/128 = 1.24218750
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
found interior intersection point x = 1119/904 = 1.23783186
found interior intersection point x = 21/17 = 1.23529412
Same slope: can't intersect, must refine!
found interior intersection point x = 465/376 = 1.23670213
found interior intersection point x = 26/21 = 1.23809524
found interior intersection point x = 129/104 = 1.24038462
Same slope: can't intersect, must refine!
found interior intersection point x = 31/25 = 1.24000000
found interior intersection point x = 1455/1174 = 1.23935264
found interior intersection point x = 555/448 = 1.23883929
found interior intersection point x = 471/380 = 1.23947368
found interior intersection point x = 327/266 = 1.22932331
found interior intersection point x = 501/412 = 1.21601942
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
found interior intersection point x = 447/364 = 1.22802198
found interior intersection point x = 543/446 = 1.21748879
found interior intersection point x = 465/382 = 1.21727749
found interior intersection point x = 28/23 = 1.21739130
Merging  JL[-2] !
found interior intersection point x = 561/458 = 1.22489083
Same slope: can't intersect, must refine!
found interior intersection point x = 93/76 = 1.22368421
found interior intersection point x = 39/32 = 1.21875000
found interior intersection point x = 555/454 = 1.22246696
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
Same slope: can't intersect, must refine!
found interior intersection point x = 471/386 = 1.22020725
found interior intersection point x = 75/61 = 1.22950820
Same slope: can't intersect, must refine!
found interior intersection point x = 123/100 = 1.23000000
Same slope: can't intersect, must refine!
found interior intersection point x = 69/56 = 1.23214286
found interior intersection point x = 16/13 = 1.23076923
Same slope: can't intersect, must refine!
found interior intersection point x = 453/368 = 1.23097826
time = 657 ms.
%1202 =
c_list_1201355 = [{
[19683/16384, 533311/663552, 446/243], 
[435/362,   7861/9774,   614/729], [  77/64, 18791/23328,   998/729], 
[465/386,  8417/10422,     68/81], [147/122,    887/1098,   124/729], 
[ 111/92,   3011/3726,  -112/243], [  93/77,   1439/1782,   196/243], 
[561/464,   7597/9396,   124/729], [501/412, 13511/16686, -1112/729], 
[465/382, 24995/30942, -2258/729], [  28/23, 27077/33534,   226/729], 
[  39/32,     349/432,   406/243], [471/386,   2815/3474,    20/243], 
[555/454,  9935/12258,  -848/729], [  93/76,   7471/9234,   976/729], 
[561/458, 10025/12366,    20/243], 
[447/364,   3985/4914,  -304/729], [ 75/61, 24023/29646,  428/729], 
[123/100,   3283/4050,    76/243], [ 16/13,   5123/6318,  206/243], 
[453/368, 36263/44712,   250/729], [ 69/56,   5521/6804, -142/729], [21/17, 2233/2754, 266/729], 

[465/376, 12355/15228,  -538/243], [  26/21, 8249/10206,  806/243], 
[555/448, 44129/54432,  1970/729], [471/380,  3751/4617,   70/729], [31/25, 5923/7290, 970/729], 
[129/104, 10273/12636,   254/243], 
[561/452, 22346/27459,   310/729],
[159/128,     469/576,      2/27],
\\ and we insert by hand the exact value of the right border of the interval list
\\ (excluded from the above table):
[  c5 = (3*sqrt(15889)-219)/128           /* = 1.24339632997...* /,
   k5 \\ = 469/576 + (c5 - 159/128) * 2/27 
      = (sqrt(15889) + 343)/576           /* = 0.81432565407186490398649501776778727208 * /, 
   0] \\ right border of interval list.
]}
*/
/*** end of file ***/
