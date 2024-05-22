/* intervals 1.25 - 1.5.gp :
 *
 * (re)Computed 01-Dec-2023, MFH
 *
 * See ReadMe file for explanation of file format.
 *
 **********  LAST PART : c = 1.25 - 1.333     **********
 *
 *  We have a total of 194 points in [5/4 = 1.25 ,  3/2 = 1.5], including these endpoints.
 *
 *  Assume that there are 'o' earlier points ('o' remains to be determined...),
 *  so the n-th point in *this* list is c[o+n].
 *  For example,
 *  - the first point in the list below is c[o+1] = 5/4 = 1.25
 *  - The 174-th point is c[o+174] = 4/3 = 1.33333...
 *
 *  The global maximum of the function  k  is reached at the point just before,
 *  c[o+173] = 21/16 = 1.3125, where  k( c = 21/16 ) = 529/648 = 0.816(35802469)...
 *
 *  We have a transition from  p+q/c  to  m*c+k°  at the 169-th point :
 *  c[o+169] = 189123/8192 - 3/90112*sqrt(429430429945) ~ 1.26980595..., 
 *  k[o+169] = 924875/331776 - sqrt(429430429945)/331776 ~ 0.81249271971
 *  (This can be computed as  p + q/c  with parameters of the preceding interval J[o+168],
 *   or as  k'[o+169]*(c - c[o+170]) + k[o+170]  with the slope k' of the interval J[o+169] 
 *   and starting point (c, k(c)) of the next interval J[o+170] (= endpoint of J[o+169]).)
 */
c_list_125 = {[\\	(01:21) gp >  pq_list(10^-5,    5/4, )
[/* 1.25000000 =*/          5/4,               13/16,              [23/32, 15/128]],   \\dist = 0.e-19
[/* 1.26136364 =*/       111/88,            961/1184,              [25/96, 89/128]],   \\dist = 0.01136364
[/* 1.26253056 =*/    4131/3272,       107227/132192,     [-9007/1536, 34519/4096]],   \\dist = 0.001166926
[/* 1.26258005 =*/    1380/1093,     4583507/5652480,     [-2449/1536, 12439/4096]],   \\dist = 4.949255 e-5
\\Note: possible miss near c = 1.262660000 - double checked with step 1e-6:
\\(01:22) gp >  pq_list(10^-6,    1380/1093, )
[/* 1.26258005 =*/    1380/1093,     4583507/5652480,     [-2449/1536, 12439/4096]],   \\dist = 0.e-19
[/* 1.26264167 =*/    5793/4588,     1603165/1977344,      [-1581/512, 20163/4096]],   \\dist = 6.161904 e-5
[/* 1.26265244 =*/    8283/6560,       104189/128512,     [-8023/1536, 31207/4096]],   \\dist = 1.076509 e-5
[/* 1.26265535 =*/    8331/6598,     4610803/5687296,      [-1887/256, 42315/4096]],   \\dist = 2.911082 e-6
[/* 1.26267880 =*/      971/769,     3223803/3977216,      [-1467/512, 19011/4096]],   \\dist = 2.345353 e-5
[/* 1.26271552 =*/    5859/4640,       810409/999936,     [-6721/1536, 26823/4096]],   \\dist = 3.671360 e-5
[/* 1.26272727 =*/    1389/1100,     1152679/1422336,       [-121/1536, 4599/4096]],   \\dist = 1.175549 e-5
[/* 1.26274714 =*/    2427/1922,     4028083/4970496,        [-541/768, 7835/4096]],   \\dist = 1.986567 e-5
[/* 1.26278163 =*/    5829/4616,     2418469/2984448,       [-565/256, 15607/4096]],   \\dist = 3.449072 e-5
[/* 1.26279497 =*/    8241/6526,   13676321/16877568,     [-5523/1024, 32089/4096]],   \\dist = 1.334483 e-5
[/* 1.26279915 =*/    8337/6602,   13835287/17074176,    [-23171/3072, 43205/4096]],   \\dist = 4.177822 e-6
[/* 1.26286765 =*/      687/544,       189907/234496,     [-3587/3072, 10229/4096]],   \\dist = 6.849529 e-5
[/* 1.26294028 =*/    8247/6530,   13676359/16889856,     [-6691/1536, 26723/4096]],   \\dist = 7.262859 e-5
[/* 1.26305651 =*/    5901/4672,     2445035/3021312,      [-3009/512, 34591/4096]],   \\dist = 0.0001162312
[/* 1.26315789 =*/        24/19,         79501/98304,     [-2693/1024, 17791/4096]],   \\dist = 0.0001013879
\\Note: possible miss near c = 1.263261055 - double checked with step 1e-7:
\\(01:23) gp >  pq_list(10^-7,    24/19, )
\\[/* 1.26315789 =*/      24/19,         79501/98304,     [-2693/1024, 17791/4096]],   \\dist = 0.e-19
[/* 1.26325920 =*/    5907/4676,       305633/378048,    [-12755/3072, 25667/4096]],   \\dist = 0.0001013012
[/* 1.26326087 =*/    5811/4600,     4810615/5950464,    [-26555/3072, 48911/4096]],   \\dist = 1.673671 e-6
[/* 1.26326743 =*/     1214/961,   12059333/14917632,      [-2959/384, 44055/4096]],   \\dist = 6.560195 e-6
[/* 1.26333042 =*/    5781/4576,         37367/46248,      [-4103/384, 59471/4096]],   \\dist = 6.298982 e-5
[/* 1.26334951 =*/     1041/824,       430549/532992,      [-8927/768, 64329/4096]],   \\dist = 1.909498 e-5
[/* 1.26336375 =*/      969/767,     3205607/3969024,         [277/768, 2313/4096]],   \\dist = 1.424033 e-5
\\(01:27) gp >  pq_list(10^-6,    969/767, )
\\[/* 1.26336375 =*/    969/767,     3205607/3969024,         [277/768, 2313/4096]],   \\dist = 0.e-19
[/* 1.26343381 =*/      964/763,    9566905/11845632,      [1283/384, -13111/4096]],   \\dist = 7.005900 e-5
[/* 1.26346655 =*/    5817/4604,       687305/850944,      [-1063/256, 25669/4096]],   \\dist = 3.273693 e-5
[/* 1.26352182 =*/    8223/6508,     6799267/8420352,     [-8005/1536, 31151/4096]],   \\dist = 5.526847 e-5
[/* 1.26354454 =*/    1376/1089,   13651277/16908288,     [-9151/6144, 11887/4096]],   \\dist = 2.271697 e-5
[/* 1.26356707 =*/    8289/6560,   13705147/16975872,     [-5237/2048, 17413/4096]],   \\dist = 2.253690 e-5
[/* 1.26361689 =*/    8259/6536,   13653301/16914432,       [-8505/2048, 3209/512]],   \\dist = 4.981789 e-5
[/* 1.26380461 =*/    8331/6592,   13759757/17061888,    [-32107/6144, 15613/2048]],   \\dist = 0.0001877206
\\Note: possible miss near c = 1.263877755 - double checked with step 1e-7:
\\(01:27) gp >  pq_list(10^-7,   8331/6592, )
\\[ 1.26380461= */    8331/6592,   13759757/17061888,    [-32107/6144, 15613/2048]],   \\dist = 0.e-19
[/* 1.26387625 =*/    1389/1099,       382191/474112,    [-18919/6144, 10057/2048]],   \\dist = 7.163949 e-5
[/* 1.26387684 =*/    5829/4612,    9623267/11937792,    [-46591/6144, 21715/2048]],   \\dist = 5.918808 e-7
[/* 1.26394785 =*/    8337/6596,   13755751/17074176,    [-17729/2048, 12247/1024]],   \\dist = 7.100416 e-5
[/* 1.26395706 =*/    8241/6520,   13596191/16877568,    [-24249/2048, 32735/2048]],   \\dist = 9.208035 e-6
[/* 1.26397919 =*/      972/769,     1603187/1990656,      [-5793/2048, 9407/2048]],   \\dist = 2.213854 e-5
[/* 1.26400862 =*/     1173/928,        84109/104448,    [-30371/6144, 14881/2048]],   \\dist = 2.942693 e-5
[/* 1.26402944 =*/    1374/1087,       755243/937984,     [21805/6144, -7103/2048]],   \\dist = 2.081813 e-5
\\Note: possible miss near c = 1.264081412. 	\\(01:33) gp >  pq_list(10^-7/2,  1374/1087+10^-5*7/4 , )
\\[ 1.26404694 =*/549607609/434800000,2719040714245/3376789149696,     [21805/6144, -7103/2048]],      \\dist = 0.e-19
[/* 1.26408126 =*/    1369/1083,     3386699/4205568,    [34801/6144, -12579/2048]],   \\dist = 3.431695 e-5
[/* 1.26408146 =*/    5835/4616,    9623281/11950080,       [7105/6144, -909/2048]],   \\dist = 2.000349 e-7
[/* 1.26410178 =*/    8247/6524,   13601329/16889856,       [-8505/2048, 3209/512]],   \\dist = 2.032224 e-5
\\(01:35) gp >  pq_list(10^-6/2,  8247/6524 , )
[/* 1.26410178 =*/    8247/6524,   13601329/16889856,       [-8505/2048, 3209/512]],   \\dist = 0.e-19
[/* 1.26432292 =*/      971/768,     1599693/1988608,       [-5433/2048, 1119/256]],   \\dist = 0.0002211386
[/* 1.26439069 =*/    8259/6532,     4534439/5638144,    [-11965/2048, 17211/2048]],   \\dist = 6.777531 e-5
[/* 1.26439791 =*/      483/382,       265169/329728,      [-2379/256, 52293/4096]],   \\dist = 7.213781 e-6
[/* 1.26442516 =*/    5829/4610,     3199279/3979264,    [-16579/1536, 60065/4096]],   \\dist = 2.725693 e-5
[/* 1.26446281 =*/      153/121,       503633/626688,      [-3389/768, 27017/4096]],   \\dist = 3.764723 e-5
[/* 1.26453488 =*/      435/344,       178921/222720,        [-365/48, 43547/4096]],   \\dist = 7.207380 e-5
[/* 1.26458333 =*/      607/480,       187175/233088,          [-20/3, 38691/4096]],   \\dist = 4.844961 e-5
[/* 1.26464844 =*/    1295/1024,         12473/15540,          [-37/6, 36101/4096]],   \\dist = 6.510417 e-5
[/* 1.26466505 =*/    7287/5762,     3992343/4974592,    [-46531/6144, 10847/1024]],   \\dist = 1.660936 e-5
[/* 1.26470588 =*/        43/34,       211955/264192,     [-36739/6144, 8783/1024]],   \\dist = 4.083549 e-5
[/* 1.26477795 =*/    8259/6530,     4521171/5638144,    [-23267/3072, 43391/4096]],   \\dist = 7.206558 e-5
[/* 1.26477795 =*/    8259/6530,     4521171/5638144,    [-23267/3072, 43391/4096]],   \\dist = 7.206558 e-5
[/* 1.26531943 =*/    1941/1534,     3173399/3975168,   [-73889/3072, 128795/4096]],   \\dist = 0.0005414784
[/* 1.26536313 =*/      453/358,       246609/309248,   [-38287/1536, 133325/4096]],   \\dist = 4.370216 e-5
[/* 1.26538060 =*/    2427/1918,     3961943/4970496,    [-12209/384, 168921/4096]],   \\dist = 1.747631 e-5
[/* 1.26552323 =*/    5829/4606,    9471671/11937792,    [-24723/512, 254413/4096]],   \\dist = 0.0001426258
[/* 1.26553106 =*/     1263/998,     2051491/2586624,  [-175783/3072, 300723/4096]],   \\dist = 7.831555 e-6
[/* 1.26554868 =*/    7773/6142,   12612827/15919104,    [-52391/768, 357725/4096]],   \\dist = 1.761909 e-5
[/* 1.26555317 =*/    8259/6526,   13397291/16914432,   [-81819/1024, 418291/4096]],   \\dist = 4.490716 e-6
[/* 1.26556073 =*/    9231/7294,        82679/104448,  [-142787/1536, 485985/4096]],   \\dist = 7.562923 e-6
[/* 1.26556291 =*/    1911/1510,     3097399/3913728,    [-74791/768, 508917/4096]],   \\dist = 2.179057 e-6
[/* 1.26556686 =*/   10203/8062,   16530899/20895744,  [-343505/3072, 583739/4096]],   \\dist = 3.942952 e-6
[/* 1.26557550 =*/    2397/1894,     1293281/1636352,  [-349187/3072, 593327/4096]],   \\dist = 8.644725 e-6
[/* 1.26560454 =*/    5799/4582,     3118405/3958784,  [-362933/3072, 616523/4096]],   \\dist = 2.903792 e-5
[/* 1.26560612 =*/    6285/4966,     3379113/4290560,  [-185191/1536, 629093/4096]],   \\dist = 1.582125 e-6
\\Note: possible miss near c = 1.265610278. - double checked with step = 2.50000000 e-7.
[/* 1.26560968 =*/    7743/6118,     4161201/5285888,  [-379559/3072, 644579/4096]],   \\dist = 3.554738 e-6
[/* 1.26561058 =*/    8229/6502,       340145/432128,      [-6083/48, 661037/4096]],   \\dist = 9.049948 e-7
[/* 1.26561210 =*/    9201/7270,     4943229/6281216,  [-400217/3072, 679439/4096]],   \\dist = 1.523180 e-6
[/* 1.26561334 =*/   10173/8038,     5464549/6944768,  [-206137/1536, 699785/4096]],   \\dist = 1.232112 e-6
[/* 1.26562500 =*/        81/64,         32581/41472,    [69501/512, -699625/4096]],   \\dist = 1.166335 e-5
[/* 1.26563659 =*/   10239/8090,   16499831/20969472,  [134957/1024, -679147/4096]],   \\dist = 1.158838 e-5
[/* 1.26563780 =*/    9267/7322,   14935871/18978816,      [4103/32, -660613/4096]],   \\dist = 1.215498 e-6
[/* 1.26563930 =*/    8295/6554,   13371839/16988160,  [128019/1024, -644023/4096]],   \\dist = 1.500363 e-6
[/* 1.26564019 =*/    7809/6170,   12589787/15992832,    [62467/512, -628405/4096]],   \\dist = 8.902478 e-7
[/* 1.26564368 =*/    6351/5018,   10243523/13006848,  [122425/1024, -615703/4096]],   \\dist = 3.488253 e-6
[/* 1.26564523 =*/    5865/4634,    9461399/12011520,  [117791/1024, -592243/4096]],   \\dist = 1.548160 e-6
[/* 1.26567318 =*/    2463/1946,     3986027/5044224,  [115845/1024, -582391/4096]],   \\dist = 2.794484 e-5
[/* 1.26568311 =*/   10209/8066,   16540307/20908032,    [75793/768, -507525/4096]],   \\dist = 9.938562 e-6
[/* 1.26568502 =*/    1977/1562,     3203671/4048896,  [144557/1536, -483801/4096]],   \\dist = 1.904899 e-6
[/* 1.26568923 =*/    9237/7298,   14974163/18917376,  [248975/3072, -416063/4096]],   \\dist = 4.210720 e-6
[/* 1.26569678 =*/    1653/1306,        92459/116736,    [17755/256, -355453/4096]],   \\dist = 7.554147 e-6
[/* 1.26570127 =*/    7779/6146,   12622091/15931392,  [179257/3072, -298407/4096]],   \\dist = 4.485045 e-6
[/* 1.26571886 =*/    6321/4994,   10266719/12945408,   [75895/1536, -252053/4096]],   \\dist = 1.759352 e-5
[/* 1.26572668 =*/     1167/922,     1896187/2390016,    [12635/384, -166473/4096]],   \\dist = 7.818493 e-6
[/* 1.26585366 =*/      519/410,      846715/1062912,   [98005/3072, -161283/4096]],   \\dist = 0.0001269774
[/* 1.26586889 =*/    2433/1922,     3971147/4982784,   [25621/1024, -125599/4096]],   \\dist = 1.522804 e-5
[/* 1.26592978 =*/    1947/1538,     3182543/3987456,     [8703/1024, -39931/4096]],   \\dist = 6.089236 e-5
\\(01:40) gp >  pq_list(10^-6*2,  1947/1538, )
\\[ 1.26592978 =      1947/1538,     3182543/3987456,     [8703/1024, -39931/4096]],   \\dist = 0.e-19
[/* 1.26647257 =*/    8265/6526,   13565737/16926720,    [29377/6144, -10323/2048]],   \\dist = 0.0005427923
[/* 1.26648097 =*/    1364/1077,       216665/270336,    [42301/6144, -15779/2048]],   \\dist = 8.394392 e-6
[/* 1.26654412 =*/      689/544,     3394061/4233216,    [65149/6144, -25425/2048]],   \\dist = 6.315200 e-5
[/* 1.26655322 =*/    8187/6464,   13444421/16766976,    [17407/2048, -19967/2048]],   \\dist = 9.100175 e-6
[/* 1.26658562 =*/    7293/5758,   11979265/14936064,    [18913/3072, -27779/4096]],   \\dist = 3.240219 e-5
\\(01:41) gp >  pq_list(10^-6/2,  7293/5758, )
\\[ 1.26658562 =*/    7293/5758,   11979265/14936064,    [18913/3072, -27779/4096]],   \\dist = 0.e-19
[/* 1.26659642 =*/     1202/949,   11846891/14770176,         [85/12, -32587/4096]],   \\dist = 1.079727 e-5
[/* 1.26660156 =*/    1297/1024,           3121/3891,        [103/12, -40369/4096]],   \\dist = 5.145219 e-6
[/* 1.26661184 =*/    7701/6080,       131787/164288,        [243/32, -35235/4096]],   \\dist = 1.027961 e-5
[/* 1.26666667 =*/        19/15,         62451/77824,     [8901/1024, -40935/4096]],   \\dist = 5.482456 e-5
[/* 1.26671779 =*/    8259/6520,     2263103/2819072,     [5641/1024, -24417/4096]],   \\dist = 5.112474 e-5
[/* 1.26678933 =*/    1377/1087,       503261/626688,    [12163/1024, -57465/4096]],   \\dist = 7.153702 e-5
[/* 1.26682588 =*/    5835/4606,       640021/796672,     [7557/1024, -34125/4096]],   \\dist = 3.655086 e-5
[/* 1.26683938 =*/      489/386,       268207/333824,        [7973/2048, -501/128]],   \\dist = 1.349895 e-5
[/* 1.26684280 =*/      959/757,     1577995/1964032,      [14029/2048, -1961/256]],   \\dist = 3.422290 e-6
[/* 1.26686082 =*/    8265/6524,   13601173/16926720,      [7505/2048, -7423/2048]],   \\dist = 1.802105 e-5
[/* 1.26692708 =*/      973/768,     1601501/1992704,    [16721/2048, -19099/2048]],   \\dist = 6.626175 e-5
[/* 1.26694542 =*/    5757/4544,    9476941/11790336,    [31987/6144, -11423/2048]],   \\dist = 1.833920 e-5
[/* 1.26714636 =*/    8277/6532,     4545699/5650432,         [-673/6144, 593/512]],   \\dist = 0.0002009339
[/* 1.26715462 =*/    1385/1093,     6845683/8509440,         [12443/6144, -99/64]],   \\dist = 8.263912 e-6
[/* 1.26717391 =*/    5829/4600,    9603949/11937792,     [-15157/6144, 4245/1024]],   \\dist = 1.929273 e-5
[/* 1.26721763 =*/      460/363,       227339/282624,     [37115/6144, -6795/1024]],   \\dist = 4.371781 e-5
[/* 1.26722561 =*/    8313/6560,   13695265/17025024,       [23995/6144, -503/128]],   \\dist = 7.978902 e-6
[/* 1.26724138 =*/      147/116,        80729/100352,      [7573/2048, -7509/2048]],   \\dist = 1.576955 e-5
[/* 1.26727510 =*/      972/767,       533851/663552,    [25981/2048, -30837/2048]],   \\dist = 3.371847 e-5
[/* 1.26728886 =*/    8283/6536,     4549997/5654528,    [19445/2048, -11277/1024]],   \\dist = 1.376391 e-5
[/* 1.26733746 =*/    8187/6460,     4499125/5588992,    [51875/6144, -19825/2048]],   \\dist = 4.859961 e-5
[/* 1.26737619 =*/    5835/4604,       641505/796672,     [24251/6144, -8155/2048]],   \\dist = 3.873331 e-5
[/* 1.26740947 =*/      455/359,         32161/39936,    [37175/6144, -13615/2048]],   \\dist = 3.327614 e-5
[/* 1.26748144 =*/    8193/6464,     4505855/5593088,      [10237/2048, -2721/512]],   \\dist = 7.196489 e-5
[/* 1.26763804 =*/    1653/1304,      909675/1128448,      [7891/6144, -2483/4096]],   \\dist = 0.0001566012
[/* 1.26765799 =*/      341/269,     3377881/4190208,    [20803/6144, -13395/4096]],   \\dist = 1.995576 e-5
[/* 1.26770320 =*/    8235/6496,      906485/1124352,      [4769/2048, -7905/4096]],   \\dist = 4.520941 e-5
[/* 1.26770929 =*/    1378/1087,     4550629/5644288,     [8363/1024, -38221/4096]],   \\dist = 6.089658 e-6
[/* 1.26771533 =*/    8301/6548,     3426743/4250112,    [21815/3072, -32687/4096]],   \\dist = 6.041298 e-6
[/* 1.26772995 =*/    8187/6458,     4506649/5588992,     [5119/1024, -21771/4096]],   \\dist = 1.461443 e-5
[/* 1.26777971 =*/    5847/4612,       804805/997888,     [-7703/3072, 17209/4096]],   \\dist = 4.975776 e-5
[/* 1.26778784 =*/      980/773,    9711911/12042240,       [1573/3072, 1529/4096]],   \\dist = 8.134469 e-6
[/* 1.26788036 =*/      975/769,     1073567/1331200,    [38485/3072, -60871/4096]],   \\dist = 9.252452 e-5
[/* 1.26788793 =*/    5883/4640,     1619575/2008064,     [9735/1024, -45183/4096]],   \\dist = 7.566925 e-6
[/* 1.26790682 =*/    7293/5752,     2008067/2489344,    [26329/3072, -40321/4096]],   \\dist = 1.888399 e-5
[/* 1.26798094 =*/    5853/4616,     4837445/5993472,    [12481/3072, -16909/4096]],   \\dist = 7.412085 e-5
[/* 1.26798749 =*/     1216/959,   12060391/14942208,     [7679/1536, -21773/4096]],   \\dist = 6.551090 e-6
[/* 1.26806167 =*/    5757/4540,     1586549/1965056,      [1803/512, -14097/4096]],   \\dist = 7.418704 e-5
[/* 1.26816380 =*/      960/757,     1058537/1310720,      [3317/512, -29457/4096]],   \\dist = 0.0001021305
[/* 1.26826585 =*/    5763/4544,       794765/983552,     [7679/1536, -21773/4096]],   \\dist = 0.0001020406
[/* 1.26831137 =*/    8277/6526,     4566715/5650432,      [5569/3072, -5219/4096]],   \\dist = 4.552483 e-5
[/* 1.26838235 =*/      345/272,        95181/117760,    [25153/3072, -38339/4096]],   \\dist = 7.098304 e-5
[/* 1.26845329 =*/    8283/6530,     4572677/5654528,     [7679/1536, -21773/4096]],   \\dist = 7.093955 e-5
[/* 1.26847826 =*/     1167/920,       120809/149376,        [779/1536, 1567/4096]],   \\dist = 2.496837 e-5
[/* 1.26850886 =*/    2433/1918,     4029829/4982784,      [-1049/768, 11299/4096]],   \\dist = 3.060253 e-5
[/* 1.26851565 =*/    8187/6454,     4520035/5588992,      [-1775/512, 22215/4096]],   \\dist = 6.785810 e-6
[/* 1.26851852 =*/      137/108,       113455/140288,        [-263/512, 6871/4096]],   \\dist = 2.869309 e-6
[/* 1.26857517 =*/    5805/4576,     2403497/2972160,     [-3077/1536, 14611/4096]],   \\dist = 5.665631 e-5
[/* 1.26857888 =*/      973/767,    9668543/11956224,      [2107/384, -24309/4096]],   \\dist = 3.703923 e-6
[/* 1.26858254 =*/    5871/4628,     4861679/6011904,      [1019/256, -16481/4096]],   \\dist = 3.662306 e-6
[/* 1.26858736 =*/    1365/1076,     1130351/1397760,      [2095/256, -38321/4096]],   \\dist = 4.819540 e-6
[/* 1.26860841 =*/      392/309,     1298651/1605632,        [601/64, -44593/4096]],   \\dist = 2.105364 e-5
[/* 1.26862054 =*/    5757/4538,       501953/620544,      [4943/768, -29241/4096]],   \\dist = 1.212344 e-5
[/* 1.26862685 =*/    8241/6496,     1706575/2109696,      [3319/768, -18253/4096]],   \\dist = 6.309609 e-6
[/* 1.26865903 =*/    8193/6458,   13574575/16779264,      [3409/1536, -7329/4096]],   \\dist = 3.218027 e-5
[/* 1.26869806 =*/      458/361,     4553269/5627904,     [9907/1536, -29313/4096]],   \\dist = 3.903338 e-5
[/* 1.26875000 =*/      203/160,       252341/311808,    [12787/1536, -39057/4096]],   \\dist = 5.193906 e-5
[/* 1.26876923 =*/    8247/6500,     6835301/8444928,      [3179/512, -28061/4096]],   \\dist = 1.923077 e-5
[/* 1.26885776 =*/    2355/1856,      976393/1205760,     [7681/1536, -21781/4096]],   \\dist = 8.852785 e-5
[/* 1.26923077 =*/        33/26,         54811/67584,       [2437/2048, -985/2048]],   \\dist = 0.0003730106
[/* 1.26926648 =*/    1367/1077,     1135267/1399808,      [6745/2048, -6453/2048]],   \\dist = 3.571173 e-5
[/* 1.26930147 =*/    1381/1088,     2293981/2828288,      [14361/2048, -2015/256]],   \\dist = 3.498962 e-5
[/* 1.26933787 =*/    8205/6464,     2726465/3360768,     [30155/6144, -5325/1024]],   \\dist = 3.640070 e-5
[/* 1.26935146 =*/    2427/1912,     4032595/4970496,      [18683/6144, -1449/512]],   \\dist = 1.359315 e-5
[/* 1.26943118 =*/    5847/4606,    9716791/11974656,      [-8953/6144, 2949/1024]],   \\dist = 7.971229 e-5
[/* 1.26953125 =*/      325/256,     1619939/1996800,     [46343/6144, -8751/1024]],   \\dist = 0.0001000733
[/* 1.26959248 =*/      405/319,       224389/276480,    [69311/6144, -13611/1024]],   \\dist = 6.122649 e-5
[/* 1.26963124 =*/    5853/4610,    9732341/11986944,      [41651/6144, -3879/512]],   \\dist = 3.875995 e-5
[/* 1.26969603 =*/    8187/6448,   13618411/16766976,    [35203/6144, -12787/2048]],   \\dist = 6.479333 e-5
[/* 1.26971171 =*/    7311/5758,     4054055/4990976,      [9815/2048, -5175/1024]],   \\dist = 1.567568 e-5
[/* 1.26973974 =*/    5757/4534,     3192685/3930112,     [11309/6144, -1337/1024]],   \\dist = 2.803870 e-5
[/* 1.26976744 =*/      273/215,         21629/26624,     [24209/6144, -4067/1024]],   \\dist = 2.769771 e-5
[/* 1.26979167 =*/     1219/960,     6084851/7489536,      [29969/6144, -2643/512]],   \\dist = 2.422481 e-5
\\ the above is c[o+168]
\\ below is c[o+169] - transition from p+q/c  to  m*c+k°.
[/* 1.26980595 =*/189123/8192 - 3/90112*sqrt(429430429945),0.81249271971, [22/243]], \\dist = 1.428287 e-5
[/* 1.28225806 =*/      159/124,             227/279,                     [-34/81]],   \\dist = 0.01245211
[/* 1.28571429 =*/          9/7,             307/378,                      [50/81]],   \\dist = 0.003456221
[/* 1.28906250 =*/      165/128,             469/576,                     [22/243]],   \\dist = 0.003348214
[/* 1.31250000 =*/        21/16,             529/648,                    [-58/243]],   \\dist = 0.02343750
/*
 *
 *	*** NOTE : 	Here we reach c = 4/3, the end of the "official" domain.
 *				But we can extend the computation:
 *
 */ 
\\(01:48) gp >  pq_list(10^-4, 21/16, 5/3 )
\\[ 1.31250000 =*/        21/16,             529/648,                    [-58/243]],   \\dist = 0.e-19
[/* 1.33333333 =*/          4/3,           1183/1458,                     [32/243]],   \\dist = 0.02083333
[/* 1.34745763 =*/      159/118,           2591/3186,                     [-68/81]],   \\dist = 0.01412429
[/* 1.35000000 =*/        27/20,               73/90,                      [92/81]],   \\dist = 0.002542373
[/* 1.35245902 =*/      165/122,           2681/3294,                     [32/243]],   \\dist = 0.002459016
[/* 1.35294118 =*/        23/17,           6725/8262,                    [134/243]],   \\dist = 0.0004821601
[/* 1.35576923 =*/      141/104,           1145/1404,                      [10/81]],   \\dist = 0.002828054
[/* 1.35937500 =*/        87/64,             235/288,                     [-2/243]],   \\dist = 0.003605769
[/* 1.36363636 =*/        15/11,             727/891,                     [64/243]],   \\dist = 0.004261364
[/* 1.36764706 =*/        93/68,             125/153,                      [10/81]],   \\dist = 0.004010695
[/* 1.41000000 =*/      141/100,               37/45,                    [-70/243]],   \\dist = 0.04235294
[/* 1.41176471 =*/        24/17,           2263/2754,                    [134/243]],   \\dist = 0.001764706
[/* 1.41346154 =*/      147/104,             385/468,                      [10/81]],   \\dist = 0.001696833
[/* 1.41964286 =*/      159/112,             415/504,                    [-82/243]],   \\dist = 0.006181319
[/* 1.42105263 =*/        27/19,           2533/3078,                    [146/243]],   \\dist = 0.001409774
[/* 1.42241379 =*/      165/116,             215/261,                      [10/81]],   \\dist = 0.001361162
[/* 1.42741935 =*/      177/124,             230/279,                    [-94/243]],   \\dist = 0.005005562
[/* 1.42857143 =*/         10/7,           2803/3402,                    [158/243]],   \\dist = 0.001152074
[/* 1.42968750 =*/      183/128,             475/576,                      [10/81]],   \\dist = 0.001116071
[/* 1.43750000 =*/        23/16,             535/648,                      [14/27]],   \\dist = 0.007812500
[/* 1.43877551 =*/       141/98,            937/1134,                     [28/243]],   \\dist = 0.001275510
\\Note: possible miss near c = 1.500200000.
\\(01:53) gp > [28/243, 141/98, 937/1134].y(3/2) 
\\ = 5/6 : Using the last data to extrapolate the value at c = 3/2
[3/2, 5/6, []]];}
