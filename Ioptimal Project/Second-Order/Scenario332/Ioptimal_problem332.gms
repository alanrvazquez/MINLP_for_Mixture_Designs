option limrow=0;
* 'number of cases output in the LST file for each named var.'
option limcol=0;
* # number of decimals being displayed in the output
option decimals=8;
* # lmit time
option reslim=14400;

* # Force a compilation error if GAMS Rev 240 or higher is not used
$version 240

* Definitions:
* ------------
* Number of parameters
$setglobal NPAR          6
$setglobal NPOINT        19
$setglobal NVAR          3
* Data sets:
* ----------
SET      np      'number of design points'      /1*%NPOINT%/;
SET      nth     'number of parameters'          /1*%NPAR%/;
SET      nv      'number of variables'           /1*%NVAR%/;
alias(np1,np);
alias(nth1,nth2, nth);

* Parameters:
*------------
PARAMETER        npz             'number of design points',
                 Ring(nv)        'Ingredient availability.'
                 lambda          'Amount of mixture for one experimental run'
                 BigM            'limit'
                 par(np)         'parameters - local values';
Ring('1')        =       10.2;
Ring('2')        =       4.0;
Ring('3')        =       4.9;
lambda           =       1;
npz              =       card(np);
BigM             =       npz;

* Variables:
*-----------
VARIABLES
x(nv,np)         'variable x - positions of the design points',
dif_th(nth,np)   'model matrix columns',
fim_th(nth,nth)  'FIM'
iv(nth,nth)      'FIM^-1 matrix'
lm(nth,nth)      'matrix L'
linv(nth,nth)    'matrix L^-1'
avpI             'average prediction variance'
z(nv,np)         'variable z - products of binary and continuous variables'

* Integer variable:
*------------------
BINARY VARIABLES
y(np) 'variables y - decision of including design point or not';

EQUATIONS
vec_fim1(np)           'parameter column',
vec_fim2(np)           'parameter column',
vec_fim3(np)           'parameter column',
vec_fim4(np)           'parameter column',
vec_fim5(np)           'parameter column',
vec_fim6(np)           'parameter column',
mat_fim(nth,nth)       'terms of the FIM after summing'

* eq for D-optimality (and E-optimality)
eq_l(nth,nth)          'terms of L matrix'
eq_u1(nth,nth)         'terms of L matrix =0'
symm_fim(nth,nth)      'Symmetry of the FIM'

* eq for A-optimality
eq_innod(nth,nth)      'terms non-diagonal'
eq_india(nth,nth)      'terms of diagonal'
eq_psd(nth,nth)        'PSD requirements'
*eq_lu(nth,nth)         'condition on the matrix L'
eq_invlu(nth,nth)      'condition on the inverse of L'
chol_sta(nth,nth)      'Cholesky factorization stability criteria'
chol_sta_inv(nth,nth)  'Cholesky factorization stability criteria L^-1'
symm_iv(nth,nth)       'Symmetry of L^-1'

* eq for mixtures
const_mixt(np)         'constraints on mixtures'
ing_avail(nv)         'ingredient availability constraints'

* eq for linearizing quadratic constraint (ingredient availability const.)
lin_iav_1(nv, np)        'linearization of ingredient availability const.'
lin_iav_2(nv, np)
lin_iav_3(nv, np)

* others
ord_nw(np)             'ordering points constraint'
est_model1         'model must be estimable'
est_model2         'model must be estimable'
est_model3         'model must be estimable'
est_model4         'model must be estimable'
est_model5         'model must be estimable'
est_model6         'model must be estimable'
ord_x(np)           'ordering of the proportions in first ingredient'
eq_avp                 'prediction var. computation - I-optimality'
;

* Equations:
*-----------
* derivatives wrt parameters
vec_fim1(np)  ..       dif_th('1',np) =e=   x('1',np);
vec_fim2(np)  ..       dif_th('2',np) =e=   x('2',np);
vec_fim3(np)  ..       dif_th('3',np) =e=   x('3',np);
vec_fim4(np)  ..       dif_th('4',np) =e=   x('1',np)*x('2',np);
vec_fim5(np)  ..       dif_th('5',np) =e=   x('1',np)*x('3',np);
vec_fim6(np)  ..       dif_th('6',np) =e=   x('2',np)*x('3',np);

* FIM
mat_fim(nth,nth1).. fim_th(nth, nth1) =e=
  sum[np,(y(np)/npz)*dif_th(nth, np)*dif_th(nth1, np)];
* Symmetry of FIM
symm_fim(nth,nth1)$(nth.ord lt nth1.ord).. fim_th(nth,nth1) =e= fim_th(nth1,nth);


* Cholesky L*L^T decomposition for A-, D- and E-optimality
eq_l(nth,nth1)$(nth.ord ge nth1.ord).. fim_th(nth,nth1)  =e=
  sum[nth2,lm(nth,nth2)*lm(nth1,nth2)];
eq_u1(nth,nth1)$(nth.ord lt nth1.ord).. lm(nth,nth1) =e= 0;
* Condition for numerical stability in Cholesky factorization
chol_sta(nth,nth).. fim_th(nth,nth) =g= sum[nth1$(nth1.ord lt nth.ord),sqr[lm(nth,nth1)]];


* Compute the inverse of L (L^-1), using Cholesky descomposition
eq_innod(nth,nth1)$(nth.ord ne nth1.ord).. sum[nth2,lm(nth,nth2)*linv(nth2,nth1)] =e= 0;
eq_india(nth,nth1)$(nth.ord eq nth1.ord).. sum[nth2,lm(nth,nth2)*linv(nth2,nth1)] =e= 1;
*eq_lu(nth,nth1)$(nth.ord lt nth1.ord).. lm(nth,nth1) =e= 0;
eq_invlu(nth,nth1)$(nth.ord lt nth1.ord).. linv(nth,nth1) =e= 0;


* equation for FIM^-1 for A-optimality
eq_psd(nth,nth1).. iv(nth,nth1) =e= sum[nth2,linv(nth2,nth)*linv(nth2,nth1)];
chol_sta_inv(nth,nth).. iv(nth,nth) =g= sum[nth1$(nth1.ord lt nth.ord),sqr[linv(nth,nth1)]];
symm_iv(nth,nth1)$(nth.ord lt nth1.ord).. iv(nth,nth1) =e= iv(nth1,nth);



* Mixture runs must sum to one
const_mixt(np) .. sum[nv, x(nv,np)] =e= 1;

* Ingredient availability constraints
ing_avail(nv) .. sum[np, lambda*z(nv,np)] =l= Ring(nv);
lin_iav_1(nv, np) .. z(nv,np) =l= y(np);
lin_iav_2(nv, np) .. z(nv,np) =l= x(nv,np);
lin_iav_3(nv, np) .. z(nv,np) =g= x(nv,np) - (1 - y(np));

* ordering constraint
ord_nw(np)$(np.ord gt 1).. y(np) =l= y(np-1);

* Model must be estimable
est_model1 .. y('1') =e= 1;
est_model2 .. y('2') =e= 1;
est_model3 .. y('3') =e= 1;
est_model4 .. y('4') =e= 1;
est_model5 .. y('5') =e= 1;
est_model6 .. y('6') =e= 1;

* ordering the components constraint
ord_x(np)$(np.ord gt 1).. x('1',np) =g= x('1',np-1);


* Input Moment Matrix.
parameters Mmat(nth,nth1);
$gdxin 'in.gdx'
$load Mmat
$gdxin
;

* equation for computation of I-optimality
eq_avp.. avpI =e= sum[nth, iv('1',nth)*Mmat(nth,'1')] + sum[nth, iv('2',nth)*Mmat(nth,'2')]
                  + sum[nth, iv('3',nth)*Mmat(nth,'3')] + sum[nth, iv('4',nth)*Mmat(nth,'4')]
                  + sum[nth, iv('5',nth)*Mmat(nth,'5')] + sum[nth, iv('6',nth)*Mmat(nth,'6')];


* problem
model COMMON /vec_fim1, vec_fim2, vec_fim3, vec_fim4, vec_fim5, vec_fim6,
              mat_fim, eq_l, eq_u1, chol_sta, ing_avail, ord_nw, symm_fim, ord_x,
              const_mixt, est_model1, est_model2, est_model3, est_model4,
              est_model5, est_model6, lin_iav_1, lin_iav_2, lin_iav_3/;
* model for I-optimal designs
model I_opt /COMMON,eq_innod,eq_india,eq_psd,eq_invlu,eq_avp,chol_sta_inv,
             symm_iv/;

* parameters for solving.-------------------------------------------------------
* Absolute gap: best_estimate - best_integer
option optcr = 1.0e-5;
* Relative gap
option optca = 1.0e-6;
option iterlim = 2000000000;
* maximum number of domain errors for a nonlinear solver.
option domlim = 1000;
option rminlp = conopt4;
option minlp = scip;
* basis acceptance threshold.
option bratio = 1;
option decimals = 8;
I_opt.optfile = 1;

* initial values.---------------------------------------------------------------
x.lo('1',np)            =       0.3;
x.lo('2',np)            =       0.0;
x.lo('3',np)            =       0.2;
x.up('1',np)            =       1.0;
x.up('2',np)            =       1.0;
x.up('3',np)            =       1.0;

* Product variables
z.lo('1',np)            =       0.0;
z.lo('2',np)            =       0.0;
z.lo('3',np)            =       0.0;
z.up('1',np)            =       1.0;
z.up('2',np)            =       1.0;
z.up('3',np)            =       1.0;

* Initial parameters and bounds
dif_th.lo(nth,np1)      =       -0.1;
dif_th.up(nth,np1)      =       1.1;

* initialization for the FIM.
fim_th.lo(nth, nth1)    =       -0.1;
fim_th.up(nth, nth1)    =       1.1;

* initialization for Cholesky descomposition, matrix L.
lm.lo(nth,nth1)         =       -1;
lm.lo(nth,nth)          =       1.0e-5;
lm.up(nth,nth1)         =       1;

* initialization for the inverse of the FIM.
iv.lo(nth,nth1)         =       -1000.0;
*iv.lo(nth,nth)          =       1e-5;
iv.up(nth,nth1)         =       1000.0;
iv.scale(nth,nth1)      =       1.0e0;

* Initialization of L^-1
linv.lo(nth,nth1) = -1000.0;
linv.lo(nth,nth) = 1e-5;
linv.up(nth,nth1) = 1000.0;

* initialization for the objective function.
avpI.lo = 0.5;
avpI.up = 2000000.0;

I_opt.scaleopt = 1;

* saving results
parameters
xsol(nv, np)         'save solutions for points for all 3 criteria',
nsol(np)             'save solutions for n for all 3 criteria',
opt                  'save optimal value'
modelstat            'save the status of the model'
time                 'save the CPU time'

* solution for I-optimality
solve I_opt using minlp minimizing avpI;
xsol(nv, np)        =       x.l(nv, np);
nsol(np)            =       y.l(np);
opt                 =       avpI.l;
modelstat           =       I_opt.modelstat;
time                =       I_opt.resusd;
