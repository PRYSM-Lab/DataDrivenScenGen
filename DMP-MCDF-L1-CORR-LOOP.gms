*_______________________________________________________________________________
*
*                               DMP MODEL C
*
*              MILP formulation integrating with k-means clustering
*
*
*_______________________________________________________________________________

SETS
i       set of uncertain parameters             /i1*i2/
n       set of uncertain scenarios              /n1*n1000/
m       set of 4 statistical moments            /m1*m4/
cl      set of prespecified nodes - clusters    /c1*c20/
mw      set of models regarding CDF weight      /mw1*mw3/
cln(cl,n)   mapping matrix of data points to clusters using k-means
;

ALIAS (i,q);
ALIAS (n,nn);
ALIAS (m,k);
alias (cl,cl1,cl2);

*___________________________________________________
*
*           PARAMETERS & DATA
*___________________________________________________

PARAMETERS
PMIN        the minimum probability which a scenario n can realize
PMAX        the maximum probability which a scenario n can realize
x(i,n)      the value of the uncertain data for element n of uncertain parameter i
ecdf(i,n)   empirical cumulative density function
jecdf(n)    joint cumulative density function
labels(n)   clustering labels using k-means
MM(i,m)     mth moment calculated from the data
wdash(m)    weights in the constraints
wcdash(i,q)
w(i,m)      --
wc(i,q)
wm(i)       --
temp(m)
CC(i,q)
;


*--------------------------------------
*INPUT FOR THE MODEL ARE:
*x: real values for each uncertain parameter
*ECDF: the value of ECDF
*Labels: the output of the culstering
*--------------------------------------
*$CALL GDXXRW.EXE i=data-bicop.xlsx o=data-bicop.gdx @mapscen.txt
$GDXIN data-bicop.gdx
$LOAD x, ecdf, labels
$GDXIN
;

*---------------------------------------------------------------------------
*   Initialization of the statistical moments - Results from Python
*---------------------------------------------------------------------------

MM('i1','m1')= 15.9376700 ;
MM('i1','m2')= 2.7275301 ;
MM('i1','m3')= -0.8655849*rpower(MM('i1','m2'),1.5) ;
MM('i1','m4')= (0.9794971+3)*rpower(MM('i1','m2'),2) ;
MM('i2','m1')= 11.0479750 ;
MM('i2','m2')= 2.4704010 ;
MM('i2','m3')= 0.7169101*rpower(MM('i2','m2'),1.5) ;
MM('i2','m4')= (0.2062456+3)*rpower(MM('i2','m2'),2) ;

CC('i1','i1') = 2.72753014 ;
CC('i1','i2') = -0.01489189 ;
CC('i2','i1') = -0.01489189 ;
CC('i2','i2') = 2.47040096 ;

*limits for the probabilities
PMIN = 0.01;
PMAX = 0.80;

temp(m) = ord(m) - 1 ;

*Formulation of dynamic set CLN
cln(cl,n) = no ;
labels(n) = labels(n) + 1 ;
loop((cl,n)$[ord(cl) eq labels(n)],
    cln(cl,n) = yes ;
);

DISPLAY x, temp, MM, ecdf, labels, cln;

*___________________________________________________
*
*                   VARIABLES
*___________________________________________________
VARIABLES
yc(cl,n)        is 1 if data point n is selected as the final node cl - otherwise 0
prob(cl,n)      probability of the final node cl with the value of data point n
mplus(i,m)      difference regarding the statistical moments - splitting the variables into two non-negative variables corresponding to the positive and negative values of the original
mminus(i,m)     --
cplus(i,q)
cminus(i,q)
phi(i,n)        deviation of the cumulative probability of final node cl from the ECDF
ed(i)           absolute value of the deviation of final node cl
dx
dm
dc
Z               objective value
;

BINARY VARIABLES yc;
POSITIVE VARIABLES prob, mplus, mminus, cplus, cminus, ed;
phi.lo(i,n) = -1;
phi.up(i,n) = 1;

*WEIGTHS OF THE OBJECTIVE FUNCTION
wdash(m) = 1 ;
wcdash(i,q) = 1 ;
w(i,m) = wdash(m)/(abs(MM(i,m))) ;
wc(i,q) = wcdash(i,q)/(abs(CC(i,q))) ;

w(i,'m1') = 10 ;
w(i,'m2') = 5 ;
w(i,'m3') = 2 ;
w(i,'m4') = 1 ;
wc(i,q) = 3 ;

wm(i) = 0.1;

*A node n that does not belong to cluster cl is not selected for the cluster cl
yc.FX(cl,n)$[not cln(cl,n)] = 0;

Equations eq1   one data point is selected for each cluster - node
          eq2   each data point n is selected for at most 1 prespecified final node
          eq3   sum of probabilities is equal to 0
          eq4   minimum probability
          eq5   maximum probability
          eq6   estimation of the error regarding the mean
          eq7   estimation of the error regarding the higher moments
          eq8   estimation of the deviations from the ECDF for the selected nodes
          eq9   estimation of the maximum deviation from the ECDF for the selected nodes
          eq10  --
          eq11  estimation of the error on correlation matrix
          objm
          objn;


OBJM..          Z =e= SUM( (i,m), w(i,m)*(mplus(i,m)+mminus(i,m)) )  ;

OBJN..          Z =e= SUM((i,m), w(i,m)*(mplus(i,m)+mminus(i,m)) ) + SUM(i, wm(i)*(ed(i)) ) + SUM((i,q), wc(i,q)*(cplus(i,q)+cminus(i,q)) ) ;

eq1(cl)..       SUM((n)$[cln(cl,n)], yc(cl,n)) =e= 1 ;

eq2(n)..        SUM((cl)$[cln(cl,n)], yc(cl,n)) =l= 1 ;

eq3..           SUM((cl,n)$[cln(cl,n)], prob(cl,n)) =e= 1 ;

eq4(cl,n)$[cln(cl,n)]..     yc(cl,n)*Pmin =l= prob(cl,n) ;
eq5(cl,n)$[cln(cl,n)]..     yc(cl,n)*Pmax =g= prob(cl,n) ;

*eq6(i,m)$[ord(m)=1]..   SUM((cl,n)$[cln(cl,n)], prob(cl,n)*x(i,n)) + mplus(i,m) - mminus(i,m) =e= MM(i,m) ;
eq6(i,m)$[ord(m)=1]..   SUM((cl,n)$[cln(cl,n)], prob(cl,n)*x(i,n)) =e= MM(i,m) ;

eq7(i,m)$[ord(m)>1 AND ord(m)<5]..  SUM((n,cl,k)$[ord(k)=1 AND cln(cl,n)], power( [ x(i,n) - MM(i,k) ] , ord(m) )*prob(cl,n) ) + mplus(i,m) - mminus(i,m) =e= MM(i,m) ;

eq8(i,n)..      SUM((cl)$[cln(cl,n)], yc(cl,n)*ecdf(i,n)) - SUM((cl1,nn)$[ [cln(cl1,nn)] AND x(i,nn) le x(i,n)], prob(cl1,nn))  =e=  phi(i,n) ;

eq9(i,n)..      ed(i) =g=  phi(i,n) - ( 1 - SUM((cl)$[cln(cl,n)], yc(cl,n) ) ) ;

eq10(i,n)..     ed(i) =g= -phi(i,n) - ( 1 - SUM((cl)$[cln(cl,n)], yc(cl,n) ) ) ;

eq11(i,q)$[ord(i) lt ord(q)]..     SUM((cl,n)$[cln(cl,n)], ( x(i,n) - MM(i,'m1') )*( x(q,n) - MM(q,'m1') )*prob(cl,n) ) + cplus(i,q) - cminus(i,q) =e= CC(i,q) ;

MODEL MMP /objm,eq2,eq3,eq4,eq5,eq6,eq7/ ;
MODEL DMP /objn,eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11/ ;




PARAMETERS
logp(*,mw,cl)      logger parameters
logx(*,mw,i,cl)    values of scenarios
logm(*,mw,i,k)       values of the moments
loged(*,mw,i)        maximum error regarding the CDF
logobj(*,mw)         objective values
logobjest(*,mw)      objective values
ww(mw)

;

ww('mw1') = 0.1 ;
ww('mw2') = 1 ;
ww('mw3') = 10 ;

OPTION mip = cplex;
OPTION reslim = 1800;
OPTION threads = 8;
OPTION optcr = 0.0;

LOOP((mw)$[ord(mw)=2],

    wm(i) = ww(mw) ;

    solve DMP using mip minimising z;

    logp('DMP',mw,cl) = SUM((n)$[cln(cl,n)], yc.l(cl,n)*prob.l(cl,n) ) ;
    logx('DMP',mw,i,cl) = SUM((n)$[cln(cl,n)], x(i,n)*Yc.l(cl,n) ) ;
    logm('DATA',mw,i,k) = MM(i,k);
    logm('DMP',mw,i,k) = MM(i,k) - mplus.l(i,k) + mminus.l(i,k) ;
    loged('DMP',mw,i) = ed.l(i);
    logobj('DMP',mw) = Z.l ;
    logobjest('DMPC',mw) = DMP.objEst ;

);

option decimals = 5;
DISPLAY  logp, logx, logm, loged, logobj;

$exit

EXECUTE_UNLOAD 'sgresults-final.gdx',
logp, logx, logobj, logm, loged;

EXECUTE 'GDXXRW sgresults-final.gdx par=logx rng=MCDFSL1!B2 rdim=3 cdim=1 par=logp rng=MCDFSL1!B18 rdim=2 cdim=1 par=logobj rng=MCDFSL1!B28 par=logm rng=MCDFSL1!B34' ;
