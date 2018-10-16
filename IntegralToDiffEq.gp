/*======================================================*/ 
/*		1.0  Miscellaneous Utilities 		*/
/*======================================================*/ 

\\ Input: Polynomial function V(z,1/z)
\\ Output: degree of 1/z and z 
{zBiDegree(poly)=my(deg);
deg=[-poldegree(subst(poly,z,1/z),z),poldegree(poly,z)];
if(deg[1]>0,[0,deg[2]],deg)};

\\ Input: Trigonometric Polynomial V(Q,P)
\\ Output: Exponential Polynomial V(1,1/z)
{TrigToExp(TrigPoly)=subst(subst(TrigPoly,
	P,(1/2)*(z+1/z)), 	Q,-I*(1/2)*(z-1/z))};


\\ Input: integer n 
\\ Output: n!! = n*(n-2)*(n-4)...(2 or 1)
{DoubleFactorial(n) = prod(i=0, (n-1)\2, n - 2*i )};

/*======================================================*/ 
/*		2.a  Hyperelliptic Invariants 		*/
/*======================================================*/ 

\\ Input: Polynomial potential function V(q)
\\ Output: Matrix encoding of Hermite reduction 
{HEllGx(poly,deg)=
my(GId,GL,GR);
GId=matrix(2*deg-1,2*deg-1,j,k,
    if(j==k&&k<deg,x,0));
GL=matrix(2*deg-1,2*deg-1,j,k,
    if(0<j-k&&j-k<=deg&&k<deg,
        -2*polcoef(poly,j-k),0));
GR=matrix(2*deg-1,2*deg-1,j,k,
    if(k>=deg&&j-(k-deg)>0 && j-(k-deg)<=deg,
        -(j-(k-deg))*polcoef(poly,j-(k-deg)),0));
GId+GL+GR};

\\ Input: Polynomial potential function V(q)
\\ Output: One-form reduction matrices
{HEllRedData(poly,deg)=
my(ID,InvG,MU,MV,dMV);
ID=matrix(2*deg-1,2*deg-1,j,k,if(j==k,1,0));
InvG=matinverseimage(HEllGx(poly,deg),ID);
MU=matrix(deg-1,deg-1,j,k,InvG[j,k]);
MV=matrix(deg,deg-1,j,k,InvG[j+deg-1,k]);
dMV=matrix(deg-1,deg-1,j,k,j*MV[j+1,k]);
[MU,dMV,MV]};

/*======================================================*/ 
/*		2.b  Trigonometric Invariants 		*/
/*======================================================*/ 

\\ Input: exponential polynomial V(z,1/z)
\\ Output: Matrix encoding of Hermite reduction 
{ExpGx(poly,deg,degTot)=my(GId,GL,GRa,GRb);
GId=matrix(2*degTot-1,2*degTot-1,j,k,
    if(j-k==-deg[1] && k<=degTot,1,0));
GL=matrix(2*degTot-1,2*degTot-1,j,k,
    if(j-k<degTot && j-k >= 0 && k<=degTot,
    polcoeff(-x*poly,j-k+deg[1],z),0));
GRa=matrix(2*degTot-1,2*degTot-1,j,k,
    if(j-(k-degTot)<degTot && j-(k-degTot) >= 0 
    && k<=degTot-deg[1] && degTot< k,
    (j-(k-degTot)+deg[1])*polcoeff(
        -x*poly,(j-(k-degTot)+deg[1]),z),0)); 
GRb=matrix(2*degTot-1,2*degTot-1,j,k,
    if(j-(k-degTot+1)<degTot && j-(k-degTot+1) >= 0 
    && k>degTot-deg[1] ,
    (j-(k-degTot+1)+deg[1])*polcoeff(
        -x*poly,(j-(k-degTot+1)+deg[1]),z),0));    
GId+GL+GRa+GRb};

\\ Input: exponential polynomial P(z,1/z)
\\ Output: One-form reduction matrices 
{ExpRedData(poly,deg,degTot)=my(ID,InvG,MU,MV,dMV);
ID=matrix(2*degTot-1,2*degTot-1,j,k,if(j==k,1,0));
InvG=matinverseimage(ExpGx(poly,deg,degTot),ID);
MU=matrix(degTot,degTot,j,k,InvG[j,k-deg[1] ]);
MV=matrix(degTot-1,degTot,j,k,InvG[j+degTot,k-deg[1]]);
dMV=matrix(degTot,degTot,j,k,if(j<=1-deg[1],
	(j+deg[1]-1)*MV[j,k],(j+deg[1]-1)*MV[j-1,k]));
[MU,dMV,MV]};

/*======================================================*/ 
/*		3.0  Reduction Recursion 		*/
/*======================================================*/ 

\\ Input: One-form and One-form reduction matrices
\\ Output: Hermite-reduced one form 
{HermiteReduce(dtform,cert,dat,n,powm)=if(n==0,
[dtform~,cert~],my(pown);pown=subst(powm,m,n);
HermiteReduce((dat[1]+1/pown*dat[2])*dtform,
cert/*+1/pown/D^pown*dat[3]*dtform*/,dat,n-1,powm))};

\\ Input: One-form reduction matrices
\\ Output: Reduced one-forms.
{PeriodBasis(dat,init,powm)=
my(nOrd,decomposition,HR);
nOrd=1; decomposition=init;
while(matrank(decomposition[1])==nOrd,
HR=HermiteReduce(([1]*init[1])~,([1]*init[2])~,
	dat,nOrd,powm);
nOrd=nOrd+1;decomposition=[
matconcat([decomposition[1];HR[1]]),
matconcat([decomposition[2];HR[2]])];);
concat([nOrd],decomposition)};

/*======================================================*/ 
/*		4.a Hyperelliptic Solve 		*/
/*======================================================*/ 

\\ Input: Polynomial Potential Function V(q)
\\ Output: Annihilator Coefficients
{HEllPicardFuchs(poly)=
my(deg,redDat,init,basis,diagM,AnnCs);
deg=poldegree(poly);
redDat = HEllRedData(poly,deg);
init=[matrix(1,deg-1,j,k,if(k==1,1,0)),
	matrix(1,deg-1,j,k,0)];
basis=PeriodBasis(redDat,init,2*m-1);
diagM=matrix(basis[1],basis[1],j,k,if(j==k,
	((-1/2)^(j-1))*DoubleFactorial(2*j-3),0));
AnnCs=lindep((diagM*basis[2])~)~};

/*======================================================*/ 
/*		4.b Trigonometric Solve 		*/
/*======================================================*/ 

\\ Input: Polynomial(z,1/z) and operating mode
\\ Output: Annihilator Coefficients 
{zPolyIntToDiff(zPoly,mode)=
my(deg,degTot,redDat,init,basis,ltM,AnnCs);
deg=zBiDegree(zPoly);
degTot=deg[2]-deg[1]+1;
redDat=ExpRedData(zPoly,deg,degTot);
init=[matrix(1,degTot,j,k,if(k==1-deg[1],1,0)),
    matrix(1,degTot-1,j,k,0)];
basis=PeriodBasis(redDat,init,
    if(mode==0,m,(2*m-1)/2));
ltM=matrix(basis[1],basis[1],j,k,if(k<=j,
    if(mode==0,(j!)/j/x^(j-1),
        DoubleFactorial(2*j-3)/(2*x)^(j-1)
    )*(-1)^(j+k)*binomial(j-1,j-k),0 ));
AnnCs=lindep((ltM*basis[2])~)~};

\\ Input Heads
SimpleIntToDiff(zPoly)=zPolyIntToDiff(zPoly,0);
{QuarticPicardFuchs(TrigPoly)=
	zPolyIntToDiff(4*TrigToExp(TrigPoly),1)};
