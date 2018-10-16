default(parisize,64000000)

{BinTest(j)=my(Ann,T,vT);
T=sum(n=0,100,binomial(j*n,n)*x^n)+O(x^101);
Ann=SimpleIntToDiff( z^(j-2)*(z+1/z)^(j) );
vT=[T];for(k=1,j-1,vT=concat(vT,[vT[k]']));
[Ann*vT~,Ann]};

\\ Annihilation of C(k*n,n) G.f. for k = 2..10 . 

BinTest(2)
##

BinTest(3)
##

BinTest(4)
##

BinTest(5)
##

BinTest(6)
##

BinTest(7)
##

BinTest(8)
##

BinTest(9)
##

BinTest(10)
##
