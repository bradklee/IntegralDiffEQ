\\ signature 2, genus 1, Edwards' Curve (Cf. OEIS: A002894). 
{T2=sum(j=0,100,binomial(2*j,j)^2*(x/16)^j)+O(x^101);
PFAnn2=QuarticPicardFuchs(P^2*Q^2);
[PFAnn2*[T2,T2',T2'']~,PFAnn2]}
##

\\ signature 3, genus 2, Hyperelliptic (Cf. OEIS: A006480). 
{T3=sum(n=0,100,(3*n)!/(n!)^3*(x/27)^n)+O(x^101);
PFAnn3 = HEllPicardFuchs(-(1/4)*(18*x^2+48*x^4+32*x^6));
[PFAnn3*[T3,T3',T3'']~,-PFAnn3]}
##

\\ signature 4, genus 1, Hyperelliptic (Cf. OEIS: A000897). 
{T4=sum(j=0,100,binomial(2*j,j)*binomial(4*j,2*j)*(x/64)^j)+O(x^101);
PFAnn4a=QuarticPicardFuchs((1/4)*Q^4);
PFAnn4b=HEllPicardFuchs((1/2)*(q^2-(1/4)*q^4));
[PFAnn4a*[T4,T4',T4'']~,-PFAnn4a,PFAnn4a+PFAnn4b]}
##


\\ signature 6, genus 1, Hyperelliptic (Cf. OEIS: A113424). 
{T6=sum(n=0,100,(6*n)!/((3*n)!*(2*n)!*n!)*(1/432*x)^n)+O(x^101);
PFAnn6 = HEllPicardFuchs((1/2)*(x^2-c*x^3));
PFAnn6 = vector(3,j,9/4*(polcoef(PFAnn6[j],2,c)*(4/27)+polcoef(PFAnn6[j],0,c)));
[PFAnn6*[T6,T6',T6'']~,PFAnn6]}
##



