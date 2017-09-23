# Pisot.mpl

# This is a commented version of Doron Zeilberger's Pisot.txt. The comments are
# for my sake, and shouldn't be taken as having any particular insight into the
# mathematics at work.

# I have changed the file extension so that my text editor automatically
# highlights its syntax.

#####################################################################
##Pisot.txt: Save this file as  Pisot.txt                            #
## To use it, stay in the                                            #
##same directory, get into Maple (by typing: maple <Enter> )         #
##and then type:  read  Pisot.txt<Enter>                             #
##Then follow the instructions given there                           #
##                                                                   #
##Written by Doron Zeilberger, Rutgers University ,                  #
#zeilberg at math dot rutgers dot edu                                #
######################################################################
 
#Created: Aug. 2016


with(combinat):
with(numtheory):
with(linalg):

print(`Created:  Aug. 9-20, 2016`):
print(` This is Pisot.txt `):
print(`It accompanyies the article `):
print(` Automated Proof (or disproof) of Linear Recurrences Satisfied by Pisot Sequences`):
print(`by Shalosh B. Ekhad, N. J. A. Sloane and Doron Zeilberger`):
print(`and also available from Zeilberger's website`):
print(``):
print(`Dedicated to Richard Guy (b. 30 Sept. 2016) on his forthcoming 100-th birthday. `):
print(``):
print(`Is is based on an earlier Maple package Cfinite.txt`):
print(``):
print(`Please report bugs to zeilberg at math dot rutgers dot edu`):
print(``):
 print(`The most current version of this  package and paper`):
 print(` are  available from`):
 print(`http://www.math.rutgers.edu/~zeilberg/  .`):



print(`-----------------------`):
 print(`For a list of the Generalized procedures type ezraG();, for help with`):
 print(`a specific procedure, type ezra(procedure_name);   .`):
 print(``):
print(`------------------------------`):

print(`-----------------------`):
 print(`For a list of the MAIN procedures type ezra();, for help with`):
 print(`a specific procedure, type ezra(procedure_name);   .`):
 print(``):
print(`------------------------------`):


print(`-----------------------`):
 print(`For a list of the supporting procedures type ezra1();, for help with`):
 print(`a specific procedure, type ezra(procedure_name);   .`):
 print(``):
print(`------------------------------`):


print(`-----------------------`):
 print(`For a list of the Story procedures from Cfinite type ezraCStory();, for help with`):
 print(`a specific procedure, type ezraC(procedure_name);   .`):
 print(``):
print(`-------------------------------`):

print(`-----------------------`):
 print(`For a list of the procedures coding famous sequences, type`):
 print(`ezraCFamous();`):
 print(` for help with`):
 print(`a specific procedure, type ezraC(procedure_name);   .`):
 print(``):
print(`-------------------------------`):

print(`-----------------------`):
 print(`For a list of the supporting procedures for the procedures from Cfinite.txt, type ezraC1();, for help with`):
 print(`a specific procedure, type ezra(procedure_name);   .`):
 print(``):
print(`-------------------------------`):

print(`-----------------------`):
 print(`For a list of the MAIN procedures from Cfinite type ezraC();, for help with`):
 print(`a specific procedure, type ezra(procedure_name);   .`):
 print(``):
print(`------------------------------`):


ezraG:=proc()

if args=NULL then
 print(` The generalized procedures are: DetSeq, GPr, PisIndG, PSg, PtoRg, ShoAbs `):
else
ezra(args):
fi:

end:

ezra1:=proc()

if args=NULL then
 print(` The supporting procedures are: Boyd, BoydP, Disc2r, Epicycle, FirstDev, GuessPol1, GuessLinear, PisInd1, PtoRvNumbered `):
else
ezra(args):
fi:

end:


ezra:=proc()

if args=NULL then

 print(` The MAIN procedures are: Asy, Binet, Disc2,  DO2,   MamarSyT, MamarSyTtex ,`):
 print(`MamarT, MamarV,Pis, PisInd, PS, PtoR, PtoRb, PtoRsy, PtoRv, Sho, Tikva, UB  `):


elif nargs=1 and args[1]=Asy then
print(`Asy(C): Inputs a C-finite sequence and outputs the pair [A,alpha] such that the n-th term is asymptotic to A*alpha^n`):
print(`Try: `):
print(` Asy([[1,1],[1,1]]); `):
print(` Asy([[10,219,4796,105030],[22,-3,18,-11]]); `):
print(` Asy(DO2([[10,219,4796,105030],[22,-3,18,-11]])); `):

elif nargs=1 and args[1]=Binet then
print(`Binet(C,n): The Binet solution to the C-finite sequence C in terms of n. Try:`):
print(`Binet([[1,1],[1,1]],n);`):


elif nargs=1 and args[1]=Boyd then
print(`Boyd(a0,k,c,r,L,K): a recurrence of order <=L valid up to K terms of the sequence PS(a0,k*a0^2+c,r). Try:`):
print(`Boyd(3,1,1,0,6,1000);`):

elif nargs=1 and args[1]=BoydP then
print(`BoydP(n,a): The Boyd polynomial in a satisfying P(n) is the quotient if P(n-1)^2  by P(n-2) or this quotient`):
print(`minus P(n-2) if this is necessary to make the leading coefficient of the remainder negative. Try:`):
print(`BoydP(5,a);`):

elif nargs=1 and args[1]=DetSeq then
print(`DetSeq(L,r): Inputs a sequence L and outputs a sequence of length nops(L)-2*(r-1) whose`):
print(`i-th entry is det([[L[i], ..., L[i+r-1]], [L[i+1], ..., L[i+r+1]],  ..., [L[i+r-1], ..., L[i+2*(r-1)]].`):
print(`Try: `):
print(`DetSeq([seq(i^4,i=1..10)],3);`):

elif nargs=1 and args[1]=Disc2 then
print(`Disc2(L): inputs a sequence of numbers L and outputs a sequence shorter by 2 whose n-th entry is`):
print(`L[n+1]^2-L[n+2]*L[n]. Try:`):
print(`Disc2(SeqFromRec([[1,1],[1,1]],22),20);`):

elif nargs=1 and args[1]=Disc2r then
print(`Disc2r(L): inputs a sequence of numbers L and outputs a sequence shorter by 2 whose n-th entry is`):
print(`(L[n+1]^2-L[n+2]*L[n])/L[n]. For the sequence L to be a r-Pisot sequence it has to be >=-r and <1-r. Try`):
print(`Disc2r(SeqFromRec([[1,1],[1,1]],22),20);`):

elif nargs=1 and args[1]=DO2 then
print(`DO2(C): inputs a C-finite sequence C, and outputs the C-finite description of its Disc-2-sequence i.e.`):
print(`The linear recurrence equation satisfied by the sequence  L[n+1]^2-L[n+2]*L[n]. Try:`):
print(`DO2([[1,1],[1,1]]);`):


elif nargs=1 and args[1]=Epicycle then
print(`Epicyle(x,y,r,L,K): inputs positive integers x,y, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]`):
print(`and positive integers K1, and K2. First guesses whether there seems to be a linear recurrence equation with constant coefficients`):
print(`(given by its C-finite representation) of order L. If not, it returns FAIL.`):
print(` It is does, it checks whether it goes all the way to K. If it does it returns the recurrence. If it does not `):
print(` it return a more complicated recurrence with an indication how far it goes. Try: `):
print(`Epicycle(5,21,0,5,1000);`):

elif nargs=1 and args[1]=FirstDev then
print(`FirstDev(x,y,r,C): inputs  positive integers x,y, and a parameter r (0<=r<=1) and a C-finite sequence C`):
print(`outputs the first n that violated the recurrence. Try:`):
print(`FirstDev(10,219,1/2,[[10, 219, 4796, 105030], [22, -3, 18, -11]]);`):


elif nargs=1 and args[1]=GuessLinear then
print(`GuessLinear(L,k): Guesses a linear expression a+b*k such that L[i]=a+b*i for i=1..k . `):
print(`L should have at least four terms.  Try:`):
print(`GuessLinear([seq(2*i+5,i=1..10)],k); `):


elif nargs=1 and args[1]=GuessPol1 then
print(`GuessPol1(L,d,n): guesses a polynomial of degree d in n for`):
print(` the list L, such that P(i)=L[i+1] for i=1..nops(L)`):
print(`For example, try: `):
print(`GuessPol1([seq(i,i=1..10)],1,n);`):

elif nargs=1 and args[1]=GPr then
print(`GPr(x,n): The template for the generalized Pisot recurrence of order n, expressing x[0] in terms of x[-1], ..., x[-n]`):
print(` such that the nby n determinant det([[x[-2*(n-1)], ..., x[-(n-1)]],[x[-2*(n-1)+1], ...,x[-n+2]], ..., [x[-(n-1)], ..., x[0]]. `):
print(` Try: `):
print(` GPr(x,3); `):

elif nargs=1 and args[1]=MamarSyT then
print(`MamarSyT(X,k,r,Ord) : inputs a positive integer X, a symbol k, a positive integer r striclty between 0 and 1`):
print(`(typically 1/2), a positive integer Ord, outputs all the recurrences for the infinite families`):
print(`PS(x,y0+k*x^2,r) in terms of the symbol k, valid for k>=1 (and if y0>x, probably also for k=0)  for all y0`):
print(`from 1 to x^2-1 that are NOT multiples of x, for all 2<=x<=X. Try:`):
print(` MamarSyT(5,k,1/2,6);`):

elif nargs=1 and args[1]=MamarSyTtex then
print(`MamarSyTtex(X,k,r,Ord)  : like MamarSyT(X,k,r,Ord) (q.v.) by with TeX output. Try:`):
print(` MamarSyTtex(5,k,1/2,6);`):

elif nargs=1 and args[1]=MamarT then
print(`MamarT(K1,K2,r,Ord,HAMON): Inputs  positive integers K1<K2, a rational number r, and a positive integer Ord. Outputs`):
print(`all the pairs  2<=a<K1, a+2<=b<=K2 such that b/a is not an integer and such that the Pisot sequence with parameter r`):
print(`and initial conditions a,b, has a provable recurrence of order <=Ord. It also returns the`):
print(`successes and failues.`):
print(`It returns a list of length 4, consisting of`):
print(`the set of pairs that produce provable non-trivial recurrences`):
print(` those giving linear sequences. Those that start-out satisfying a recurrence`):
print(`of order <=Ord, but ultimately fail, and those that fail right away. Try:`):
print(`MamarT(20,20,1/2,6,60);`):



elif nargs=1 and args[1]=MamarV then
print(`MamarV(K1,K2,r,Ord): Verbose version of MamarT(K1,K2,r,Ord) (q.v.) complete with proofs. `):
print(`all the pairs  2<=a<K1, a+2<=b<=K2 such that b/a is not an integer and such that the Pisot sequence with parameter r`):
print(`and initial conditions a,b, has a provable recurrence of order <=Ord. It also returns the`):
print(`successes and failues.`):
print(`It returns a list of length 4, consisting of`):
print(`the set of pairs that produce provable non-trivial recurrences`):
print(` those giving linear sequences. Those that start-out satisfying a recurrence`):
print(`of order <=Ord, but ultimately fail, and those that fail right away. Try:`):
print(`MamarV(10,10,1/2,6);`):


elif nargs=1 and args[1]=Pis then
print(`Pis(C): Inputs a C-finite sequence and outputs the absolute value of the second-largest root`):
print(`It is a Pisot number if it is less than 1.`):
print(`Pis([[1,1],[1,1]]);`):
print(`Pis([[10,219,4796,105030],[22,-3,18,-11]]);`):
print(`Pis([[5,21,88,368,1538],[4,1,-1,0,-1]]); `):

elif nargs=1 and args[1]=PisInd then
print(`PisInd(C,r,K1,K2): inputs a C-finite description C, describing a sequence, let's call it L, and a positive integer K, outputs`):
print(`the min and max of the sequence, from n=K1 to K2, of the quantity (L[n+1]^2-L[n+2]*L[n])/L[n]+r `):
print(`They should all be >=0 and <1 respectively. Try:`):
print(`PisInd([[2,3],[1,1]],1/2,1,100); `):
print(`PisInd([[4,7,12],[2,-1,1]],1/2,50,100); `):
print(`PisInd([[5,21,88,368,1538],[4,1,-1,0,-1]],0,1,28); `):


elif nargs=1 and args[1]=PisInd1 then
print(`PisInd1(C,r,K): inputs a C-finite description C, describing a sequence, let's call it L, and a positive integer K, outputs`):
print(`the sequence, from n=1 to K, of the quantity (L[n+1]^2-L[n+2]*L[n])/L[n]+r for n from 1 to K. `):
print(`They should all be >=0 and <1 respectively. Try:`):
print(`PisInd1([[2,3],[1,1]],1/2,100); `):
print(`PisInd1([[4,7,12],[2,-1,1]],1/2,100); `):
print(`PisInd1([[5,21,88,368,1538],[4,1,-1,0,-1]],0,28); `):

elif nargs=1 and args[1]=PisIndG then
print(`PisIndG(C,k,r,K1,K2): inputs a C-finite description C, describing a sequence, let's call it L, and  positive integers K1, K2, outputs`):
print(`the min and max from n=K1 to K2 of the quantity DetSeq(L,k)/DetSeq(L,k-1).`):
print(`They should be >=0 and <1 respectively. Try:`):
print(`PisIndG([[2,3],[1,1]],2,1/2,100,200);`):
print(`PisIndG([[4,7,12],[2,-1,1]],2,1/2,100,200);`):


elif nargs=1 and args[1]=PSg then
print(`PSg(INI,r,K): the first K terms of the Generalized Pisot sequence of order n, defined by`):
print(`a[i]=INI[i], where INI is a list of even length 2*(n-1)`):
print(`PSg([2,5,9,10],1/2,20);`):

elif nargs=1 and args[1]=PtoR then
print(`PtoR(x,y,r,L,K): inputs positive integers x,y, and a rational number r STRICTLY between 0 and 1  [typically 1/2]`):
print(`and a positive integer K. First guesses whether there seems to be a linear recurrence equation with constant coefficients`):
print(` (given by its C-finite representation) of order L. If not, it returns FAIL. `):
print(` If it persists to K terms (make K large), then it returns `):
print(`the C-finite description, followed by its Pisot-index-pair for these K terms, and the ones for the second half.`):
print(`(so it returns a list of length 4 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex], PisotIndicator ]`):
print(`If the initial conjecture fails, then it returns  [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:`):
print(`PtoR(2,3,1/2,10,1000);`):
print(`PtoR(5,17,1/2,5,1000);`):
print(`PtoR(4,7,1/2,10,5000);`):
print(`PtoR(10,219,1/2,10,2000);`):


elif nargs=1 and args[1]=PtoRb then
print(`PtoRb(x,y,r,L,K): inputs positive integers x,y, and a rational number r  between 0 and 1 (inclusive)  [typically 0 or 1]`):
print(`(for 0<r<1 please use PtoR(x,y,r,L,K))`):
print(`and positive integers K. First guesses whether there seems to be a linear recurrence equation with constant coefficients`):
print(` (given by its C-finite representation) of order L. If not, it returns FAIL. `):
print(` If it persists to K terms (make K large), then it returns `):
print(`the C-finite description, followed by its Pisot-index-pair for these K terms, and the ones for K from 2K .`):
print(`(so it returns a list of length 4 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex], PisotIndicator ]`):
print(`If the initial conjecture fails, then it returns  [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:`):
print(`PtoRb(5,21,0,5,1000);`):
print(`PtoRb(8,55,1,5,1000);`):

elif nargs=1 and args[1]=PtoRg then
print(`PtoRg(INI,r,L,K): inputs a list of positive integers INI, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]`):
print(`and positive integer K. First guesses whether there seems to be a linear recurrence equation with constant coefficients`):
print(`(given by its C-finite representation) of order L for PSg(INI,r). If not, it returns FAIL.`):
print(`If it persists to K terms (make K large), then it returns`):
print(`the C-finite description, followed by its generalized`):
print(`Pisot-index-pair for these K terms, and the ones for the second half.`):
print(`(so it returns a list of length 3 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex, PisotIndicator] ]`):
print(`If the initial conjecture fails, then it returns it returns [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:`):
print(`PtoRg([2,3],1/2,10,1000);`):
print(`PtoRg([5,17],1/2,10,1000);`):
print(`PtoRg([1,4,,8,23],1/2,10,1000);`):


elif nargs=1 and args[1]=PtoRsy then
print(`PtoRsy(x,y0,r,k,ORD,Kama): inputs a positive integer x, another one, y0, between 1 and x^2-1 (not divisible by x),`):
print(`a rational number r, between 0 and 1, a positive integer, ORD, at least 2, and a positive integer Kama, at least 4`):
print(`tries (by doing Kama cases), to conjecture a  recurrence of order<=ORDER whose coefficients are linear in k`):
print(`that  for the Pisot Sequence E_r(x,x^2*k+y0,r)`):
print(` valid for all k>=0, just giving "infinitely" many facts. It returns FAIL. Try: `):
print(`PtoRsy(2,1,1/2,k,2,10);`):

elif nargs=1 and args[1]=PtoRv then
print(`PtoRv(x,y,r,L,N,K): Verbose version of PtoR(x,y,r,L,K) (q.v.).`):
print(`Outputs a theorem with a sketch of a proof, or a satement that a recurrence of order<=L does not exist.`):
print(`and also printing out the first N terms of the Pisot sequence PS(x,y,r)`):
print(` Try: `):
print(`PtoRv(2,3,1/2,10,30,1000);`):
print(`PtoRv(5,17,1/2,10,30,1000);`):
print(`PtoRv(5,21,0,5,30,1000);`):
print(`PtoRv(4,7,1/2,10,30,5000);`):
print(`PtoRv(10,219,1/2,10,30,2000);`):
print(`PtoRv(8,55,1,10,30,12000);`):

elif nargs=1 and args[1]=PtoRvNumbered then
print(`PtoRvNumbered(x,y,r,L,N,K,MISPAR): Like PtoRv(x,y,r,L,N,K) (q.v.) but with the theorem or fact numbered by MISPAR`):
print(`called-on by procedure MamarV (q.v.)`):
print(` Try: `):
print(`PtoRvNumbered(2,3,1/2,10,30,1000,41);`):
print(`PtoRvNumbered(5,17,1/2,10,30,1000,41);`):
print(`PtoRvNumbered(5,21,0,5,30,1000,41);`):
print(`PtoRvNumbered(4,7,1/2,10,30,5000,41);`):
print(`PtoRvNumbered(10,219,1/2,10,30,2000,41);`):
print(`PtoRvNumbered(8,55,1,10,30,12000,41);`):

elif nargs=1 and args[1]=PS then
print(`PS(x,y,r,K): the first K terms of the Pisot-like sequence defined by`):
print(`a[1]=x, a[2]=y, a[n]=[ a(n-1)^2/a(n-2) +r] (where [t] means integer-part). For example, for the first 20 terms op the`):
print(`Fibonacci sequence`):
print(`(without 0,1,1) type:`):
print(`PS(2,3,1/2,20);`):

elif nargs=1 and args[1]=Sho then
print(`Sho(C): all the roots of the characteristic equation of the recurrence F. Try:`):
print(`Sho([[1,1],[1,1]]);`):
print(`Sho([[10,219,4796,105030],[22,-3,18,-11]]);`):
print(`Sho([[5,21,4796,105030],[22,-3,18,-11]]);`):
print(`Sho([[5,21,88,368,1538],[4,1,-1,0,-1]]); `):

elif nargs=1 and args[1]=ShoAbs then
print(`ShoAbs(C): the sorted (from largest to smallest of the absolute values of the characteristic equation of the recurrence C`):
print(`Try:`):
print(`ShoAbs([[1,1],[1,1]]);`):

elif nargs=1 and args[1]=Tikva then
print(`Tikva(C,N):  the first N terms of the sequence (L[n+1]^2-L[n+2]*L[n])/L[n] where L[n] is the C-finite sequence given by`):
print(`C. For it to be an r-Pisot sequence EVERY term must be >=-r and <1-r. Try:`):
print(`Tikva([[5,21,88,368,1538],[4,1,-1,0,-1]],27);`):

elif nargs=1 and args[1]=UB then
print(`UB(C): Inputs a C-finite sequence and outputs the largest absolute value of the roots of the characteristic equation`):
print(` Try: `):
print(` UB([[1,1],[1,1]]); `):


else

print(`There is no ezra for`, args):
fi:
end:

###start from Cfinite.txt
ezraCFamous:=proc():
print(`The procedures coding famous sequences are`):
print(`fn,Fn, Ln, Pn, Ux,Trn, Tx`):
end:

ezraCStory:=proc()

if args=NULL then
 print(` The Story procedures are: BTv,FactorizeI1v, FindCHv, GuessNLRv, `):
 print(` KefelV, mSectV, SeferGW `):
else
ezra(args):
fi:

end:

ezraC1:=proc()

if args=NULL then
 print(` The supporting procedures are: CartProd, CartProd2 `):
 print(` CHseq, CtoR, Curt1,`):
 print(` etop, GenPol,  FindCH1, FindCHg, GuessNLR1, GuessPR1`):
 print(` GuessRec1, GuessRec1mini,IISF, ISF, Krovim, ptoe `):
 print(`  ProdIndicator, ProdIndicatorG, ProdProfile, `):
 print(` ProdProfileE, RtoC, TDB`):

else
ezraC(args):
fi:

end:


ezraC:=proc()

if args=NULL then
 print(`The main procedures are: BT,  CH, Factorize, FactorizeI1`):
 print(` FindCH,GuessNLR, GuessPR, GuessRec, IsProd, IsProdG `):
 print(` Kefel, KefelR, KefelS, KefelSm, Khibur, KhiburS`):
 print(` mSect, SeqFromRec `):
 print(` `):



elif nops([args])=1 and op(1,[args])=BT then
print(`BT(C): The Binomial Transform of the C-finite sequence C,`):
print(`in other words, the sequence add(binomial(n,i)*C[i],i=0..n):`):
print(`For example, try:`):
print(`BT([[0,1],[1,1]]);`):

elif nops([args])=1 and op(1,[args])=BTv then
print(`BTv(C): Verbose form of BT(C)`):
print(`For example, try:`):
print(`BTv([[0,1],[1,1]]);`):

elif nops([args])=1 and op(1,[args])=CartProd then
print(`CartProd(L): The Cartesian product of the lists`):
print(`in the list of lists L`):
print(`for example, try:`):
print(`CartProd([[2,5],[4,9],[-2,-7]]);`):

elif nops([args])=1 and op(1,[args])=CartProd2 then
print(`CartProd2(L1,L2): The Cartesian product of lists L1 and L2,`):
print(`for example, try:`):
print(`CartProd2([2,5],[4,9]);`):

elif nops([args])=1 and op(1,[args])=CH then
print(`CH(L,n,j): inputs  a list of pairs [C,e], where C is a C-finite`):
print(`sequence, and e is an affine-linear expression in the symbols`):
print(`n and j, of the form a*n+b*j+c, where a and a+b are non-negative`):
print(`integers. Outputs the C-finite description of`):
print(`Sum(L[1][1](L[1][2])*L[2][1](L[2][2]),j=0..n-1)`):
print(`Implements the sequences discussed by Curtis Greene and`):
print(`Herb Wilf ("Closed form summation of C-finite sequences",`):
print(`Trans. AMS.)`):
print(`For example, try:`):
print(`CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j);`):

elif nops([args])=1 and op(1,[args])=CHseq then
print(`CHseq(L,n,j,K): inputs  a list of pairs [C,e], where C is a C-finite`):
print(`sequence, and e is an affine-linear expression in the symbols`):
print(`n and j, of the form a*n+b*j+c, where a and a+b are non-negative`):
print(`integers. Outputs the first K terms of`):
print(`Sum(L[1][1](L[1][2])*L[2][1](L[2][2]),j=0..n-1)`):
print(`using up to K values. For example, try:`):
print(`CHseq([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j,20);`):

elif nops([args])=1 and op(1,[args])=CtoR then
print(`CtoR(S,t), the rational function, in t, whose coefficients`):
print(`are the members of the C-finite sequence S. For example, try:`):
print(`CtoR([[1,1],[1,1]],t);`):

elif nops([args])=1 and op(1,[args])=Curt1 then
print(`Curt1(j,n,M): outputs `):
print(`all {a1n+b1j} with 0<=a1<=M and 0<=a1+b1<=M`):
print(`Type: Curt1(j,n,4);`):


elif nops([args])=1 and op(1,[args])=etop then
print(`etop(e): If e is a list of numbers (or expressions)`):
print(`representing the elementary symmetric functions`):
print(`outputs the power-symmetric functions`):
print(`For example, try`):
print(`etop([1,3,6]);`):

elif nops([args])=1 and op(1,[args])=Factorize then
print(`Factorize(C,L): Given a C-finite recurrence `):
print(`(without the initial conditions)`):
print(`tries to find a factorization into C-finite sequences of order`):
print(`L[i],i=1..nops(L).`):
print(`It returns all the factoriziations, if found, followed by `):
print(`the free parameters and FAIL otherwise`):
print(` For example, try:`):
print(`Factorize([-3,106,-72,-576],[2,2]);`):
print(`Factorize([2,7,2,-1],[2,2]);`):

elif nops([args])=1 and op(1,[args])=FactorizeI1 then
print(`FactorizeI1(C,k): Given a C finite integer sequence,C`):
print(`tries to factorize it into a C finite sequence of`):
print(`order k and order nops(C[1])-k. For example, try:`):
print(`FactorizeI1([[2,10,36,145],[2,7,2,-1]],2);`):
print(`FactorizeI1([[1,8,95,1183],[15,-32,15,-1]],2);`):
print(`FactorizeI1(RtoC(TDB(t)[3],t),2);`):

elif nops([args])=1 and op(1,[args])=FactorizeI1v then
print(`FactorizeI1v(C,k): verbose version of FactorizeI1(C,k)`):
print(` For example, try:`):
print(`FactorizeI1v([[2,10,36,145],[2,7,2,-1]],2);`):
print(`FactorizeI1v([[1,8,95,1183],[15,-32,15,-1]],2);`):
print(`FactorizeI1v(RtoC(TDB(t)[3],t),2);`):

elif nops([args])=1 and op(1,[args])=Fn then
print(`Fn(): The C-finite description of the (usual) Fibonacci sequence`):

elif nops([args])=1 and op(1,[args])=fn then
print(`fn(): The C-finite description of the fibonacci sequence `):
print(`fibonacci(n+1)`):

elif nops([args])=1 and op(1,[args])=GenPol then
print(`GenPol(Resh,a,deg): Given a list of variables Resh, a symbol a,`):
print(`and a non-negative integer deg, outputs a generic polynomial`):
print(`of total degree deg,in the variables Resh, indexed by the letter a of`):
print(`followed by the set of coeffs. For example, try:`):
print(`GenPol([x,y],a,1);`):


elif nops([args])=1 and op(1,[args])=FindCH then
print(`FindCH(C1,r,d, x,L,n,j): Inputs a C-finite sequence C1 and pos. integers`):
print(`r and d, and a letter x, and a list L of affine-linear expressions`):
print(`in n,j, and`):
print(`and outputs a polynomial P(x[1],..., x[r]) of degree<= d such that `):
print(`P(C1[n],C1[n-1],C1[n-r])=C2[n]`):
print(`where C2(n) is the sum(C1[L[1]]*C1[L[2]*...C1[L[nops(L)],j=0..n-1):`):
print(`If none found it returns FAIL. For example to reproduce`):
print(`Eq. (12) of the Greene-Wilf article.`):
print(`If there is none of degree <=d then it returns FAIL`):
print(` type:`):
print(`FindCH(Fn(),1,3,x,[j,j,2*n-j],n,j);`):

elif nops([args])=1 and op(1,[args])=FindCHv then
print(`FindCHv(C1,r,d, x,L,n,j):  verbose version of`):
print(`FindCH(C1,r,d, x,L,n,j)`):
print(` type:`):
print(`FindCHv(Fn(),1,3,[j,j,2*n-j],n,j);`):

elif nops([args])=1 and op(1,[args])=FindCHg then
print(`FindCHg(C1,r,d,x,C2): Inputs a C-finite sequence C1 and`):
print(` a pos. integer`):
print(`d and outputs a polynomial P(x[1],..., x[r]) of minimal`):
print(` degree  such that `):
print(`P(C1[n],C1[n-1],C1[n-r])=C2[n]`):
print(`for all n. If none found it returns FAIL. For example to reproduce`):
print(`Eq. (12) of the Greene-Wilf article,  type:`):
print(`FindCHg([[0,1],[1,1]],1,3,x,CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j));`):

elif nops([args])=1 and op(1,[args])=GuessCH1 then
print(`GuessCH1(C1,r,d,x,C2): Inputs a C-finite sequence C1 and`):
print(` a pos. integer`):
print(`d and outputs a polynomial P(x[1],..., x[r]) of degree d such that `):
print(`P(C1[n],C1[n-1],C1[n-r])=C2[n]`):
print(`for all n. If none found it returns FAIL. For example to reproduce`):
print(`Eq. (12) of the Greene-Wilf article,  type:`):
print(`GuessCH1([[0,1],[1,1]],1,3,x,CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j));`):

elif nops([args])=1 and op(1,[args])=GuessNLRv then
print(`GuessNLRv(C1,r,d): Verbose form of GuessNLR(C1,r,d,x)`):
print(`Inputs a C-finite sequence C and `):
print(`a pos. integer d and outputs a polynomial relation`):
print(` of todal degree<= d such that P(C1[n],C1[n-1],C1[n-r])=0`):
print(`for all n. If none found it returns FAIL. For example to solve`):
print(`the Hugh Montgomery proposed (but not selected) Putnam problem,`):
print(`(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,`):
print(`where F(n) is the Fibonacci sequence type`):
print(`GuessNLRv([[0,1],[1,1]],1,10);`):

elif nops([args])=1 and op(1,[args])=GuessNLR then
print(`GuessNLR(C1,r,d,x): Inputs a C-finite sequence C and `):
print(`a pos. integer d and outputs a polynomial P(x[0],..., x[r])`):
print(` of todal degree<= d such that P(C1[n],C1[n-1],C1[n-r])=0`):
print(`for all n. If none found it returns FAIL. For example to solve`):
print(`the Hugh Montgomery proposed (but not selected) Putnam problem,`):
print(`(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,`):
print(`where F(n) is the Fibonacci sequence type`):
print(`GuessNLR([[0,1],[1,1]],1,10,x);`):

elif nops([args])=1 and op(1,[args])=GuessNLR1 then
print(`GuessNLR1(C1,r,d,x): Inputs a C-finite sequence C and `):
print(`a pos. integer d and outputs a polynomial P(x[0],..., x[r])`):
print(` of degree d such that P(C1[n],C1[n-1],C1[n-r])=0`):
print(`for all n. If none found it returns FAIL. For example to solve`):
print(`the Hugh Montgomery proposed (but not selected) Putnam problem,`):
print(`(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,`):
print(`where F(n) is the Fibonacci sequence type`):
print(`GuessNLR1([[0,1],[1,1]],1,4,x);`):

elif nops([args])=1 and op(1,[args])=GuessPR then
print(`GuessPR(C1,C2,d,x,y): `):
print(`Inputs two C-finite sequences and a pos. integer`):
print(`d and outputs a polynomial P(x,y), of lowest possible degree`):
print(` <=d such that P(C1[n],C2[n])=0`):
print(`for all n. If none found it returns FAIL. For example to solve`):
print(`the Hugh Montgomery proposed (but not selected) Putnam problem,`):
print(`(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,`):
print(`where F(n) is the Fibonacci sequence type`):
print(`GuessPR([[0,1],[1,1]],[[1,1],[1,1]],10,x,y);`):

elif nops([args])=1 and op(1,[args])=GuessPR1 then
print(`GuessPR1(C1,C2,d,x,y): `):
print(`Inputs two C-finite sequences and a pos. integer`):
print(`d and outputs a polynomial P(x,y)`):
print(` of degree d such that P(C1[n],C2[n])=0`):
print(`for all n. If none found it returns FAIL. For example to solve`):
print(`the Hugh Montgomery proposed (but not selected) Putnam problem,`):
print(`(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,`):
print(`where F(n) is the Fibonacci sequence type`):
print(`GuessPR1([[0,1],[1,1]],[[1,1],[1,1]],4,x,y);`):

elif nops([args])=1 and op(1,[args])=GuessRec then
print(`GuessRec(L): inputs a sequence L and tries to guess`):
print(`a recurrence operator with constant cofficients`):
print(`satisfying it. It returns the initial  values and the operator`):
print(` as a list. For example try:`):
print(`GuessRec([1,1,1,1,1,1]);`):

elif nops([args])=1 and op(1,[args])=GuessRec1 then
print(`GuessRec1(L,d): inputs a sequence L and tries to guess`):
print(`a recurrence operator with constant cofficients of order d`):
print(`satisfying it. It returns the initial d values and the operator`):
print(` as a list. For example try:`):
print(`GuessRec1([1,1,1,1,1,1],1);`):

elif nops([args])=1 and op(1,[args])=GuessRec1mini then
print(`GuessRec1mini(L,d): inputs a sequence L and tries to guess`):
print(`a recurrence operator with constant cofficients of order d`):
print(`satisfying it. It returns the initial d values and the operator`):
print(`It is like GuessRec1(L,d) but allows shorter lists`):
print(` as a list. For example try:`):
print(`GuessRec1mini([1,1,1,1,1,1],1);`):


elif nops([args])=1 and op(1,[args])=IISF then
print(`IISF(L,k): given a sequence of integers finds all its`):
print(`Increasing integer`):
print(`factors of length k. For example try:`):
print(`IISF([1,4,9,16,32],3);`):

elif nops([args])=1 and op(1,[args])=ISF then
print(`ISF(L,k): given a sequence of integers finds all its`):
print(`factors of length k. For example try:`):
print(`ISF([1,4,9,16,32],3);`):

elif nops([args])=1 and op(1,[args])=IsProd then
print(`IsProd(f,t,m,n,eps): Given a rational function f of t`):
print(`with NUMERICAL coefficients`):
print(`decides, empirically,`):
print(` whether the C-finite sequence of its coeeficients`):
print(`is a product of a C-finite sequences of order m and n.`):
print(`For example, try:`):
print(`IsProd(1/((1-t)*(1-2*t)*(1-3*t)*(1-6*t)),t,2,2,1/1000);`):


elif nops([args])=1 and op(1,[args])=IsProdG then
print(`IsProdG(f,t,L,eps): Given a rational function f of t`):
print(`with NUMERICAL coefficients`):
print(`and a list of positive integers, L`):
print(`decides, empirically,`):
print(`decides whether the C-finite sequence of its coeeficients`):
print(`is a product of a C-finite sequences of orders given by the`):
print(`members of L`):
print(`For example, try:`):
print(`IsProdG(1/((1-t)*(1-2*t)*(1-3*t)*(1-6*t)),t,[2,2],1/1000);`):
print(`IsProdG(CtoR([[2,10,36,145],[2,7,2,-1]],t),t,[2,2],1/100000);`):
print(`IsProdG(TDB(t)[1],t,[2,2],1/100000);`):
print(`IsProdG(TDB(t)[2],t,[2,2,2],1/100000);`):
print(`IsProdG(TDB(t)[3],t,[2,2,2,2],1/100000);`):
print(`IsProdG(TDB(t)[4],t,[2,2,2,2,2],1/100000);`):

elif nops([args])=1 and op(1,[args])=Kefel then
print(`Kefel(S1,S2): inputs two C-finite sequences S1,S2, and  outputs`):
print(`their product, For example, try:`):
print(`Kefel([[1,2],[3,4]],[[2,5],[-1,6]]);`):

elif nops([args])=1 and op(1,[args])=KefelV then
print(`KefelV(C1,C2): Verbose version of Kefel(S1,S2);`):
print(`For example try:`):
print(`KefelV([[1,2],[3,4]],[[2,5],[-1,6]]);`):

elif nops([args])=1 and op(1,[args])=KefelR then
print(`KefelR(R1,R2,t): inputs rational functions of t,`):
print(`R1 and R2,`):
print(`and outputs the rational function `):
print(`whose coeffs. is their Hadamard product, for example try:`):
print(`KefelR(1/(1-t-t^2),1/(1-t-t^3),t);`):


elif nops([args])=1 and op(1,[args])=KefelS then
print(`KefelS(a,b,d1,d2): The recurrence operator`):
print(`satisfied by the product of any solution of`):
print(`the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]`):
print(`and any solution of the recurrence `):
print(` g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1] `):
print(`For example, try:`):
print(`KefelS(a,b,2,2);`):

elif nops([args])=1 and op(1,[args])=KefelSm then
print(`KefelSm(a,L): The recurrence operator`):
print(`satisfied by the product of nops(L) recurrences`):
print(`the recurrence f[i][n]=a[i,1]*f[n-1]+...+a[i,L[i]]*f[n-L[i]]`):
print(`It also returns the set of generic coefficients.`):
print(`For example, try:`):
print(`KefelSm(a,[2,2,2]);`):

elif nops([args])=1 and op(1,[args])=Khibur then
print(`Khibur(S1,S2): inputs two C-finite sequences S1,S2, and  outputs`):
print(`their sum, For example, try:`):
print(`Khibur([[1,2],[3,4]],[[2,5],[-1,6]]);`):


elif nops([args])=1 and op(1,[args])=KhiburS then
print(`KhiburS(a,b,d1,d2): The recurrence operator`):
print(`satisfied by the sum of any solution of`):
print(`the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]`):
print(`and any solution of the recurrence `):
print(` g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1]`):
print(`For example, try:`):
print(`KhiburS(a,b,2,2);`):

elif nops([args])=1 and op(1,[args])=Krovim then
print(`Krovim(L,i,eps): Given a list L and an index between`):
print(`1 and nops(L), finds the set of indices (including i)`):
print(`such that abs(L[j]-L[i])<eps. For example, try:`):
print(`Krovim([1.0001,1.01,1.00011],1,1/1000);`):

elif nops([args])=1 and op(1,[args])=Ln then
print(`Ln(): The C-finite description of the Lucas number sequence`):

elif nops([args])=1 and op(1,[args])=mSect then
print(`mSect(C1,m,i): Given a C-finite sequence C1(n), finds`):
print(`the C-finite description of C(m*n+i) where m is a pos. integer`):
print(`and 0<=i<m. For example, try:`):
print(`mSect([[0,1],[1,1]],2,0);`):

elif nops([args])=1 and op(1,[args])=mSectV then
print(`mSectV(C1,m,i): Verbose form of mSect(C1,m,i)`):
print(`Given a C-finite sequence C1(n), finds`):
print(`the C-finite description of C(m*n+i) where m is a pos. integer`):
print(`and 0<=i<m. For example, try:`):
print(`mSectV([[0,1],[1,1]],2,0);`):

elif nops([args])=1 and op(1,[args])=Pn then
print(`Pn(): The C-finite description of the Pell Numbers `):

elif nops([args])=1 and op(1,[args])=ptoe then
print(`ptoe(p): If p is a list of numbers (or expressions)`):
print(`representing the power-sum symmetric functions`):
print(`outputs the elementary symmetic functions`):
print(`For example, try`):
print(`ptoe([1,3,6]);`):

elif nops([args])=1 and op(1,[args])=ProdIndicator then
print(`ProdIndicator(m,n): The profile of multiplicity of`):
print(`the (m*n)^2 ratios of a list of m*n numbers to`):
print(`be the Cartesian product of a set of m numbers with`):
print(`a set of n numbers. For example, try:`):
print(`ProdIndicator(2,2);`):

elif nops([args])=1 and op(1,[args])=ProdIndicatorG then
print(`ProdIndicatorG(L): The profile of multiplicity of`):
print(`the convert(L,`*`)^2 ratios of a list of `):
print(`be the Cartesian product of lists of size given by L`):
print(`For example, try:`):
print(`ProdIndicatorG([2,2]);`):

elif nops([args])=1 and op(1,[args])=ProdProfile then
print(`ProdProfile(L,eps): Given a list of complex numbers finds its`):
print(`(approximate) profile to within eps. For example, try:`):
print(`ProdProfile([1,2,3,6],1/1000);`):

elif nops([args])=1 and op(1,[args])=ProdProfileE then
print(`ProdProfileE(L): Given a list of numbers (or expressions) finds its`):
print(`(exact) profile to within eps. For example, try:`):
print(`ProdProfileE([1,2,3,6]);`):



elif nops([args])=1 and op(1,[args])=RtoC then
print(`RtoC(R,t): Given a rational function R(t) finds the`):
print(`C-finite description of its sequence of coefficients`):
print(`For example, try:`):
print(`RtoC(1/(1-t-t^2),t);`):


elif nops([args])=1 and op(1,[args])=SeferGW then
print(`SeferGW(m,d,M): outputs a book of Fibonacci identities with`):
print(`right hand sides of degree at most d in F(n),F(n+1)`):
print(`for Greene-Wilf type sums where the summand is`):
print(`of the form F(j)*F(a1n+b1j)*... with m terms, for example, try:`):
print(`SeferGW(2,5,2);`):


elif nops([args])=1 and op(1,[args])=SeqFromRec then
print(`SeqFromRec(S,N0): Given a C-finite sequence S, written as`):
print(`S=[INI,ope], where INI is the list of initial conditions`):
print(`ope is the  recurrence operator (linear with constant coefficients)`):
print(`codes a list of the same size as INI`):
print(`(the recurrence is: f(n)=ope[1]*f(n-1)+ope[2]*f(n-2)+ ...) , `):
print(` finds the first N0 terms of the sequence satisfying`):
print(`the recurrence f(n)=ope[1]*f(n-1)+...ope[L]*f(n-L).`):
print(`For example, for the first 20 Fibonacci numbers,`):
print(`try: SeqFromRec([[1,1],[1,1]],20);`):

elif nops([args])=1 and op(1,[args])=TDB then
print(`TDB(t): the list of rational functions in t`):
print(`of length 4, consisting of the generating functions`):
print(`for the number of perfect matchings in a `):
print(`4 by n, 6 by n, 8 by n and 10 by n`):
print(`Do TDB(t);`):

elif nops([args])=1 and op(1,[args])=Trn then
print(`Trn(): The C-finite description of the Triboanacci numbers Numbers`):

elif nops([args])=1 and op(1,[args])=Tx then
print(`Tx(x): The C-finite description of the Chebychev polynomials of`):
print(`the first kind. Try: Tx(x);`):

elif nops([args])=1 and op(1,[args])=Ux then
print(`Ux(x): The C-finite description of the Chebychev polynomials of`):
print(`the second kind. Try: Ux(x);`):

else
print(`There is no ezra for`,args):
fi:
 
end:


TDB:=proc(t):
[-(1+t)/(-5*t^2-t+1-t^3+t^4)*(t-1), -(t^6-8*t^4+2*t^3+8*t^2-1)/(1+t)/(t-1)/(t^6+t^5-19*t^4+11*t^3+19*t^2+t-1), -(-1-1033*t^8+110*t^9+26*t^3+43*t^2-360*t^4-26*t^11+t
^14+1033*t^6-43*t^12-110*t^5+360*t^10)/(t^16-t^15-76*t^14-69*t^13+921*t^12+584*t^11-4019*t^10-829*t^9+7012*t^8-829*t^7-4019*t^6+584*t^5+921*t^4-69*t^3-76*t^2-t+1),
-(-1+t^30-2064705*t^8+156848*t^9+10324398*t^13+214*t^3+197*t^2-54699758*t^16-15137114*t^15-9741*t^4-2972710*t^11+54699758*t^14+202037*t^6+56736*t^7-32425754*t^12-\
7262*t^5+11058754*t^10-197*t^28+214*t^27+9741*t^26-202037*t^24-7262*t^25+10324398*t^17-2972710*t^19+32425754*t^18-11058754*t^20+2064705*t^22+156848*t^21+56736*t^23)
/(1-t+t^32+411*t^29-285*t^30+t^31+6149853*t^8+471319*t^9-58545372*t^13-411*t^3-285*t^2+429447820*t^16+123321948*t^15+18027*t^4+10402780*t^11-335484428*t^14-472275*t
^6-271027*t^7+157353820*t^12+20689*t^5-42303393*t^10+18027*t^28-20689*t^27-472275*t^26+6149853*t^24+271027*t^25-123321948*t^17+58545372*t^19-335484428*t^18+
157353820*t^20-42303393*t^22-10402780*t^21-471319*t^23)]:
end:

#SeqFromRec(S,N): Inputs S=[INI,ope]
#where INI is the list of initial conditions, a ope a list of
#size L, say, and a recurrence operator ope, codes a list of
#size L, finds the first N0 terms of the sequence satisfying
#the recurrence f(n)=ope[1]*f(n-1)+...ope[L]*f(n-L).
#For example, for the first 20 Fibonacci numbers,
#try: SeqFromRec([[1,1],[1,1]],20);
SeqFromRec:=proc(S,N) local gu,L,n,i,INI,ope:
INI:=S[1]:ope:=S[2]:
if not type(INI,list) or not type(ope,list) then
 print(`The first two arguments must be lists `):
 RETURN(FAIL):
fi:
L:=nops(INI):
if nops(ope)<>L then
 print(`The first two arguments must be lists of the same size`):
RETURN(FAIL):
fi:

if not type(N,integer) then
 print(`The third argument must be an integer`, L):
 RETURN(FAIL):
fi:

if N<L then
 RETURN([op(1..N,INI)]):
fi:

gu:=INI:

for n from 1 to N-L do
 gu:=[op(gu),expand(add(gu[nops(gu)-i+1]*ope[i],i=1..L))]:
od:

gu:
end:



# I think that I pretty much wrote this too. It isn't too tricky. Grab enough
# of the terms, form the linear system that the thing would have to satisfy to
# be of degree d, then look for solutions. If you find one, you have a possible
# match. If you don't, then the recurrence can't be linear of that degree.

#GuessRec1(L,d): inputs a sequence L and tries to guess
#a recurrence operator with constant cofficients of order d
#satisfying it. It returns the initial d values and the operator
# as a list. For example try:
#GuessRec1([1,1,1,1,1,1],1);
GuessRec1:=proc(L,d) local eq,var,a,i,n:

if nops(L)<=2*d+2 then
 print(`The list must be of size >=`, 2*d+3 ):
 RETURN(FAIL):
fi:

var:={seq(a[i],i=1..d)}:

eq:={seq(L[n]-add(a[i]*L[n-i],i=1..d),n=d+1..2*d+1)}:

var:=solve(eq,var):

if var=NULL then
 RETURN(FAIL):
else
if expand({seq(L[n]-add(subs(var,a[i])*L[n-i],i=1..d),n=2*d+1..nops(L))})<>{0}  then
  RETURN(FAIL):
else
 RETURN([[op(1..d,L)],[seq(subs(var,a[i]),i=1..d)]]):
fi:
fi:

end:




#GuessRec(L): inputs a sequence L and tries to guess
#a recurrence operator with constant cofficients 
#satisfying it. It returns the initial  values and the operator
# as a list. For example try:
#GuessRec([1,1,1,1,1,1]);
GuessRec:=proc(L) local gu,d:

for d from 1 to trunc(nops(L)/2)-4 do
 gu:=GuessRec1(L,d):
 if gu<>FAIL then
   RETURN(gu):
 fi:
od:
FAIL:

end:

# Gvul is the lowest starting value.
#GuessRecAd(L,Gvul): inputs a sequence L and tries to guess
#a recurrence operator with constant cofficients 
#satisfying it. It returns the initial  values and the operator
# as a list. For example try:
#GuessRec([1,1,1,1,1,1]);
GuessRecAd:=proc(L,Gvul) local gu,d:

for d from 1 to min(Gvul,trunc(nops(L)/2)-4) do
 gu:=GuessRec1(L,d):
 if gu<>FAIL then
   RETURN(gu):
 fi:
od:
FAIL:

end:


#Khibur(S1,S2): inputs two C-finite sequences S1,S2, and  outputs
#their sum, For example, try:
#Khibur([[1,2],[3,4]],[[2,5],[-1,6]]);
Khibur:=proc(S1,S2) local d1,d2,L1,L2,L,i,K:

if nops(S1[1])<>nops(S1[2]) or nops(S2[1])<>nops(S2[2]) then
print(`Bad input`):
RETURN(FAIL):
fi:

d1:=nops(S1[1]): d2:=nops(S2[1]):
K:=2*(d1+d2)+4:

L1:=SeqFromRec(S1,K):
L2:=SeqFromRec(S2,K):

L:=[seq(L1[i]+L2[i],i=1..K)]:

GuessRec(L):

end:



#Kefel(S1,S2): inputs two C-finite sequences S1,S2, and  outputs
#their product, For example, try:
#Kefel([[1,2],[3,4]],[[2,5],[-1,6]]);
Kefel:=proc(S1,S2) local d1,d2,L1,L2,L,i,K:

if nops(S1[1])<>nops(S1[2]) or nops(S2[1])<>nops(S2[2]) then
print(`Bad input`):
RETURN(FAIL):
fi:

d1:=nops(S1[1]): d2:=nops(S2[1]):
K:=2*(d1*d2)+4:

L1:=SeqFromRec(S1,K):
L2:=SeqFromRec(S2,K):

L:=[seq(L1[i]*L2[i],i=1..K)]:

GuessRec(L):

end:



#KhiburSslow(a,b,d1,d2): The recurrence operator
#satisfied by the sum of any solution of
#the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]
#and g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1]
#For example, try:
#KhiburSslow(a,b,2,2);
KhiburSslow:=proc(a,b,d1,d2) local S1,S2,i:

S1:=[[seq(ithprime(5+i),i=1..d1)],[seq(a[i],i=1..d1)]]:

S2:=[[seq(ithprime(9+d1+i),i=1..d2)],[seq(b[i],i=1..d2)]]:

Khibur(S1,S2)[2]:

end:



#KefelSslow(a,b,d1,d2): The recurrence operator
#satisfied by the product of any solution of
#the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]
#and g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1]
#For example, try:
#KefelSslow(a,b,2,2);
KefelSslow:=proc(a,b,d1,d2) local S1,S2,i:

S1:=[[seq(ithprime(5+i),i=1..d1)],[seq(a[i],i=1..d1)]]:

S2:=[[seq(ithprime(9+d1+i),i=1..d2)],[seq(b[i],i=1..d2)]]:

Kefel(S1,S2)[2]:

end:


#KhiburS(a,b,d1,d2): The recurrence operator
#satisfied by the sum of any solution of
#the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]
#and g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1]
#For example, try:
#KhiburS(a,b,2,2);
KhiburS:=proc(a,b,d1,d2) local S1,S2,i,x,y,L,C,var,eq,L1,L2,gu:

S1:=[[seq(x[-d1+i],i=0..d1-1)],[seq(a[i],i=1..d1)]]:

S2:=[[seq(y[-d2+i],i=0..d2-1)],[seq(b[i],i=1..d2)]]:

L1:=SeqFromRec(S1,d1+d2+1):
L2:=SeqFromRec(S2,d1+d2+1):

L:=expand([seq(L1[i]+L2[i],i=1..d1+d2+1)]):

gu:=expand(L[d1+d2+1]-add(C[i]*L[d1+d2+1-i],i=1..d1+d2)):

var:={seq(C[i],i=1..d1+d2)}:
eq:={seq(coeff(gu,x[-i],1),i=1..d1),seq(coeff(gu,y[-i],1),i=1..d2)}:
var:=solve(eq,var):

[seq(subs(var,C[i]),i=1..d1+d2)]:


end:




#KefelS(a,b,d1,d2): The recurrence operator
#satisfied by the product of any solution of
#the recurrence f[n]=a[1]*f[n-1]+...+a[d1]*f[n-d1]
#and g[n]=b[1]*g[n-1]+...+b[d1]*g[n-d1]
#For example, try:
#KefelS(a,b,2,2);
KefelS:=proc(a,b,d1,d2) local S1,S2,i,x,y,L,C,var,eq,L1,L2,gu,j:

S1:=[[seq(x[-d1+i],i=0..d1-1)],[seq(a[i],i=1..d1)]]:

S2:=[[seq(y[-d2+i],i=0..d2-1)],[seq(b[i],i=1..d2)]]:

L1:=SeqFromRec(S1,d1*d2+1):
L2:=SeqFromRec(S2,d1*d2+1):

L:=expand([seq(L1[i]*L2[i],i=1..d1*d2+1)]):

gu:=expand(L[d1*d2+1]-add(C[i]*L[d1*d2+1-i],i=1..d1*d2)):

var:={seq(C[i],i=1..d1*d2)}:
eq:={seq(seq(coeff(coeff(gu,x[-i],1),y[-j],1),j=1..d2),i=1..d1)}:
var:=solve(eq,var):

[seq(subs(var,C[i]),i=1..d1*d2)]:


end:


#ProdIndicatorN(m,n): The profile of multiplicity of
#the (m*n)^2 ratios of a list of m*n numbers to
#be the Cartesian product of a set of m numbers with
#a set of n numbers. For example, try:
#ProdIndicatorN(2,2);
ProdIndicatorN:=proc(m,n) local a,b,i,j,gu,x,lu,mu:
option remember:
gu:=[seq(seq(a[i]*b[j],j=1..n),i=1..m)]:
ProdProfileE(gu):
end:

#ProdIndicator(m,n): The profile of multiplicity of
#the (m*n)^2 ratios of a list of m*n numbers to
#be the Cartesian product of a set of m numbers with
#a set of n numbers. For example, try:
#ProdIndicator(2,2);
ProdIndicator:=proc(m,n) local a,b,i,j,gu,x,lu,mu:
option remember:
gu:=[seq(seq(a[i]*b[j],j=1..n),i=1..m)]:
lu:=[seq(seq(gu[i]/gu[j],j=1..m*n),i=1..m*n)]:
lu:=add(x[lu[i]],i=1..nops(lu)):

mu:=[]:

for i from 1 to nops(lu) do
 if type(op(1,op(i,lu)),integer) then
   mu:=[op(mu),op(1,op(i,lu))]:
 else:
   mu:=[op(mu),1]:
 fi:
od:

sort(mu):

end:


#ProdProfile(L,eps): Given a list of complex numbers finds its
#(approximate profile). For example, try:
#ProdProfile([1,2,3,6],10^(-4));
ProdProfile:=proc(L,eps) local i,j,lu,mu,S,gu:

gu:=[seq(seq(L[i]/L[j],j=1..nops(L)),i=1..nops(L))]:
mu:=[]:
S:={seq(i,i=1..nops(gu))}:

while S<>{} do
 lu:=Krovim(gu,S[1],eps):
 mu:=[op(mu),nops(lu)]:
 S:=S minus lu:
od:
sort(mu):
end:

#IsProd(f,t,m,n,eps): Given a rational function f of t
#decides whether the C-finite sequence of its coeeficients
#is a product of a C-finite sequences of order m and n.
#For example, try:
#IsProd(1/((1-t)*(1-2*t)*(1-3*t)*(1-6*t)),t,2,2);
IsProd:=proc(f,t,m,n,eps) local lu,gu:
lu:=denom(f):

if degree(lu,t)<>m*n then
 print(`The degree of the denominator of`, f, `should be`, m*n):
 RETURN(FAIL):
fi:

gu:=evalf([solve(lu,t)]):


if ProdProfile(gu,eps)=ProdIndicator(m,n) then
 RETURN(true):
else

 RETURN(false,gu):
fi:

end:


 

#CtoR(S,t), the rational function, in t, whose coefficients
#are the members of the C-finite sequence S. For example, try:
#CtoR([[1,1],[1,1]],t);

CtoR:=proc(S,t) local D1,i,N1,L1,f,f1,L:
if not (type(S,list) and  nops(S)=2 and type(S[1],list) and type(S[2],list)
        and nops(S[1])=nops(S[2]) and type(t, symbol) ) then
   print(`Bad input`):
   RETURN(FAIL):
fi:

D1:=1-add(S[2][i]*t^i,i=1..nops(S[2])):
N1:=add(S[1][i]*t^(i-1),i=1..nops(S[1])):
L1:=expand(D1*N1):
L1:=add(coeff(L1,t,i)*t^i,i=0..nops(S[1])-1):
f:=L1/D1:
L:=degree(D1,t)+10:
f1:=taylor(f,t=0,L+1):
if expand([seq(coeff(f1,t,i),i=0..L)])<>expand(SeqFromRec(S,L+1)) then
print([seq(coeff(f1,t,i),i=0..L)],SeqFromRec(S,L+1)):
 RETURN(FAIL):
else
 RETURN(f):
fi:
end:


#RtoC(R,t): Given a rational function R(t) finds the
#C-finite description of its sequence of coefficients
#For example, try:
#RtoC(1/(1-t-t^2),t);
RtoC:=proc(R,t) local S1,S2,D1,R1,d,f1,i,a:
R1:=normal(R):
D1:=denom(R1):
d:=degree(D1,t):
if degree(numer(R1),t)>=d then
print(`The degree of the top must be less than the degree of the bottom`):
 RETURN(FAIL):
fi:

a:=coeff(D1,t,0):
S2:=[seq(-coeff(D1,t,i)/a,i=1..degree(D1,t))]:
f1:=taylor(R,t=0,d+1):
S1:=[seq(coeff(f1,t,i),i=0..d-1)]:
[S1,S2]:
end:



#KefelR(R1,R2,t): inputs rational functions of t,
#R1 and R2,
#and outputs the rational function 
#whose coeffs. is their Hadamard product, for example try:
#KefelR(1/(1-t-t^2),1/(1-t-t^3),t);
KefelR:=proc(R1,R2,t) local S1,S2,S:
S1:=RtoC(R1,t):
S2:=RtoC(R2,t):
S:=Kefel(S1,S2):
normal(CtoR(S,t)):
end:





#Krovim(L,i,eps): Given a list L and an index between
#1 and nops(L), finds the set of indices (including i)
#such that abs(L[j]-L[i])<eps. For example, try:
#Krovim([1.0001,1.01,1.00011],1,1/1000);
Krovim:=proc(L,i,eps) local gu,j,a:
gu:={}:
a:=L[i]:

for j from 1 to nops(L) do
 if abs(L[j]-a)<eps then
  gu:=gu union {j}:
 fi:
od:
gu:
end:

#KrovimE(L,i): Given a list L and an index between
#1 and nops(L), finds the set of indices (including i)
#such that L[i]=L[j] For example, try:
#KrovimE([1,2,1,3],1);
KrovimE:=proc(L,i) local gu,j,a:
gu:={}:
a:=L[i]:

for j from 1 to nops(L) do
 if L[j]=a then
  gu:=gu union {j}:
 fi:
od:
gu:
end:


#ProdProfileE(L): Given a list of numbers (or symbols) finds its
#exact profile. For example, try:
#ProdProfileE([1,2,3,6]);
ProdProfileE:=proc(L) local i,j,lu,mu,S,gu:

gu:=[seq(seq(normal(L[i]/L[j]),j=1..nops(L)),i=1..nops(L))]:
mu:=[]:
S:={seq(i,i=1..nops(gu))}:

while S<>{} do
 lu:=KrovimE(gu,S[1]):
 mu:=[op(mu),nops(lu)]:
 S:=S minus lu:
od:
sort(mu):
end:


#etop(e): If e is a list of numbers (or expressions)
#representing the elementary symmetric functions
#outputs the power-symmetric functions
#For example, try
#etop([1,3,6]);
etop:=proc(e) local k,p,i:

for k from 1 to nops(e) do
 p[k]:=expand(add((-1)^(i-1)*e[i]*p[k-i],i=1..k-1)+(-1)^(k-1)*k*e[k]):
od:

[seq(p[k],k=1..nops(e))]:

end:


#ptoe(p): If p is a list of numbers (or expressions)
#representing the power-sum symmetric functions
#outputs the elementary symmetric functions
#For example, try
#ptoe([1,3,6]);
ptoe:=proc(p) local k,e,i:

for k from 1 to nops(p) do
 e[k]:=expand(add((-1)^(i-1)*p[i]*e[k-i]/k,i=1..k-1)+(-1)^(k-1)*p[k]/k):
od:

[seq(e[k],k=1..nops(p))]:

end:




#IsProdN(f,t,m,n): Given a rational function f of t
#decides whether the C-finite sequence of its coeeficients
#is a product of a C-finite sequences of order m and n.
#For example, try:
#IsProd(1/((1-t)*(1-2*t)*(1-3*t)*(1-6*t)),t,2,2);
IsProdN:=proc(f,t,m,n) local lu,k,e,p,i:
lu:=denom(f):
k:=degree(lu,t):
if k<>m*n then
print(`The denominator of `, f, `should have been of degree`, m*n ):
print(`but it is in fact of degree`, k):
RETURN(FAIL):
fi:

lu:=expand(lu/coeff(lu,t,k)):
e:=[seq((-1)^i*coeff(lu,t,k-i),i=1..k)]:

p:=etop(e):
RETURN(p):
ProdProfileE(p), ProdIndicator(m,n):
end:



#CartProd2(L1,L2): The Cartesian product of lists L1 and L2,
#for example, try:
#CartProd2([2,5],[4,9]);
CartProd2:=proc(L1,L2) local i,j:
[seq(seq(L1[i]*L2[j],i=1..nops(L1)),j=1..nops(L2))]:
end:


#CartProd(L): The Cartesian product of the lists in the list
#of lists L
#for example, try:
#CartProd([[2,5],[4,9],[4,7]]);
CartProd:=proc(L) local L1,L2,i:
if not type(L,list) then
 print(`Bad input`):
 RETURN(FAIL):
fi:

if {seq(type(L[i],list),i=1..nops(L))}<>{true} then
 print(`Bad input`):
 RETURN(FAIL):
fi:

if nops(L)=1 then
RETURN(L[1]):
fi:


L1:=CartProd([op(1..nops(L)-1,L)]):
L2:=L[nops(L)]:
CartProd2(L1,L2):

end:



#ProdIndicatorG(L): The profile of multiplicity of
#the convert(L,`*`)^2 ratios of a list of 
#be the Cartesian product of lists of size given by L
#For example, try:
#ProdIndicatorG([2,2]);
ProdIndicatorG:=proc(L) local a,i,j,gu,x,lu,mu,M:
option remember:
M:=[seq([seq(a[i][j],j=1..L[i])],i=1..nops(L))]:
gu:=CartProd(M):

lu:=[seq(seq(gu[i]/gu[j],j=1..nops(gu)),i=1..nops(gu))]:
lu:=add(x[lu[i]],i=1..nops(lu)):

mu:=[]:

for i from 1 to nops(lu) do
 if type(op(1,op(i,lu)),integer) then
   mu:=[op(mu),op(1,op(i,lu))]:
 else:
   mu:=[op(mu),1]:
 fi:
od:

sort(mu):

end:



#IsProdG(f,t,L,eps): Given a rational function f of t
#decides whether the C-finite sequence of its coeeficients
#is a product of a C-finite sequences of the orders given
#by the list of positive integers L
#For example, try:
#IsProdG(1/((1-t)*(1-2*t)*(1-3*t)*(1-6*t)),t,[2,2],10^(-10));
IsProdG:=proc(f,t,L,eps) local lu,gu,ku1,ku2:
lu:=denom(f):

if degree(lu,t)<>convert(L,`*`) then
 print(`The degree of the denominator of`, f, `should be`, convert(L,`*`)):
 RETURN(FAIL):
fi:

gu:=evalf([solve(lu,t)]):


ku1:=ProdProfile(gu,eps):
ku2:=ProdIndicatorG(L):

if ku1=ku2 then
 RETURN(true):
else

 RETURN(false,[ku1,ku2]):
fi:

end:




#KefelSm(a,L): The recurrence operator
#satisfied by the product of nops(L) recurrences
#the recurrence f[i][n]=a[i,1]*f[n-1]+...+a[i,L[i]]*f[n-L[i]]
#It also returns the set of generic coefficients
#For example, try:
#KefelSm(a,[2,2,2]);
KefelSm:=proc(a,L) local gu,L1,mu,A,B,i,i1,j:

if not (type(a,symbol) and type(L,list) ) then
print(`Bad input`):
 RETURN(FAIL):
fi:

if not {seq(type(L[i],integer),i=1..nops(L))}={true} then
print(`Bad input`):
 RETURN(FAIL):
fi:

if not {seq(evalb(L[i]>=2),i=1..nops(L))}={true} then
print(`Bad input`):
 RETURN(FAIL):
fi:

if nops(L)=1 then
 RETURN([seq(a[1,j],j=1..L[1])]):
fi:

L1:=[op(1..nops(L)-1,L)]:
gu:=expand(KefelSm(a,L1)):


mu:=KefelS(A,B,convert(L1,`*`),L[nops(L)]):
mu:=subs({seq(A[i1]=gu[i1],i1=1..nops(gu)),
seq(B[i1]=a[nops(L),i1],i1=1..L[nops(L)])},mu):
mu:=expand(mu):
mu,{seq(seq(a[i,j],j=1..L[i]),i=1..nops(L))}:

end:

#Factorize(C,L): Given a C-finite recurrence (without the initial conditions)
#tries to find a factorization into C-finite sequences of order
#L[i],i=1..nops(L). 
#It returns the factoriziation, if found, followed by the free parameters
#and FAIL otherwise
#For example, try:
#Factorize([-3,106,-72,-576],[2,2]);
Factorize:=proc(C,L) local gu,a,eq,var,i,var1,j,lu,k,lu1,mu:
if nops(C)<>convert(L,`*`) then
 print(`The number of elements of `, C):
 print(`should be`, convert(L,`*`)):
  ERROR(`Bad input`):
fi:

gu:=KefelSm(a,L):
var:=gu[2]:
gu:=gu[1]:
eq:={seq(gu[i]=C[i],i=1..nops(gu))}:
var1:=[solve(eq,var)]:

if var1=[] then
print(`There is no  solution`):
RETURN(FAIL):
fi:

for k from 1 to nops(var1) do
var1:=var1[k]:
lu:={}:

for i from 1 to nops(var1) do
 if op(1,op(i,var1))=op(2,op(i,var1)) then
   lu:=lu union {op(1,op(i,var1))}:
 fi:
od:



if nops(lu)>0 then
 mu:=[[seq([seq(subs(var1,a[i,j]),j=1..L[i])],i=1..nops(L))]]:
 mu:=subs({seq(lu1=1,lu1 in lu)},mu):
 RETURN(mu):
fi:
od:

FAIL:
end:




#GuessRec1mini(L,d): inputs a sequence L and  guesses
#a recurrence operator with constant cofficients of order d
#satisfying it. It returns the initial d values and the operator
# as a list. For example try:
#GuessRec1mini([1,1,1,1,1,1],1);
GuessRec1mini:=proc(L,d) local eq,var,a,i,n:

if nops(L)<2*d then
 print(`The list must be of size >=`, 2*d ):
 RETURN(FAIL):
fi:

var:={seq(a[i],i=1..d)}:

eq:={seq(L[n]-add(a[i]*L[n-i],i=1..d),n=d+1..nops(L))}:

var:=solve(eq,var):

if var=NULL then
 RETURN(FAIL):
else
 RETURN([[op(1..d,L)],[seq(subs(var,a[i]),i=1..d)]]):
fi:

end:




#ISF(L,k): given a sequence of integers finds all its
#factors of length k. For example try:
#ISF([1,4,9,16,32],3);
ISF:=proc(L,k) local gu,gu1,lu,lu1:
if k<0 then
 RETURN({}):
fi:

if k>nops(L) then
 print(`Bad input`):
 RETURN(FAIL):
fi:

if k=0 then
 RETURN({[]}):
fi:
gu:=ISF(L,k-1):
lu:=divisors(L[k]):

{seq(seq([op(gu1),lu1],lu1 in lu),gu1 in gu)}:


end:



#IISF(L,k): given a sequence of integers finds all its
#increasing sequences
#factors of length k. For example try:
#IISF([1,4,9,16,32],3);
IISF:=proc(L,k) local gu,gu1,lu,lu1,mu:
if k<0 then
 RETURN({}):
fi:

if k>nops(L) then
 print(`Bad input`):
 RETURN(FAIL):
fi:

if k=0 then
 RETURN({[]}):
fi:
lu:=divisors(L[k]):
if k=1 then
 RETURN({seq([lu1],lu1 in lu)}):
fi:

gu:=IISF(L,k-1):


mu:={}:
for gu1 in gu do
 for lu1 in lu do
  if gu1[nops(gu1)]<lu1 then
    mu:=mu union {[op(gu1),lu1]}:
  fi:
 od:
od:
mu:

end:


#FactorizeI1(C,k): Given a C finite integer sequence,C
#tries to factorize it into a C finite sequence of
#order k and order nops(C[1])-k. For example, try:
#FactorizeI1(RtoC(TDB(t)[1],t),2);
FactorizeI1:=proc(C,k) local gu,mu,mu1,lu,gu1,gu11,C1,i1:
gu:=SeqFromRec(C,4*k):
mu:=ISF(gu,2*k) minus {[1$(2*k)]}:

for mu1 in mu do
 lu:=GuessRec1mini(mu1,k):

  if lu<>FAIL then
   gu1:=SeqFromRec(lu,4*k):
     if {seq(type(gu1[i1],integer),i1=1..nops(gu1))}={true} then
      gu11:=[seq(gu[i1]/gu1[i1],i1=1..4*k)]:
         if {seq(type(gu11[i1],integer),i1=1..nops(gu11))}={true} then
          gu:=SeqFromRec(C,10*k):
          gu1:=SeqFromRec(lu,10*k):
          gu11:=[seq(gu[i1]/gu1[i1],i1=1..10*k)]:
           C1:=GuessRec1(gu11,nops(C[2])/k):
           if C1<>FAIL then
              RETURN(lu,C1):
           fi:
         fi:
     fi:
   fi:

  od:
FAIL:
end:




#GuessPR1(C1,C2,d,x,y): Inputs two C-finite sequences and a pos. integer
#d and outputs a polynomial P of degree d such that P(C1[n],C2[n])=0
#for all n. If none found it returns FAIL. For example to solve
#the Hugh Montgomery proposed (but not selected) Putnam problem,
#(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,
#where F(n) is the Fibonacci sequence type
#GuessPR1([[0,1],[1,1]],[[1,1],[1,1]],4,x,y);
GuessPR1:=proc(C1,C2,d,x,y) local eq,var,gu,a,S1,S2,N,var1,pol,n1,mu,v1:
 gu:=GenPol([x,y],a,d):
 var:=gu[2]:
 gu:=gu[1]:
 N:=nops(var):
 S1:=SeqFromRec(C1,N+10):
 S2:=SeqFromRec(C2,N+10):
 eq:={seq(subs({x=S1[n1],y=S2[n1]},gu),n1=1..N+10)}:
 var1:=solve(eq,var):
 pol:=subs(var1,gu):
  if pol=0 then
    RETURN(FAIL):
  fi:

 S1:=SeqFromRec(C1,N+30):
 S2:=SeqFromRec(C2,N+30):

if {seq(normal(subs({x=S1[n1],y=S2[n1]},pol)),n1=N+11..N+30)}<>{0} then
print(`The conjecture`, pol , `did not materialize`):
RETURN(FAIL):
fi:

mu:={}:

for v1 in var1 do
 if op(1,v1)=op(2,v1) then
  mu:=mu union {op(1,v1)}:
 fi:
od:

if nops(mu)>1 then
 print(`There are than one solution`):
fi:

factor(coeff(pol,mu[1],1)):

end:



#GuessPR(C1,C2,d,x,y): Inputs two C-finite sequences and a pos. integer
#d and outputs a polynomial P of degree <=d such that P(C1[n],C2[n])=0
#for all n. If none found it returns FAIL. For example to solve
#the Hugh Montgomery proposed (but not selected) Putnam problem,
#(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,
#where F(n) is the Fibonacci sequence type
#GuessPR([[0,1],[1,1]],[[1,1],[1,1]],6,x,y);
GuessPR:=proc(C1,C2,d,x,y) local gu,d1:
for d1 from 1 to d do
 gu:=GuessPR1(C1,C2,d1,x,y):
 if gu<>FAIL then
  RETURN(gu):
 fi:
od:
FAIL:
end:



#GenPol(Resh,a,deg): Given a list of variables Resh, a symbol a,
#and a non-negative integer deg, outputs a generic polynomial
#of total degree deg,in the variables Resh, indexed by the letter a of
#followed by the set of coeffs. For example, try:
#GenPol([x,y],a,1);
GenPol:=proc(Resh,a,deg)
local var,POL1,i,x,Resh1,POL:

if not type(Resh,list) or nops(Resh)=0 then
 ERROR(`Bad Input`):
fi:

x:=Resh[1]:
if nops(Resh)=1 then
RETURN(convert([seq(a[i]*x^i,i=0..deg)],`+`),{seq(a[i],i=0..deg)}):
fi:

Resh1:=[op(2..nops(Resh),Resh)]:
var:={}:
POL:=0:

for i from 0 to deg do
POL1:=GenPol(Resh1,a[i],deg-i):
var:=var union POL1[2]:
POL:=expand(POL+POL1[1]*x^i):
od:
POL,var:
end:


#GuessNLR1(C1,r,d,x): Inputs a C-finite sequence C and a pos. integer
#d and outputs a polynomial P(x[1],..., x[r]) of degree d such that 
#P(C1[n],C1[n-1],C1[n-r])=0
#for all n. If none found it returns FAIL. For example to solve
#the Hugh Montgomery proposed (but not selected) Putnam problem,
#(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,
#where F(n) is the Fibonacci sequence type
#GuessNLR1([[0,1],[1,1]],1,4,x);
GuessNLR1:=proc(C1,r,d,x) local eq,var,gu,a,S1,N,var1,pol,n1,mu,v1,i,i1:
 gu:=GenPol([seq(x[i],i=0..r)],a,d):
 var:=gu[2]:
 gu:=gu[1]:
 N:=nops(var)+r:
 S1:=SeqFromRec(C1,N+10):
 eq:={seq(subs({seq(x[i1]=S1[n1+i1],i1=0..r)},gu),n1=1..N+10-r)}:
 var1:=solve(eq,var):
 pol:=subs(var1,gu):
  if pol=0 then
    RETURN(FAIL):
  fi:

 S1:=SeqFromRec(C1,N+30):

if 
{seq(
normal(subs({seq(x[i1]=S1[n1+i1-1],i1=0..r)},pol)),n1=N+11..N+30-r-1)}<>{0} then
print(`The conjecture`, pol , `did not materialize`):
RETURN(FAIL):
fi:

mu:={}:

for v1 in var1 do
 if op(1,v1)=op(2,v1) then
  mu:=mu union {op(1,v1)}:
 fi:
od:

if nops(mu)>1 then
 print(`There are than one solution`):
fi:

factor(coeff(pol,mu[1],1)):

end:



#GuessNLR(C1,r,d,x): Inputs a C-finite sequence C and a pos. integer
#d and outputs a polynomial P(x[1],..., x[r]) of degree<= d such that 
#P(C1[n],C1[n-1],C1[n-r])=0
#for all n. If none found it returns FAIL. For example to solve
#the Hugh Montgomery proposed (but not selected) Putnam problem,
#(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,
#where F(n) is the Fibonacci sequence type
#GuessNLR([[0,1],[1,1]],1,4,x);
GuessNLR:=proc(C1,r,d,x) local gu,d1:

for d1 from 1 to d do
gu:=GuessNLR1(C1,r,d1,x):
if gu<>FAIL then
 RETURN(gu):
fi:
od:
FAIL:
end:




#GuessNLRv(C1,r,d): Verbose form of GuessNLR(C1,r,d,x), i.e.
#Inputs a C-finite sequence C and a pos. integer
#d and outputs a polynomial P(x[1],..., x[r]) of degree<= d such that 
#P(C1[n],C1[n-1],C1[n-r])=0
#for all n. If none found it returns FAIL. For example to solve
#the Hugh Montgomery proposed (but not selected) Putnam problem,
#(finding a polynomial of degree 4 such that P(F(n),F(n+1))=0,
#where F(n) is the Fibonacci sequence type
#GuessNLRv([[0,1],[1,1]],1,4,x);
GuessNLRv:=proc(C1,r,d) local gu,x,F,n,i1:
gu:=GuessNLR(C1,r,d,x):

if gu=FAIL then
  RETURN(FAIL):
fi:

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(`Then we have the following identity valid for every non-negative integer n`):
print(subs({seq(x[i1]=F(n+i1),i1=0..r)},gu)=0):
print(``):
print(`Proof: Routine! (since everything in sight is C-finite).`):
end:



#FactorizeI1v(C1,k): Verbose form of FactorizeI1(C1,k), i.e.
FactorizeI1v:=proc(C1,k) local gu,C2,C3,i1,n:
gu:=FactorizeI1(C1,k): 

C2:=gu[1]:
C3:=gu[2]:

if gu=FAIL then
  RETURN(FAIL):
fi:

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):

print(``):
print(`Let G(n) be the sequence defined on non-negative integers by the recurrence`):
print(G(n)=add(C2[2][i1]*G(n-i1),i1=1..nops(C2[2]))):
print(`subject to the initial conditions`):
print(seq(G(i1-1)=C2[1][i1],i1=1..nops(C2[1]))):

print(``):
print(`Let H(n) be the sequence defined by the recurrence`):
print(H(n)=add(C3[2][i1]*H(n-i1),i1=1..nops(C3[2]))):
print(`subject to the initial conditions`):
print(seq(H(i1-1)=C3[1][i1],i1=1..nops(C3[1]))):

print(`Then the following is true for every non-negative integer n`):
print(F(n)=G(n)*H(n)):
print(``):
print(`Proof: Routine! (since everything in sight is C-finite) .`):
end:


#mSect(C1,m,i): Given a C-finite sequence C1(n), finds
#the C-finite description of C(m*n+i) where m is a pos. integer
#and 0<=i<m. For example, try:
#mSect([[0,1],[1,1]],2,0);
mSect:=proc(C1,m,i) local S,n1:
S:=SeqFromRec(C1,2*nops(C1[1])*(m+1)+3*m+10):
S:=[seq(S[i+1+m*n1],n1=0..2*nops(C1[1])+3)]:
GuessRec(S):
end:



#mSectV(C1,m,i): Verbose form of mSect(C1,m,i)
#For example, try:
#mSectV([[0,1],[1,1]],2,0);
mSectV:=proc(C1,m,i) local C2,F,G,n,i1:
C2:=mSect(C1,m,i):

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(``):
print(`Define a new sequence`):
print(G(n)=F(i+m*n)):
print(`then G(n) satisfies the recurrence`):
print(G(n)=add(C2[2][i1]*G(n-i1),i1=1..nops(C2[2]))):
print(`subject to the initial conditions`):
print(seq(G(i1-1)=C2[1][i1],i1=1..nops(C2[1]))):
print(``):
end:



#BT(C1): The Binomial Transform of the C-finite sequence C,
#in other words, the sequence add(binomial(n,i)*C1[i],i=0..n):
#For example, try:
#BT([[0,1],[1,1]]);
BT:=proc(C1) local S,i1,n1:
S:=SeqFromRec(C1,nops(C1[1])*3+5):
S:=[seq(add(S[i1+1]*binomial(n1,i1),i1=0..n1),n1=0..nops(S)-1)]:
GuessRec(S):
end:


#BTv(C1): Verbose form of mSect(C1,m,i)
#For example, try:
#BTv([[0,1],[1,1]]);
BTv:=proc(C1) local C2,F,G,n,i1,i:
C2:=BT(C1):

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(``):
print(`Let G(n) be the sequence that is its Binomial Transform `):
print(`In other words, the sequence defined by`):
print(G(n)=Sum(binomial(n,i)*F(i),i=0..n)):
print(`then G(n) satisfies the recurrence`):
print(G(n)=add(C2[2][i1]*G(n-i1),i1=1..nops(C2[2]))):
print(`subject to the initial conditions`):
print(seq(G(i1-1)=C2[1][i1],i1=1..nops(C2[1]))):
print(``):
end:

#Ux(x): The C-finite description of the Chebychev polynomials of
#the second kind. Try: Ux(x);
Ux:=proc(x):[[1,2*x],[2*x,-1]]:end:

#Tx(x): The C-finite description of the Chebychev polynomials of
#the first kind
Tx:=proc(x):[[1,x],[2*x,-1]]:end:

#Fn(): The C-finite description of the Fibonacci sequence
Fn:=proc():[[0,1],[1,1]]:end:
#Fn(): The C-finite description of fibonacci(n+1)
fn:=proc():[[1,1],[1,1]]:end:
#Ln(): The C-finite description of the Lucas Numbers
Ln:=proc():[[2,1],[1,1]]:end:
#Pn(): The C-finite description of the Pell Numbers
Pn:=proc():[[0,1],[2,1]]:end:

#Trn(): The C-finite description of the Triboanacci numbers Numbers
Trn:=proc():[[1,1,2],[1,1,1]]:end:




#KefelV(C1,C2): Verbose version of Kefel(S1,S2);
#For example try:
#KefelV([[1,2],[3,4]],[[2,5],[-1,6]]);
KefelV:=proc(C1,C2) local F,G,H,n,t,C,i1,S1,S2,S,m:
C:=Kefel(C1,C2):
print(`A proof in the Spirit of Zeilberger of A C-finite Identity`):

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(``):
print(`Equivalently, the ordinary generating function of F(n) w.r.t. to t is`):
print(Sum(F(n)*t^n,n=0..infinity)=CtoR(C1,t)):
print(`Also Let G(n) be the sequence defined by the recurrence`):
print(G(n)=add(C2[2][i1]*G(n-i1),i1=1..nops(C2[2]))):
print(`subject to the initial conditions`):
print(seq(G(i1-1)=C2[1][i1],i1=1..nops(C2[1]))):
print(`Equivalently, the ordinary generating function of G(n) w.r.t. to t is`):
print(Sum(G(n)*t^n,n=0..infinity)=CtoR(C2,t)):
print(`Let H(n)=F(n)G(n), then`):
print(Sum(H(n)*t^n,n=0..infinity)=CtoR(C,t)):
print(`or equivalently, H(n) is the sequence defined by the recurrence`):
print(H(n)=add(C[2][i1]*H(n-i1),i1=1..nops(C[2]))):
print(`subject to the initial conditions`):
print(seq(H(i1-1)=C[1][i1],i1=1..nops(C[1]))):

print(``):
m:=nops(C1[1])*nops(C2[1])+nops(C[1]):
print(`Proof: F(n) is a C-finite sequence of order`, nops(C1[2])):
print(`and G(n) is a C-finite sequence of order`, nops(C2[2])):
print(`Hence their product is a C-finite sequence of order`):
print( nops(C1[1])*nops(C2[1])):
print(`Since the claimed sequence for their product is a C-finite sequence`):
print(`of order`, nops(C[1]),` it suffices to check the equality for`):
print(`the first`, m, ` values `):

S1:=expand(SeqFromRec(C1,m)):
S2:=expand(SeqFromRec(C2,m)):
S:=expand(SeqFromRec(C,m)):

print(`The first`, m , `members of the sequence F(n) are`):
print(S1):

print(`The first`, m,  `members of the sequence G(n) are`):
print(S2):
print(`The first`, m, `members of the proposed   product are`):
print(S):

for i1 from 1 to m do
if expand(S1[i1]*S2[i1]-S[i1])<>0 then
 print(`Oops something is wrong `):
 RETURN(FAIL):
else
 print(S1[i1], `times `, S2[i1], `is indeed `, S[i1]):
fi:
od:
print(`QED!`):

end:




#CHseq(L,n,j,K): inputs  a list of pairs [C,e], where C is a C-finite
#sequence, and e is an affine-linear expression in the symbols
#n and j, of the form a*n+b*j+c, where a and a+b are non-negative
#integers. Outputs the first K terms of
#Sum(L[1][1](L[1][2])*L[2][1](L[2][2]),j=0..n-1)
#using up to K values. For example, try:
#CHseq([[Fn(),j],[Fn(),j],[Fn(),2*n-j],20);
CHseq:=proc(L,n,j,K) local i1,S,j1,n1,Max1:

Max1:=0:
for i1 from 1 to nops(L) do
 for n1 from 0 to K do
   for j1 from 0 to n1-1 do
   Max1:=max(Max1,subs({j=j1,n=n1},L[i1][2])):
    if subs({j=j1,n=n1},L[i1][2])<0 then
      RETURN(FAIL):
    fi:
   od:
 od:
od:

for i1 from 1 to nops(L) do
S[i1]:=SeqFromRec(L[i1][1],Max1+1):
od:

[
seq(
add(mul(S[i1][subs({j=j1,n=n1},L[i1][2])+1],i1=1..nops(L)),j1=0..n1-1),
n1=0..K)]:

end:


#CH(L,n,j): 
#Find a sequence of the style of Curtis Greene and
#Herb Wilf ("Closed form summation of C-finite sequences",
#Trans. AMS)
#inputs  a list of pairs [C,e], where C is a C-finite
#sequence, and e is an affine-linear expression in the symbols
#n and j, of the form a*n+b*j+c, where a and a+b are non-negative
#integers. Outputs the C-finite description of
#Sum(L[1][1](L[1][2])*L[2][1](L[2][2]),j=0..n-1)
#using up to K values. For example, try:
#CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]);
CH:=proc(L,n,j) local K, C:
C:=GuessRec(CHseq(L,n,j,20)):

if C<>FAIL then
 RETURN(C):
fi:
for K from 25 by 5 while C=FAIL do
C:=GuessRec(CHseq(L,n,j,K)):
od:

C:

end:





#GuessPol1(L,d,n): guesses a polynomial of degree d in n for
# the list L, such that P(i)=L[i+1] for i=1..nops(L)
#For example, try: 
#GuessPol1([seq(i,i=1..10)],1,n);
GuessPol1:=proc(L,d,n) local P,i,a,eq,var:
if d>nops(L)-2 then
 ERROR(`the list is too small`):
fi:

P:=add(a[i]*n^i,i=0..d):
var:={seq(a[i],i=0..d)}:
eq:={seq(subs(n=i,P)-L[i],i=1..nops(L))}:

var:=solve(eq,var):

if var=NULL then
 RETURN(FAIL):
fi:

subs(var,P):

end:

#GuessPol(L,n): guesses a polynomial of degree d in n for
# the list L, such that P(i)=L[i] for i=1..nops(L) for example, try: 
#GuessPol([seq(i,i=1..10)],n);
GuessPol:=proc(L,n) local d,gu:

for d from 0 to nops(L)-2 do
 gu:=GuessPol1(L,d,n):
 if gu<>FAIL then
    RETURN(gu):
 fi:
od:

FAIL:

end:


#####New stuff



#GuessCH1(C1,r,d,x,C2): Inputs a C-finite sequence C1 and a pos. integer
#d and outputs a polynomial P(x[1],..., x[r]) of degree d such that 
#P(C1[n],C1[n-1],C1[n-r])=C2[n]
#for all n. If none found it returns FAIL. For example to reproduce
#Eq. (12) of the Greene-Wilf article
# type
#GuessCH1([[0,1],[1,1]],1,3,x,CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j));
GuessCH1:=proc(C1,r,d,x,C2) local eq,var,gu,a,S1,N,var1,pol,n1,mu,v1,i,i1,S2:
 gu:=GenPol([seq(x[i],i=0..r)],a,d):
 var:=gu[2]:
 gu:=gu[1]:
 N:=nops(var)+r:
 S1:=SeqFromRec(C1,N+11):
 S2:=SeqFromRec(C2,N+11):
 eq:={seq(subs({seq(x[i1]=S1[n1+i1+1],i1=0..r)},gu)-S2[n1+1],n1=1..N+10-r)}:
 var1:=solve(eq,var):
 pol:=subs(var1,gu):
  if pol=0 then
    RETURN(FAIL):
  fi:


 S1:=SeqFromRec(C1,N+30):
 S2:=SeqFromRec(C2,N+30):

if 
{seq(
normal(subs({seq(x[i1]=S1[n1+i1+1],i1=0..r)},pol)-S2[n1+1]),n1=N+11..N+30-r-2)
}<>{0} then
RETURN(FAIL):
fi:

mu:={}:

for v1 in var1 do
 if op(1,v1)=op(2,v1) then
  mu:=mu union {op(1,v1)}:
 fi:
od:

if mu<>{} then
factor(subs({seq(mu1=0,mu1 in mu)},pol)):
else factor(pol):
fi:

end:




#FindCHg(C1,r,d, x,C2): Inputs a C-finite sequence C1 and a pos. integer
#d and outputs a polynomial P(x[1],..., x[r]) of degree<= d such that 
#P(C1[n],C1[n-1],C1[n-r])=C2[n]
#for all n. If none found it returns FAIL. For example to reproduce
#Eq. (12) of the Greene-Wilf article.
#If there is none of degree <=d then it returns FAIL
# type
#FindCHg([[0,1],[1,1]],2,4,x,CH([[Fn(),j],[Fn(),j],[Fn(),2*n-j]],n,j));
FindCHg:=proc(C1,r,d,x,C2) local gu,d1:

for d1 from 1 to d do
gu:=GuessCH1(C1,r,d1,x,C2):
if gu<>FAIL then
 RETURN(gu):
fi:
od:
FAIL:
end:

###end new stuff



#FindCH(C1,r,d, x,L,n,j): Inputs a C-finite sequence C1 and pos. integers
#r and d, and a letter x, and a list L of affine-linear expressions
#in n,j, and
#and outputs a polynomial P(x[1],..., x[r]) of degree<= d such that 
#P(C1[n],C1[n-1],C1[n-r])=C2[n]
#where C2(n) is the sum(C1[L[1]]*C1[L[2]*...C1[L[nops(L)],j=0..n-1):
#If none found it returns FAIL. For example to reproduce
#Eq. (12) of the Greene-Wilf article.
#If there is none of degree <=d then it returns FAIL
# type
#FindCH(Fn(),1,3,x,[j,j,2*n-j],n,j);
FindCH:=proc(C1,r,d,x,L,n,j) local C2,i1,L1:
L1:=[seq([C1,L[i1]],i1=1..nops(L))]:
FindCHg([[0,1],[1,1]],r,d,x,CH(L1,n,j)):
end:



#FindCHv(C1,r,d, L,n,j): verbose version of FindCH(C1,r,d,x, L,n,j): 
# type
#FindCHv(Fn(),1,3,[j,j,2*n-j],n,j);
FindCHv:=proc(C1,r,d,L,n,j) local gu,x,F,i1,t:
gu:=FindCH(C1,r,d,x,L,n,j):

if gu=FAIL then
  RETURN(gu):
fi:

print(`Theorem: Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(``):
print(`Equivalently, the ordinary generating function of F(n) w.r.t. to t is`):
print(Sum(F(n)*t^n,n=0..infinity)=CtoR(C1,t)):
print(`Also Let G(n) be the sequence defined by `):
print(G(n)=Sum(mul(F(L[i1]),i1=1..nops(L)),j=0..n-1)):
print(`then we have the following closed-form expression for G(n)`):
print(G(n)=subs({seq(x[i1]=F(n+i1),i1=0..r)},gu)):
print(``):
print(`Proof: Both sides are C-finite, so it is enough to check the`):
print(`first few terms, and this is left to the dear reader.`):

end:


#FindCHvTN(C1,r,d, L,n,j,mispar): verbose version of FindCH(C1,r,d,x, L,n,j): 
# type
#FindCHvTN(Fn(),1,3,[j,j,2*n-j],n,j,8);
FindCHvTN:=proc(C1,r,d,L,n,j,mispar) local gu,x,F,i1,t:
gu:=FindCH(C1,r,d,x,L,n,j):

if gu=FAIL then
  RETURN(gu):
fi:

print(`Theorem Number`, mispar):
print( `Let F(n) be the sequence defined by the recurrence`):
print(F(n)=add(C1[2][i1]*F(n-i1),i1=1..nops(C1[2]))):
print(`subject to the initial conditions`):
print(seq(F(i1-1)=C1[1][i1],i1=1..nops(C1[1]))):
print(``):
print(`Equivalently, the ordinary generating function of F(n) w.r.t. to t is`):
print(Sum(F(n)*t^n,n=0..infinity)=CtoR(C1,t)):
print(`Also Let G(n) be the sequence defined by `):
print(G(n)=Sum(mul(F(L[i1]),i1=1..nops(L)),j=0..n-1)):
print(`then we have the following closed-form expression for G(n)`):
print(G(n)=subs({seq(x[i1]=F(n+i1),i1=0..r)},gu)):
print(``):
print(`Proof: Both sides are C-finite, so it is enough to check the`):
print(`first few terms, and this is left to the dear reader.`):

end:




#Curt1(j,n,M): outputs 
#all {a1n+b1j} with 0<=a1<=M and 0<=a1+b1<=M
#Type: Curt1(j,n,4);
Curt1:=proc(j,n,M)  local a1,b1:
[seq(seq(a1*n+b1*j, b1=-a1..-1),a1=0..M),
seq(seq(a1*n+b1*j, b1=1..M-a1 ),a1=0..M)]:
end:

#IV1(m,K): all the lists of weakly-increasing positive integers of size m
#with largest element K. Do IV(3,5);
IV1:=proc(m,K) local gu,i,gu1,gu11:
option remember:
if m=1 then
 RETURN([seq([i],i=1..K)]):
fi:
gu:={}:
for i from 1 to K do
gu1:=IV1(m-1,i):
gu:=gu union {seq([op(gu11),K],gu11 in gu1)}:
od:
gu:
end:



#IV(m,K): all the lists of weakly-increasing positive integers of size m
#with largest element <=K. Do IV(3,5);
IV:=proc(m,K) local i:
{seq(op(IV1(m,i)),i=1..K)}:
end:


Khevre:=proc(m,j,n,M) local gu,mu,mu1,i1:
gu:=Curt1(j,n,M):
mu:=IV(m-1,nops(gu)):
[seq([j,seq(gu[mu1[i1]],i1=1..m-1)],mu1 in mu)]:
end:

#SeferGW(m,d,M): outputs a book of Fibonacci identities with
#right hand sides of degree at most d in F(n),F(n+1)
#for Greene-Wilf type sums where the summand is
#of the form F(j)*F(a1n+b1j)*... with m terms, for example, try:
#SeferGW(2,5,2);
SeferGW:=proc(m,d,M) local n,j,gu,mu,mu1,x,C1,co:
co:=0:
print(`A Book of Definite Summation Fibonacci Identities`):
print(`in the style of Curtis Greene and Herbert Wilf`):
print(``):
print(`By Shalosh B. Ekhad `):

C1:=Fn():
mu:=Khevre(m,j,n,M):
for mu1 in mu do
gu:=FindCH(C1,1,d,x,mu1,n,j):

 if gu<>FAIL then
co:=co+1:
print(`------------------------------------------------`):
   FindCHvTN(C1,1,d,mu1,n,j,co):
 fi:
od:

print(`------------------------------------------------`):
print(`I hope that you enjoyed, dear reader, these`, co, `beautiful and`):
print(`deep theorems. `):
end:
###end from Cfinite.txt

#PS(x,y,r,K): the first K terms of the Pisot-like sequence defined by
#a[1]=x, a[2]=y, a[n]=[ a(n-1)^2/a(n-2) +r] (where [t] means integer-part). For example, for the first 20 terms op the
#Fibonacci sequence
#(without 0,1,1) type:
#PS(2,3,1/2,20);
PS:=proc(x,y,r,K) local gu,i:

gu:=[x,y]:

for i from 1 to K-2 do
 gu:=[op(gu), trunc(gu[nops(gu)]^2/gu[nops(gu)-1]+r)]:
od:

gu:

end:

#Disc2(L): inputs a sequence of numbers L and outputs a sequence shorter by 2 whose n-th entry is
#L[n+1]^2-L[n+2]*L[n]. Try:
#Disc2(SeqFromRec([[1,1],[1,1]],22),20);
Disc2:=proc(L) local i:
[seq(-L[i+2]*L[i]+L[i+1]^2,i=1..nops(L)-2)]:
end:

#Disc2r(L): inputs a sequence of numbers L and outputs a sequence shorter by 2 whose n-th entry is
#(L[n+1]^2-L[n+2]*L[n])/L[n]. For the sequence L to be a r-Pisot sequence it has to be >=-r and <1-r. Try
#Disc2r(SeqFromRec([[1,1],[1,1]],22),20);
Disc2r:=proc(L) local i:
[seq( (-L[i+2]*L[i]+L[i+1]^2)/L[i],i=1..nops(L)-2)]:
end:

#DO2(C): inputs a C-finite sequence C, and outputs the C-finite description of its Disc-2-sequence i.e.
#The linear recurrence equation satisfied by the sequence  L[n+1]^2-L[n+2]*L[n]. Try:
#DO2([[1,1],[1,1]]);
DO2:=proc(C) local d1,L,M:
d1:=nops(C[1]):
L:=SeqFromRec(C,6*d1+10):
M:=Disc2(L):
GuessRec(M):
end:



#PisInd1(C,r,K): inputs a C-finite description C, describing a sequence, let's call it L, and a positive integer K, outputs
#the the sequence, from n=1 to K, of the quantity (L[n+1]^2-L[n+2]*L[n])/L[n]+r for n from 1 to K. 
#They should all be >=0 and <1 respectively. Try:
#PisInd1([[2,3],[1,1]],1/2,100). 
#PisInd1([[4,7,12],[2,-1,1]],1/2,100). 
PisInd1:=proc(C,r,K) local L,n:
L:=SeqFromRec(C,K+2):
[seq((L[n+1]^2-L[n+2]*L[n])/L[n]+r,n=1..K)]:
end:


# PisInd: Pisot index.
#PisInd(C,r,K1,K2): inputs a C-finite description C, describing a sequence, let's call it L, and  positive integers K1, K2, outputs
#the min and max from n=K1 to K2 of the quantity (L[n+1]^2-L[n+2]*L[n])/L[n]+r 
#They should be >=0 and <1 respectively. Try:
#PisInd([[2,3],[1,1]],1/2,100,200). 
#PisInd([[4,7,12],[2,-1,1]],1/2,100,200). 
PisInd:=proc(C,r,K1,K2) local L,n,lu:
L:=SeqFromRec(C,K2+2):
lu:=[seq((L[n+1]^2-L[n+2]*L[n])/L[n]+r,n=K1..K2)]:
[min(op(lu)),max(op(lu))]:
end:

# 2017-09-12
# Robert's commentary on PtoR:

# This function takes positive integers x, y, and a rational r in [0, 1]. It
# then does some computations to guess whether or not the Pisot sequence E_r(x,
# y) is a C-finite sequence. Success or failure to live up to this expectation
# is signaled through various return values. These values are commented on
# below.  According to the original comments, the final branch of this
# procedure indicates that the Pisot sequence _is_ a C-finite sequence, namely
# the one that is returned.


#PtoR(x,y,r,L,K): inputs positive integers x,y, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]
#and a positive integer K. First guesses whether there seems to be a linear recurrence equation with constant coefficients
#(given by its C-finite representation) of order L. If not, it returns FAIL.
#If it persists to K terms (make K large), then it returns
#the C-finite description, followed by its Pisot-index-pair for these K terms, and the ones for the second half.
#(so it returns a list of length 3 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex, PisotIndicator] ]
#If the initial conjecture fails, then it returns it returns [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:
#PtoR(2,3,1/2,10,1000);
#PtoR(5,17,1/2,10,1000);
#PtoR(5,21,0,5,1000);
#PtoR(4,7,1/2,10,5000);
#PtoR(10,219,1/2,10,2000);
#PtoR(8,55,1,10,12000);
PtoR:=proc(x,y,r,L,K) local C,gu1,gu2,i,ku1,ku2:

if r<=0 or r>=1 then
 print(`r must be strictly between 0 and 1 for this procedure`):
 RETURN(FAIL):
fi:

# A possible C-finite sequence representation is guessed using linear algebra.
# If no such guess can be made, either we didn't look for enough degrees, or
# the Pisot sequence isn't a C-finite sequence. Return FAIL if no guess was
# made.
C:=GuessRecAd(PS(x,y,r,2*L+10),L):
if C=FAIL then
 RETURN(FAIL):
fi:

# Compute K terms of the conjectured sequence and the Pisot sequence. If they
# don't match, then return
#   [FAIL, conjectured sequence, index of first difference,
#       absolute value of second largest root for conjectured sequence].
gu1:=SeqFromRec(C,K):

gu2:=PS(x,y,r,K):

if gu1<>gu2 then
 for i from 1 to K while gu1[i]=gu2[i] do od:
 RETURN([FAIL,C,i,Pis(C)]):
fi:

# `PisInd` computes the Pisot-index of the sequence in the given range. The
# Pisot-index over the given indicies is the minimum and maximum values of the
# fractional expression on page 4, then adding r. If this quantity is in [0,
# 1), then the conjecture holds for all checked values because the floor
# relation holds. `ku1` is a list: the first element is the min, the second is
# the max.
ku1:=evalf(PisInd(C,r,1,K),10):


# By the very first check in this procedure, 0 < r < 1, so this doesn't make a
# lot of sense.
if (r=0 and ku1[1]<1/10000) or (r=1 and ku1[2]>1-1/10000) then
 RETURN( [FAIL, C,ku1]):
fi:

# This is the Pisot-index for the "second-half" of 1, 2, ..., K. I'm not sure
# why we care about this, but Zeilberger thought it was important to see.
ku2:=evalf(PisInd(C,r,trunc(K/2),K),10):

if C[2]=[2,-1] then

 [C, ku1, ku2, LINEAR ]:

# `Pis(C)` computes the absolute value of the second-largest root for the
# C-finite sequence. (It's, uh, pretty computationally expensive. I don't know
# why we don't save it. Maybe Maple uses memeoization?)

# If the root is larger than 1, then we _know_ that the conjecture can't hold.
# In fact, using logarithms we can discover exactly when it will break down. (I
# follow the gist of their argument, but am not fluent enough to derive it
# completely without breaking open an aysmptotic analysis reference and some
# paper.)

# If the conjecture doesn't hold, then return
#   [FAIL, C-finite, spot of first difference,
#       absolute value of second largest root.]
elif Pis(C)>1 then

 [FAIL,C, FirstDev(x,y,r,C), Pis(C) ]:

else

# And now we are pretty sure that it works. The second-largest root has
# magnitude less than one, so the requisite ratio goes to zero exponentially
# quickly after a certain point. (This point could be made explicit if need
# be.)

# C-finite representation, Pisot-index up to K, Pisot-index for [K/2, K], and
# absolute value of second largest root of the characteristic equation of the
# linear recurrence.
 [C, ku1, ku2, Pis(C) ]:
fi:

end:



#Binet(C,n): The Binet solution to the C-finite sequence C in terms of n. Try:
#Binet([[1,1],[1,1]],n);
Binet:=proc(C,n) local eq,var,F,i,lu,x,A:
lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])),x)]:

if nops(lu)<>nops(C[1]) then
 RETURN(FAIL):
fi:

F:=add(A[i]*lu[i]^n,i=1..nops(lu)):
var:={seq(A[i],i=1..nops(lu))}:
eq:={seq(subs(n=i-1,F)-C[1][i],i=1..nops(C[1]))}:
var:=solve(eq,var):
subs(var,F):


end:

#Asy(C): Inputs a C-finite sequence and outputs the pair [A,alpha] such that the n-th term is asymptotic to A*alpha^n
#Try:
#Asy([[1,1],[1,1]]);
Asy:=proc(C) local eq,var,F,i,lu,x,A,aluf,n:
lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])))]:

if nops(lu)<>nops(C[1]) then
 RETURN(FAIL):
fi:

aluf:=1:

for i from 2 to nops(lu) do
 if abs(evalf(lu[i]))>abs(evalf(lu[aluf])) then
  aluf:=i:
 fi:
od:

lu:=[lu[aluf],op(1..aluf-1,lu),op(aluf+1..nops(lu),lu)]:

if max(seq(evalf(abs(lu[i])),i=2..nops(lu)))=evalf(abs(lu[1])) then
RETURN(FAIL):
fi:



F:=add(A[i]*lu[i]^n,i=1..nops(lu)):
var:={seq(A[i],i=1..nops(lu))}:
eq:={seq(subs(n=i-1,F)-C[1][i],i=1..nops(C[1]))}:
var:=solve(eq,var):
[subs(var,A[1]),lu[1]]:


end:

#UB(C): Inputs a C-finite sequence and outputs the largest absolute value of the roots of the characteristic equation
#Try:
#UB([[1,1],[1,1]]);
UB:=proc(C) local x,lu,i:
lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])))]:

if nops(lu)<>nops(C[1]) then
 RETURN(FAIL):
fi:

max(seq(abs(evalf(lu[i])),i=1..nops(lu))):


end:


#Tikva(C,N):  the first N terms of the sequence (L[n+1]^2-L[n+2]*L[n])/L[n] where L[n] is the C-finite sequence given by
#C. For it to be an r-Pisot sequence EVERY term must be >=-r and <1-r. Try:
#Tikva([[5,21,88,368,1538],[4,1,-1,0,-1]],27);
Tikva:=proc(C,N) Disc2r(SeqFromRec(C,N+2)):end:

#Sho(C): all the roots of the characteristic equation of the recurrence F. Try:
#Sho([[1,1],[1,1]]);
#Sho([[10,219,4796,105030],[22,-3,18,-11]]);
Sho:=proc(C) local x,P,C2,i:
C2:=C[2]:
P:=x^nops(C2)-add(C2[i]*x^(nops(C2)-i),i=1..nops(C2)):
[solve(P,x)]:
end:



# How does this procedure handle only one root that is not equal to one? I
# believe that it would fail, but I have not tried it.

# Sketch of result for single root:
#   The conjecture is that a_n = b_n, where b_n = p^n x for some real p.
#   This is true iff
#       p^n x = floor(p^n x + r),
#   which holds iff 0 <= r < 1.

# That is, if it looks like the sequence is a trivial geometric sequence, then
# it probably is, as long as r < 1.

# More formally: If y / x = p, p^n x is an integer for all nonnegative integers
# n, and 0 <= r < 1, then E_r(x, y) is given by a_n = p^n x.

#Pis(C): Inputs a C-finite sequence and outputs the absolute value of the second-largest root
#It is a Pisot number if it is less than 1.
#Fis([[1,1],[1,1]]);
#Pis([[10,219,4796,105030],[22,-3,18,-11]]);
Pis:=proc(C) local x,lu,i,aluf,mu:
lu:=[solve(x^nops(C[2])-add(C[2][i]*x^(nops(C[2])-i),i=1..nops(C[2])))]:

if nops(lu)<>nops(C[1]) then
 RETURN(FAIL):
fi:

if member(1,lu) then
lu:=convert({op(lu)} minus {1},list):
if lu=[] then
 RETURN(FAIL):
fi:
fi:

aluf:=1:

for i from 2 to nops(lu) do
 if abs(evalf(lu[i]))>abs(evalf(lu[aluf])) then
  aluf:=i:
 fi:
od:

mu:=evalf([op(1..aluf-1,lu),op(aluf+1..nops(lu),lu)]):

max(seq(abs(mu[i]),i=1..nops(mu))):


end:



#PtoRv(x,y,r,L,N,K): Verbose version of PtoR(x,y,r,L,K) (q.v.)
# and also printing out the first N terms of the Pisot sequence PS(x,y,r)
#Outputs a theorem with a sketch of a proof. Try
#PtoRv(2,3,1/2,10,1000);
#PtoRv(5,17,1/2,10,1000);
#PtoRv(5,21,0,5,1000);
#PtoRv(4,7,1/2,10,5000);
#PtoRv(10,219,1/2,10,2000);
#PtoRv(8,55,1,10,12000);
PtoRv:=proc(x,y,r,L,N,K) local gu,mu,a,n,C,lu1,lu2,i,beta,ku,i1,t,aluf,lu,A,ru:

if r<=0 or r>=1 then
 print(`r must be strictly between 0 and 1 for this procedure`):
 RETURN(FAIL):
fi:
gu:=PtoR(x,y,r,L,K):


# I think we only have to check a few terms for trivial linearity. If it really
# is linear, then that proves it. (I think.)
if gu[-1]=LINEAR then
print(`The Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`whose `, N, `first terms `):
print(``):
print(PS(x,y,r,N)):
print(``):
print(`is a trivial linear sequence`):
RETURN(gu):
fi:

# We couldn't find a linear recurrence. Like the output says, this proves that
# the Pisot sequence can't satisfy any linear recurrence relation with constant
# coefficients for any of the degrees that we checked. (This is due to some
# linear algebra work. If it did, then the sequence would satisfy a certain set
# of linear equations, which it didn't.)
if gu=FAIL then
print(`The Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`whose `, N, `first terms `):
print(``):
print(PS(x,y,r,N)):
print(``):
print(`does not satisfy`):
print(`a linear recurrence equation with constant coefficients of order<=`, L):
RETURN(gu):
fi:

mu:=PS(x,y,r,N):

# We found a linear recurrence, but the Pisot indicator (absolute value of
# second largest root) is greater than 1, so the linear recurrence can't be a
# Pisot sequence. From the paper, we can use logarithms to determine exactly
# when it breaks down. (This is done later.)
if gu<>FAIL and gu[1]=FAIL and nops(gu)=4 then
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
print(`Fact: Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`At first sight it seems to satisfy the following linear recurrence: `):
C:=gu[2]:
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):

print(`Alas, it breaks down at the`, gu[3], `-th term. `):
print(``):
lu2:=SeqFromRec(C,gu[3])[gu[3]]:
lu1:=PS(x,y,r,gu[3])[gu[3]]:
print(a(gu[3]), `equals `, lu1 , `while the corresponding term  for the solution of the recurrence is `, lu2):
print(``):
print(`So the difference of the former from the latter is`, lu1-lu2):
print(` `):
print(`Note that the Pisot Indicator is`, gu[4] ):

if gu[4]>1 and r>0 and r<1 then
print(`Since it is larger than 1, it is not at all suprising, that it does not go on for ever.`):
fi:
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
RETURN(gu):
fi:


# As far as I can tell, PtoR never actually returns a list of length 3, so this
# branch should never execute.
if gu<>FAIL and gu[1]=FAIL and nops(gu)=3 then
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
print(`Fact: Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`It  satisfies the following linear recurrence for n from 1 to `, K):
C:=gu[2]:
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):

# It is not at all clear that this should run...
print(`but it is not at all clear that this goes for ever, and it probably does not. `):

if r=1 then
print(`Note that the maximum of the Pisot index is`, gu[3][2] ):
fi:

if r=0 then
print(`Note that the minimum of the Pisot index is`, gu[3][1] ):
fi:

print(`You are welcome to try to make`, K , ` bigger. `):
RETURN(gu):
fi:





if r=0 and gu[2][1]<1/100 or r=1 and gu[2][2]>99/100 then
 RETURN(FAIL):
fi:

# We are in the "success" stages.
if gu<>FAIL and gu[1]<>FAIL and r<1 and r>0 then


# This is the Pisot indicator (absolute value of second-largest root.)
beta:=gu[4]:

# This is the conjectured C-finite sequence.
C:=gu[1]:



print(``):
print(`-------------------------------------------------------------------------`):
print(``):

# This should have already been caught, but perhaps not.
if beta>1 then
print(`Fact: Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`At first sight it seems to satisfy the following linear recurrence: `):
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):
print(`However, this can't go on for ever, since the Pisot indicator`, beta, `is >= 1. `):

else

print(`Theorem: Let a(n) be the Pisot Sequence with parameter`, r, `and with initial conditions`, a(1)=x, a(2)=y , ` i.e. for n>1`):
print(``):
print(a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(``):
print(`[Here are the first`, N, `terms `):
print(PS(x,y,r,N)):
print(`]`):
print(`The sequence a(n) satisfies, for n>=`, nops(C[2])+1, ` the linear recurrence equation with constant coefficient `):
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2]))):
print(``):
print(`with initial conditions`, seq(a(i)=mu[i],i=1..nops(C[1]))  , ` . ` ):

print(``):
print(`Proof: Let b(n) be the solution of the above recurrence, i.e. the sequence defined by `):
print(``):
print(b(n)=add(C[2][i]*b(n-i),i=1..nops(C[2]))):
print(``):
print(`with initial conditions`, seq(b(i)=mu[i],i=1..nops(C[1])) , ` . `):
print(``):
print(`We have to prove that a(n)=b(n), but it is easier to prove the equivalent statement that b(n)=a(n). `):
print(``):
print(`It is readily checked that b(n)=a(n) for n from 1 to`, nops(C[1]), ` . `):
print(``):
print(`It remains to prove that b(n) satisfies the Pisot-recurrence `):
print(``):
print(b(n+2)=trunc(r+b(n+1)^2/b(n)) , ` . `)  :
print(``):
print(`But this is equivalent to the INEQALITIES:`):
print(``):
print(-r <=b(n+1)^2/b(n) -b(n+2) ,  b(n+1)^2/b(n)-b(n+2) <1-r ):
print(``):
print(`In other words`):
print(``):
print(-r <=normal(b(n+1)^2/b(n) -b(n+2)) ,  normal(b(n+1)^2/b(n)-b(n+2)) <=1-r ):
print(``):
print(`The characteristic equation in t,  of the recurrence satisfied by b(n) is`):

ku:=t^nops(C[2])- add(C[2][i1]*t^(nops(C[2])-i1),i1=1..nops(C[2])):
print(``):
print(sort(ku)=0):

print(``):

print(`whose roots are`):
print(``):
lu:=[solve(ku)]:
print(lu):

if member(1, {op(lu)}) then
 print(`Since 1 is a root, let's remove it, and the remaining roots are`):
 lu:={op(lu)} minus {1}:
 lu:=convert(lu,list):
fi:

print(`In floating-point`):
print(evalf(lu)):
aluf:=1:

for i from 2 to nops(lu) do
 if abs(evalf(lu[i]))>abs(evalf(lu[aluf])) then
  aluf:=i:
 fi:
od:

print(`The largest root is`, evalf(lu[aluf])):
print(``):
mu:=evalf([op(1..aluf-1,lu),op(aluf+1..nops(lu),lu)]):
print(`and the remaining roots are `):
print(``):
print(mu):
print(``):
print(`whose absolute values are`):
print(``):
print([seq(abs(mu[i1]),i1=1..nops(mu))]):
print(``):

if beta<1 then
print(`so the largest absolute value is`, gu[4]):
print(`that is less than 1,`):
else
print(`so the largest absolute value is`, gu[4]):
print(`that equals 1.`):

fi:

print(``):
print(`It follows that the sequence `):
print(``):
print(c(n)=normal(b(n+1)^2/b(n) -b(n+2))):
print(`satisfies the inequality `):
print(abs(c(n))<=A*gu[4]^n ):
print(``):
if beta<1 then
# c(n) is as in the paper, so this tells us that the numerator of the requisite
# fraction goes to zero. (And knowing that beta < 1 actually tells us that the
# whole thing goes to zero.)
print(`for some fixed constant, A, (independent of n), that can be easily determined, if desired. In particular the sequence c(n) tends to 0`):
else
# This else should never run. We only use beta < 1 this far down. (It's been
# checked like 10 times.)
print(`for some fixed constant, A, (independent of n), that can be easily determined, if desired. In particular the sequence c(n) is bounded.`):
fi:

print(`we have to show that for all n, c(n) is between `, -r, `and `, 1-r ):
print(`For the first`,  3*N-2,  `terms the sequence c(n) is`):
print(``):
mu:=PS(x,y,r,3*N):
ru:=evalf([seq(mu[i1+1]^2/mu[i1] -mu[i1+2],i1=1..nops(mu)-2)]):
print(ru):
if max(op(ru))>1-r or min(op(ru))<-r then
 RETURN(FAIL):
fi:

# This part is a little subtle. Yes, it tends to zero, and yes, we could find
# that N0. However, what guarantees that it will hold to N0? We only know
# _after_ N0, and up to what we've checked "by hand." I assume that there is
# some property or step that I missed that guarantees this, but I will take it
# for granted for now.
print(`The largest is`):
print(max(op(ru))):
print(`The smallest is`):
print(min(op(ru))):
print(`and as n goes to infinity, c(n) tends to zero (and it is possible (but a waste of time) to find an N0 such that |c(n)| is guaranteed to be`):
print(`less than `, min(r,1-r), ` for n>=N0 (or any epsilon for that matter), and check that for  for n<=N0. `):
print(`c(n)>= `, -r, ` and  c(n) < `, 1-r):
print(`QED.`):
RETURN(gu):
fi:
fi:
print(`-------------------------------------------------------------------------`):
print(``):
end:




#Epicyle(x,y,r,L,K): inputs positive integers x,y, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]
#and positive integers K1, and K2. First guesses whether there seems to be a linear recurrence equation with constant coefficients
#(given by its C-finite representation) of order L. If not, it returns FAIL.
#It is does, it checks whether it goes all the way to K. If it does it returns the recurrence. If it does not
#it return a more complicated recurrence with an indication how far it goes. Try
#Epicycle(5,21,0,5,1000);
Epicycle:=proc(x,y,r,L,K) local C,gu1,gu2,i,k:
C:=GuessRec(PS(x,y,r,2*L+10)):
if C=FAIL then
 RETURN(FAIL):
fi:

gu1:=SeqFromRec(C,K):

gu2:=PS(x,y,r,K):

if gu1=gu2 then
 RETURN(C):
fi:

while gu1<>gu2 do

 for i from 1 to K while gu1[i]=gu2[i] do od:

k:=gu1[i]-gu2[i]:

C:=[[op(1..i,gu2)],[op(C[2]),0$(i-nops(C[2])-1),-k]]:

if SeqFromRec(C,i+3)<>[op(1..i+3,gu2)] then
 RETURN(FAIL):
fi:

gu1:=SeqFromRec(C,K):

od:

C:
end:


#MamarT(K1,K2,r,Ord,HAMON): Inputs a positive integer K, a rational number r, and a positive integer Ord. Outputs
#all the pairs  2<=a<=K1, a+2<b<=K2 such that b/a is not an integer and such that the Pisot sequence with parameter r
#and initial conditions a,b, has a provable recurrence of order <=Ord. It also returns the
#successes and failues.  
#It returns, the pairs that produce provable non-trivial recurrences, those giving linear sequences, 
#and the set of false recurrences
#and the set of failures. It all list, at the end, all those false alarams whose longevity is>=HAMON
#Try
#MamarT(10,10,1/2,6);

MamarT:=proc(K1,K2,r,Ord,HAMON) local gu,a,b,kvS,kvF,PISOT,kvL,kvFA,fu,lu,aluf,si,zman0,METU,i1	:

zman0:=time():

print(`A list of all Pisot sequences with initial terms a<`,K1,`a+2<=b<=`,K2, `that satisfy `):
print(` Non-Trivial Linear Recurrences With Constant Coefficients of Order at most`, Ord):
print(``):
print(`By Shalosh B. Ekhad `):
print(``):
kvS:={}:
kvF:={}:
kvL:={}:
kvFA:={}:
METU:={}:
for a from 2 to K1 do
 for b from a+2 to K2 do
   if not type(b/a,integer) then
    gu:=PtoR(a,b,r,Ord,10*Ord+100):
    if gu<>FAIL and gu[1]<>FAIL and nops(gu)=4 and gu[4]=LINEAR  then
     kvL:=kvL union {[a,b]}:

elif gu<>FAIL and gu[1]=FAIL and nops(gu)=4 and gu[4]>1 then
 kvFA:=kvFA union {[a,b]}:

    elif gu<>FAIL and gu[1]<>FAIL and nops(gu)=4 and gu[4]<>FAIL then
     kvS:=kvS union {[a,b]}:
        lprint(PISOT(a,b)=gu[1]):
       else
       kvF:=kvF union {[a,b]}:
  fi:
 fi:
od:
od:

print(`out of the`, nops(kvF)+nops(kvS)+nops(kvL)+nops(kvFA), `cases, the following`, nops(kvS), `had (rigorously!) provable recurrences that are listed above`):
print(`Here they are`):
print(``):
lprint(kvS):
print(``):
print(`while the following `, nops(kvL), `give  linear trivial (linear) sequences. `):
print(``):
lprint(kvL):

print(``):
print(`while the following `, nops(kvFA), ` satisfy recurrences of order <=`,Ord):
print(` that eventually break down. Here they are where `):
print(`each line lists the recurrences, the first time it breaks down, and the Pisot index`):
print(`that is always larger than 1, explaining the break-down. `):
print(``):

aluf:=lu[1]:
si:=1:

for lu in kvFA do
 fu:=PtoR(op(lu),r,Ord,10*Ord+100):

    if fu[3]-1>HAMON then
        METU:=METU union {[lu,fu[3]-1,fu[4]]}:
     fi:
   if fu[3]-1>si then
      si:=fu[3]-1:
     aluf:=lu:
   fi:
print(PISOT(op(lu)) <> fu[2], fu[3]-1, fu[4]):
od:

print(`Finally,  the following `, nops(kvF), `do definitely not have a recurrence of order <= `,Ord):
print(``):
lprint(kvF):

print(``):
print(`To sum-up out of the`, nops(kvS)+nops(kvL)+nops(kvFA)+nops(kvF), `pairs of integers a<=`,K1, `a+2<=b<=`, K2):
print(``):
print(nops(kvS), `satisy provable recurrences of order <=`,Ord):
print(``):
print(nops(kvL), `are trivial, linear, sequences`):
print(``):

print(nops(kvF), `do not satisy linear recurrences of order <=`, Ord):
print(``):

print(nops(kvFA), `gets your hope up but then FAIL`):
print(``):

if METU<>{} then
print(`Out of them the following pairs lead to false recurrences with longevity>=`, HAMON, `together with their longevity and Pisot index`):
print(METU):

 aluf:=METU[1]:

 si:=aluf[2]:


  for i1 from 2 to nops(METU) do
   if METU[i1][2]>si then
     aluf:=METU[i1]:
     si:=METU[i1][2]:
   fi:
  od:

print(`The Pisot sequence with largest longevity is`, aluf, `that has longevity`, si):	
fi:
print(``):
print(`This concludes this article that took `, time()-zman0, `to generate`):
print(`I hope, dear readers, that you enjoyed reading it as much as I did generating it. `):
print(``):
print(`tam v'nishlam shevach le-el bore olam . `):
print(``):

[kvS,kvL,kvFA,kvF]:

end:


#FirstDev(x,y,r,C): inputs  positive integers x,y, and a parameter r (0<=r<=1) and a C-finite sequence C
#outputs the first n that violated the recurrence. Try:
#FirstDev(10,219,1/2,[[10, 219, 4796, 105030], [22, -3, 18, -11]]);
FirstDev:=proc(x,y,r,C) local gu,i,i1,R:

if Pis(C)<1 then
 RETURN(FAIL):
fi:


gu:=PS(x,y,r,nops(C[2])+2):
R:=nops(gu):
for i from R to 10^9 do

 gu:=[op(2..nops(gu),gu), trunc(gu[nops(gu)]^2/gu[nops(gu)-1]+r)]:
 if gu[nops(gu)]-add(C[2][i1]*gu[nops(gu)-i1],i1=1..nops(C[2]))<>0 then
   RETURN(i+1):
 fi:
od:
print(`The first deviation is after`, 10^9):
RETURN(FAIL):

end:





#PtoRvNumbered(x,y,r,L,N,K, MISPAR):  Like PtoRv(x,y,r,L,N,K) but with the theorem (or fact) numbered by MISPAR
#Try
#PtoRvNumbered(2,3,1/2,10,1000,10);
#PtoRvNumbered(5,17,1/2,10,1000,11);
PtoRvNumbered:=proc(x,y,r,L,N,K,MISPAR) local gu,mu,a,n,C,lu1,lu2,i,beta,ku,i1,t,aluf,lu,A,ru:

gu:=PtoR(x,y,r,L,K):


if gu[-1]=LINEAR then
print(`Fact`, MISPAR, `:The Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`whose `, N, `first terms `):
print(``):
print(PS(x,y,r,N)):
print(``):
print(`is a trivial linear sequence`):
RETURN(gu):
fi:

if gu=FAIL then
print(`Fact`, MISPAR, `:The Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`whose `, N, `first terms `):
print(``):
print(PS(x,y,r,N)):
print(``):
print(`does not satisfy`):
print(`a linear recurrence equation with constant coefficients of order<=`, L):
RETURN(gu):
fi:

mu:=PS(x,y,r,N):

if gu<>FAIL and gu[1]=FAIL and nops(gu)=4 then
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
print(`Fact`, MISPAR, `: Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`At first sight it seems to satisfy the following linear recurrence: `):
C:=gu[2]:
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):

print(`Alas, it breaks down at the`, gu[3], `-th term. `):
print(``):
lu2:=SeqFromRec(C,gu[3])[gu[3]]:
lu1:=PS(x,y,r,gu[3])[gu[3]]:
print(a(gu[3]), `equals `, lu1 , `while the corresponding term  for the solution of the recurrence is `, lu2):
print(``):
print(`So the difference of the former from the latter is`, lu1-lu2):
print(` `):
print(`Note that the Pisot Indicator is`, gu[4] ):

if gu[4]>1 and r>0 and r<1 then
print(`Since it is larger than 1, it is not at all suprising, that it does not go on for ever.`):
fi:
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
RETURN(gu):
fi:


if gu<>FAIL and gu[1]=FAIL and nops(gu)=3 then
print(``):
print(`-------------------------------------------------------------------------`):
print(``):
print(`Fact`, MISPAR,  `:Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`It  satisfies the following linear recurrence for n from 1 to `, K):
C:=gu[2]:
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):

print(`but it is not at all clear that this goes for ever, and it probably does not. `):

if r=1 then
print(`Note that the maximum of the Pisot index is`, gu[3][2] ):
fi:

if r=0 then
print(`Note that the minimum of the Pisot index is`, gu[3][1] ):
fi:

print(`You are welcome to try to make`, K , ` bigger. `):
RETURN(gu):
fi:





if r=0 and gu[2][1]<1/100 or r=1 and gu[2][2]>99/100 then
 RETURN(FAIL):
fi:

if gu<>FAIL and gu[1]<>FAIL and r<1 and r>0 then


beta:=gu[4]:

C:=gu[1]:



print(``):
print(`-------------------------------------------------------------------------`):
print(``):

if beta>1 then
print(`Fact`, MISPAR, `: Consider the Pisot Sequence a(n), defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(`[BTW, the  first`,N, `terms are:`, op(mu), ` ]. `):
print(``):
print(`At first sight it seems to satisfy the following linear recurrence: `):
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2])), ` . `):
print(``):
print(`However, this can't go on for ever, since the Pisot indicator`, beta, `is >= 1. `):

else

print(`Theorem `, MISPAR, `: Let a(n) be the Pisot Sequence with parameter`, r, `and with initial conditions`, a(1)=x, a(2)=y , ` i.e. for n>1`):
print(``):
print(a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(``):
print(`[Here are the first`, N, `terms `):
print(PS(x,y,r,N)):
print(`]`):
print(`The sequence a(n) satisfies, for n>=`, nops(C[2])+1, ` the linear recurrence equation with constant coefficient `):
print(``):
print(a(n)=add(C[2][i]*a(n-i),i=1..nops(C[2]))):
print(``):
print(`with initial conditions`, seq(a(i)=mu[i],i=1..nops(C[1]))  , ` . ` ):

print(``):
print(`Proof: Let b(n) be the solution of the above recurrence, i.e. the sequence defined by `):
print(``):
print(b(n)=add(C[2][i]*b(n-i),i=1..nops(C[2]))):
print(``):
print(`with initial conditions`, seq(b(i)=mu[i],i=1..nops(C[1])) , ` . `):
print(``):
print(`We have to prove that a(n)=b(n), but it is easier to prove the equivalent statement that b(n)=a(n). `):
print(``):
print(`It is readily checked that b(n)=a(n) for n from 1 to`, nops(C[1]), ` . `):
print(``):
print(`It remains to prove that b(n) satisfies the Pisot-recurrence `):
print(``):
print(b(n+2)=trunc(r+b(n+1)^2/b(n)) , ` . `)  :
print(``):
print(`But this is equivalent to the INEQALITIES:`):
print(``):
print(-r <=b(n+1)^2/b(n) -b(n+2) ,  b(n+1)^2/b(n)-b(n+2) <1-r ):
print(``):
print(`In other words`):
print(``):
print(-r <=normal(b(n+1)^2/b(n) -b(n+2)) ,  normal(b(n+1)^2/b(n)-b(n+2)) <=1-r ):
print(``):
print(`The characteristic equation in t,  of the recurrence satisfied by b(n) is`):

ku:=t^nops(C[2])- add(C[2][i1]*t^(nops(C[2])-i1),i1=1..nops(C[2])):
print(``):
print(sort(ku)=0):

print(``):

print(`whose roots are`):
print(``):
lu:=[solve(ku)]:
print(lu):

if member(1, {op(lu)}) then
 print(`Since 1 is a root, let's remove it, and the remaining roots are`):
 lu:={op(lu)} minus {1}:
 lu:=convert(lu,list):
fi:

print(`In floating-point`):
print(evalf(lu)):
aluf:=1:

for i from 2 to nops(lu) do
 if abs(evalf(lu[i]))>abs(evalf(lu[aluf])) then
  aluf:=i:
 fi:
od:

print(`The largest root is`, evalf(lu[aluf])):
print(``):
mu:=evalf([op(1..aluf-1,lu),op(aluf+1..nops(lu),lu)]):
print(`and the remaining roots are `):
print(``):
print(mu):
print(``):
print(`whose absolute values are`):
print(``):
print([seq(abs(mu[i1]),i1=1..nops(mu))]):
print(``):

if beta<1 then
print(`so the largest absolute value is`, gu[4]):
print(`that is less than 1,`):
else
print(`so the largest absolute value is`, gu[4]):
print(`that equals 1.`):

fi:

print(``):
print(`It follows that the sequence `):
print(``):
print(c(n)=normal(b(n+1)^2/b(n) -b(n+2))):
print(`satisfies the inequality `):
print(c(n)<=A*gu[4]^n ):
print(``):
if beta<1 then
print(`for some fixed constant, A, (independent of n), that can be easily determined, if desired. In particular the sequence c(n) tends to 0`):
else
print(`for some fixed constant, A, (independent of n), that can be easily determined, if desired. In particular the sequence c(n) is bounded.`):
fi:

print(`we have to show that for all n, c(n) is between `, -r, `and `, 1-r ):
print(`For the first`,  3*N-2,  `terms the sequence c(n) is`):
print(``):
mu:=PS(x,y,r,3*N):
ru:=evalf([seq(mu[i1+1]^2/mu[i1] -mu[i1+2],i1=1..nops(mu)-2)]):
print(ru):
if max(op(ru))>1-r or min(op(ru))<-r then
 RETURN(FAIL):
fi:

print(`The largest is`):
print(max(op(ru))):
print(`The smallest is`):
print(min(op(ru))):
print(`and as n goes to infinity, c(n) tends to zero (and it is possible (but a waste of time) to find an N0 such that |c(n)| is guaranteed to be`):
print(`less than `, min(r,1-r), ` for n>=N0 (or any epsilon for that matter), and check that for  for n<=N0. `):
print(`c(n)>= `, -r, ` and  c(n) < `, 1-r):
print(`QED.`):
RETURN(gu):
fi:
fi:
print(`-------------------------------------------------------------------------`):
print(``):
end:


#MamarV(K1,K2,r,Ord): Inputs a positive integer K, a rational number r, and a positive integer Ord. Outputs
#an article about
#all the pairs  2<=a<=K1, a+2<b<=K2 such that b/a is not an integer and such that the Pisot sequence with parameter r
#and initial conditions a,b, has a provable recurrence of order <=Ord. 
#It gives complete proof (modulo routine verifications) ot all its correct statements.
#It returns, the pairs that produce provable non-trivial recurrences, those giving linear sequences, 
#and the set of false recurrences
#and the set of failures
#Try
#MamarV(10,10,1/2,6);

MamarV:=proc(K1,K2,r,Ord) local gu,a,b,kvS,kvF,PISOT,kvL,kvFA,fu,lu,aluf,si,MI1,MI2,zman0:

zman0:=time():

print(`Theorems and Proofs about  Pisot sequences with initial terms a<`,K1,`a+2<=b<=`,K2, `that satisfy `):
print(` Non-Trivial Linear Recurrences With Constant Coefficients of Order at most`, Ord):
print(``):
print(`By Shalosh B. Ekhad `):
print(``):
kvS:={}:
kvF:={}:
kvL:={}:
kvFA:={}:

MI1:=0:
MI2:=0:

for a from 2 to K1 do
 for b from a+2 to K2 do
   if not type(b/a,integer) then
    gu:=PtoR(a,b,r,Ord,10*Ord+100):
    if gu<>FAIL and gu[1]<>FAIL and nops(gu)=4 and gu[4]=LINEAR  then
     kvL:=kvL union {[a,b]}:
     MI2:=MI2+1:
   PtoRvNumbered(a,b,r,Ord,20,10*Ord+100,MI2):
   print(``):
elif gu<>FAIL and gu[1]=FAIL and nops(gu)=4 and gu[4]>1 then
 kvFA:=kvFA union {[a,b]}:
     MI2:=MI2+1:
   PtoRvNumbered(a,b,r,Ord,20,10*Ord+100,MI2):
 


    elif gu<>FAIL and gu[1]<>FAIL and nops(gu)=4 and gu[4]<>FAIL then
     kvS:=kvS union {[a,b]}:
      MI1:=MI1+1:
      PtoRvNumbered(a,b,r,Ord,20,10*Ord+100,MI1):
       else
       kvF:=kvF union {[a,b]}:
      MI2:=MI2+1:
   PtoRvNumbered(a,b,r,Ord,20,10*Ord+100,MI2):
 
  fi:
 fi:
od:
od:

print(`To sum up,`):
print(`out of the`, nops(kvF)+nops(kvS)+nops(kvL)+nops(kvFA), `cases, the following`, nops(kvS), `had (rigorously!) provable recurrences that are listed above`):
print(`Here they are`):
print(``):
lprint(kvS):
print(``):
if kvL<>{} then
print(`while the following `, nops(kvL), `give  linear trivial (linear) sequences. `):
print(``):
lprint(kvL):
fi:

if kvFA<>{} then
print(``):
print(`while the following `, nops(kvFA), ` satisfy recurrences of order <=`,Ord):
print(` that eventually break down. Here they are where `):
print(`each line lists the recurrences, the first time it breaks down, and the Pisot index`):
print(`that is always larger than 1, explaining the break-down. `):
print(``):

aluf:=lu[1]:
si:=1:


for lu in kvFA do
 fu:=PtoR(op(lu),r,Ord,10*Ord+100):

   if fu[3]-1>si then
      si:=fu[3]-1:
     aluf:=lu:
   fi:
print(PISOT(op(lu)) <> fu[2], fu[3]-1, fu[4]):
od:

print(``):
print(`Among the false recurrences the PISOT sequence `, PISOT(op(aluf),r), `lasted longest, up to n=`,si):
print(``):
print(`Finally,  the following `, nops(kvF), ` pairs do definitely not produce (as initial values), Pisot sequences with parameter`, r):
print(` that have a recurrence of order <= `,Ord):
print(``):
lprint(kvF):
fi:
print(``):
print(`To sum-up out of the`, nops(kvS)+nops(kvL)+nops(kvFA)+nops(kvF), `pairs of integers a<=`,K1, `a+2<=b<=`, K2):
print(``):
print(nops(kvS), `satisy provable recurrences of order <=`,Ord):
print(``):
print(nops(kvL), `are trivial, linear, sequences`):
print(``):
print(nops(kvFA), `gets your hope up but then FAIL`):
print(``):
print(nops(kvF), `do not satisy linear recurrences of order <=`, Ord):
print(``):

print(`This concludes this article that contained`, MI1, `theorems and`, MI2, `facts, and took `, time()-zman0, `to generate`):
print(`I hope, dear readers, that you enjoyed reading it as much as I did generating it. `):
print(``):
print(`tam v'nishlam shevach le-el bore olam . `):
print(``):

[kvS,kvL,kvFA,kvF]:

end:


#PtoRb(x,y,r,L,K): inputs positive integers x,y, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]
#and positive integers K1, and K2. First guesses whether there seems to be a linear recurrence equation with constant coefficients
#(given by its C-finite representation) of order L. If not, it returns FAIL.
#If it persists to K terms (make K large), then it returns
#the C-finite description, followed by its Pisot-index-pair for these K terms, and the ones for the second half.
#(so it returns a list of length 3 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex, PisotIndicator] ]
#If the initial conjecture fails, then it returns it returns [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:
#PtoRb(5,21,0,6,1000);
PtoRb:=proc(x,y,r,L,K) local C,gu1,gu2,i,ku1,ku2:

if r< 0 or r>1 then
 print(`r must be weakly between 0 and 1 for this procedure`):
 RETURN(FAIL):
fi:

C:=GuessRec(PS(x,y,r,2*L+10)):
if C=FAIL then
 RETURN(FAIL):
fi:

gu1:=SeqFromRec(C,K):

gu2:=PS(x,y,r,K):

ku1:=PisInd(C,r,1,K):
if gu1<>gu2 then
 for i from 1 to K while gu1[i]=gu2[i] do od:
 RETURN([FAIL,C,i,ku1,Pis(C)]):
fi:


if (r=0 and ku1[1]<0) or (r=1 and ku1[2]>=1) then
 RETURN( [FAIL, C, i,ku1,Pis(C)]):
fi:

ku2:=evalf(PisInd(C,r,trunc(K/2),K),10):

if C[2]=[2,-1] then

 [C, evalf(ku1), ku2, LINEAR ]:

elif Pis(C)>1 then
 for i from 1 to K while gu1[i]=gu2[i] do od:
 RETURN([FAIL,C, i , ku1,ku2,Pis(C) ]):

else


 [C, evalf(ku1), ku2, Pis(C) ]:


fi:

end:


#Boyd(a0,k,c,r,L,K): a recurrence of order <=L valid up to K terms of the sequence PS(a0,k*a0^2+c,r). Try:
#Boyd(3,1,1,0,6,1000);
Boyd:=proc(a0,k,c,r,L,K)  local C,gu1,gu2,i,x,y:

x:=a0:
y:=k*a0^2+c:

C:=GuessRec(PS(x,y,r,2*L+10)):
if C=FAIL then
 RETURN(FAIL):
fi:

gu1:=SeqFromRec(C,K):

gu2:=PS(x,y,r,K):

if gu1<>gu2 then
 for i from 1 to K while gu1[i]=gu2[i] do od:
 RETURN([FAIL,C,i,Pis(C)]):
else
RETURN(C):
fi:


end:


#BoydP(n,a): The Boyd polynomial in a satisfying P(n) is the quotient if P(n-1)^2  by P(n-2) or this quotient
#minus P(n-2) if this is necessary to make the leading coefficient of the remainder negative. Try:
#BoydP(5,a);
BoydP:=proc(n,a) local R,Q:
option remember:

if n=0 then
 RETURN(a):
fi:

if n=1 then
 RETURN(a^2-1):
fi:

Q:=quo(BoydP(n-1,a)^2,BoydP(n-2,a),a):
R:=rem(BoydP(n-1,a)^2,BoydP(n-2,a),a):

if lcoeff(R,a)<0 then
RETURN(Q):
else
 RETURN(Q-R):
fi:

end:


#GPr(x,n): The template for the generalized Pisot recurrence of order n, expressing x[0] in terms of x[-1], ..., x[-n]
#such that the nby n determinant det([[x[-2*(n-1)], ..., x[-(n-1)]],[x[-2*(n-1)+1], ...,x[-n+2]], ..., [x[-(n-1)], ..., x[0]].
#Try:
#GPr(x,3);
GPr:=proc(x,n) local M,i,j,gu:
option remember:
M:= [seq([seq(x[-2*(n-1)+i+j],j=0..n-1)],i=0..n-1)]:

gu:=expand(linalg[det](M)):

normal(-coeff(gu,x[0],0)/coeff(gu,x[0],1)):
end:




#PSg(INI,r,K): the first K terms of the Generalized Pisot sequence of order n, defined by
#a[i]=INI[i], where INI is a list of even length 2*(n-1)
#PSg([2,5,9,10],1/2,20);
PSg:=proc(INI,r,K) local mu,gu,i,x,i1,n,mu1,mu2:

if not (nops(INI) mod 2=0 and nops(INI)>=2) then
 print(`The length of INI must be an even positive integer`):
fi:
n:=nops(INI)/2+1:


mu:= GPr(x,n):

mu2:=denom(mu):

gu:=INI:
for i from 1 to K-2*(n-1) do
if subs({seq(x[-i1]=gu[nops(gu)-i1+1],i1=1..2*(n-1))},mu2)=0 then
 RETURN(FAIL):
else
 gu:=[op(gu), trunc(subs({seq(x[-i1]=gu[nops(gu)-i1+1],i1=1..2*(n-1))},mu)+r)]:
fi:
od:



end:

#ShoAbs(C): the sorted (from largest to smallest of the absolute values of the characteristic equation of the recurrence C
#Try:
#ShoAbs([[1,1],[1,1]]);
ShoAbs:=proc(C) local mu,i:
mu:=evalf(Sho(C)):
mu:=[seq(abs(mu[i]),i=1..nops(mu))]:
mu:=sort(mu):
[seq(mu[nops(mu)-i+1],i=1..nops(mu))]:
end:


#DetSeq(L,r): Inputs a sequence L and outputs a sequence of length nops(L)-2*(r-1) whose
#i-th entry is det([[L[i], ..., L[i+r-1]], [L[i+1], ..., L[i+r+1]],  ..., [L[i+r-1], ..., L[i+2*(r-1)]].
#Try:
#DetSeq([seq(i^4,i=1..10)],3);
DetSeq:=proc(L,r) local i1,j1,n1:

[
seq(
det( [ seq([seq(L[n1+i1+j1],j1=0..r-1)],i1=0..r-1)]),
n1=1..nops(L)-2*(r-1))]:
end:




#PisIndG(C,k,r,K1,K2): inputs a C-finite description C, describing a sequence, let's call it L, and  positive integers K1, K2, outputs
#the min and max from n=K1 to K2 of the quantity DetSeq(L,k)/DetSeq(L,k-1).
#They should be >=0 and <1 respectively. Try:
#PisIndG([[2,3],[1,1]],2,1/2,100,200);
#PisIndG([[4,7,12],[2,-1,1]],2,1/2,100,200);
PisIndG:=proc(C,k,r,K1,K2) local L,lu,L1,L2,i1:
L:=SeqFromRec(C,K2+2*k):
L1:=DetSeq(L,k):
L2:=DetSeq(L,k-1):
lu:=[]:

for i1 from K1 to K2 do
 if L2[i1]=0 then
 RETURN(FAIL):
 else
lu:=[op(lu),  -L1[i1]/L2[i1]+r]:
fi:
od:
[min(op(lu)),max(op(lu))]:
end:

#PtoRg(INI,r,L,K): inputs a list of positive integers INI, and a rational number r between 0 and 1 (inclusive) [typically 0,1/2,1]
#and positive integer K. First guesses whether there seems to be a linear recurrence equation with constant coefficients
#(given by its C-finite representation) of order L for PSg(INI,r). If not, it returns FAIL.
#If it persists to K terms (make K large), then it returns
#the C-finite description, followed by its generalized
#Pisot-index-pair for these K terms, and the ones for the second half.
#(so it returns a list of length 3 [C, [TotalMinPisIndex,TotalMaxPisIndex], [SecondHalfMinPisIndex,SecondHalfMaxPisIndex, PisotIndicator] ]
#If the initial conjecture fails, then it returns it returns [FAIL,C,FirstPlaceOfDisagreement, PisotIndicator]. Try:
#PtoRg([2,3],1/2,10,1000);
#PtoRg([5,17],1/2,10,1000);
#PtoRg([1,4,,8,23],1/2,10,1000);

PtoRg:=proc(INI,r,L,K) local C,gu1,gu2,i,ku1,ku2,dor:

if r<=0 or r>=1 then
 print(`r must be strictly between 0 and 1 for this procedure`):
 RETURN(FAIL):
fi:

C:=GuessRec( PSg(INI,r,2*L+10+nops(INI)) ) :
if C=FAIL then
 RETURN(FAIL):
fi:

gu1:=SeqFromRec(C,K):

gu2:=PSg(INI,r,K):


if gu1<>gu2 then
 for i from 1 to K while gu1[i]=gu2[i] do od:
 RETURN([FAIL,C,i,ShoAbs(C)[nops(INI)/2+1] ] ):
fi:



ku1:=evalf(PisIndG(C,nops(INI)/2+1,r,1,K),10):



if (r=0 and ku1[1]<1/10000) or (r=1 and ku1[2]>1-1/10000) then
 RETURN( [FAIL, C,ku1]):
fi:



ku2:=evalf(PisIndG(C,nops(INI)/2+1,r,trunc(K/2),K),10):

if nops(ShoAbs(C))<=2 then
 dor:=0:
else
dor:= ShoAbs(C)[nops(INI)/2+1]:
fi:

if C[2]=[2,-1] then

 [C, ku1, ku2, LINEAR ]:


elif  dor >1 then

 [FAIL,C, FirstDev(x,y,r,C), dor ]:

else

 [C, ku1, ku2, dor ]:
fi:

end:





#GuessLinear(L,k): Guesses a linear expression a+b*k such that L[i]=a+b*i for i=1..k . 
#L should have at least four terms.  Try:
#GuessLinear([seq(2*i+5,i=1..10)],k); 
GuessLinear:=proc(L,k) local mu,a,b,i:
if (not type(L,list) and nops(L)>=4) then
print(L, `should have been a list of at least four elements`):
RETURN(FAIL):
fi:

mu:={seq(L[i]-L[i-1],i=2..nops(L))}:

if nops(mu)<>1 then
 RETURN(FAIL):
fi:
b:=mu[1]:
a:=L[1]-b:
a+b*k:
end:



#PtoRsy(x,y0,r,k,ORD,Kama): inputs a positive integer x, another one, y0, between 1 and x^2-1 (not divisible by x),
#a rational number r, between 0 and 1, a positive integer, ORD, at least 2, and a positive integer Kama, at least 4
#tries (by doing Kama cases), to conjecture a  recurrence of order<=ORDER whose coefficients are linear in k
#that  for the Pisot Sequence E_r(x,x^2*k+y0,r)
#valid for all k>=0, just giving "infinitely" many facts. It returns FAIL. Try:
#PtoTsy(2,1,1/2,k,2,10);
PtoRsy:=proc(x,y0,r,k,ORD,Kama) local k0, yeda,seder,i1,j1,mu,hal,kha,yofee,INI,bdok:

if not (type(x,integer) and x>=2 and y0>=1 and y0<=x^2-1 and y0 mod x<>0 and Kama>=4) then
 print(`Bad input`):
 RETURN(FAIL):
fi:


yeda:=[seq(PtoR(x,y0+k0*x^2,r,ORD,10*ORD+10) ,k0=1..Kama)  ]:

if member(FAIL,{op(yeda)}) then
 RETURN(FAIL):
fi:

seder:={seq(nops(yeda[i1][1][2]),i1=1..nops(yeda))}:

if nops(seder)<>1 then
 print(`No uniform order`):
 RETURN(FAIL):
fi:

seder:=seder[1]:

if seder>Kama-1 then
 print(Kama, `is too small, make it bigger`):
 RETURN(FAIL):
fi:

mu:=[]:

for  j1 from 1 to seder do
 hal:=GuessLinear([seq(yeda[i1][1][2][j1],i1=1..Kama)],k):
 if hal=FAIL then
  RETURN(FAIL):
 else
 mu:=[op(mu),hal]:
fi:

od:

INI:=[x,x^2*k+y0]:

for j1 from 3 to seder do

kha:=GuessPol1([seq(yeda[i1][1][1][j1],i1=1..nops(yeda))],j1-1,k):


if kha=FAIL then
 RETURN(FAIL):
else
INI:=[op(INI),kha]:
fi:
od:

yofee:=[INI,mu]:

if [seq(yeda[i1][1],i1=1..Kama)]<>[seq(subs(k=i1,yofee),i1=1..Kama)] then
# print(yofee, `did not work out`):
 RETURN(FAIL):
fi:

bdok:=PtoR(x,y0+1000*x^2,r,ORD,1000):

if subs(k=1000,yofee)<>bdok[1] then
 print(yofee, `did not work out for k=1000`):
fi:

if bdok[nops(bdok)]>0.97 then
 print(`Take it with a grain of salt, for k=1000 the Pisot index is`, bdok[nops(bdok)]):
fi:

yofee:

end:









#MamarSyT(X,k,r,Ord)  inputs a positive integer X, a symbol k, a positive integer r striclty between 0 and 1
#(typically 1/2), a positive integer Ord, outputs all the recurrences for the infinite families
#PS(x,y0+k*x^2,r) in terms of the symbol k, valid for k>=1 (and if y0>x, probably also for k=0)  for all y0
#from 1 to x^2-1 that are NOT multiples of x, for all 2<=x<=X. Try:
# MamarSyT(5,k,1/2,6);
MamarSyT:=proc(X,k,r,Ord) local zman0,kishalon,PISOT,x,y0,gu,a,b:


zman0:=time():
kishalon:={}:
print(`A list of Succesful Recurrences or Order<=`, Ord,  `for the Infinite Families of Pisot Sequence with parameter`, 1/2):
print(`with initial conditions x, x^2*k+y, for x from 2 to `, X, `and y from 1 to x^2-1 that is not a multiple of x`):
print(`For GENERAL (symbolic!) k`):
print(``):
print(`By Shalosh B. Ekhad `):
print(``):


print(`The Pisot Sequence a(n)=PS(x,y)(n), is  defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(``):
print(`In this  gripping article, we will list PROVED recurrences for the Pisot Sequenes`):
print(``):
print(PISOT(x,k*x^2+y)):
print(``):
print(`for all (numeric) integers x between 2 and`, X, `and all numeric integers, y, that are not divisible by x from 1 to x^2-1.` ):
print(`We will use the convention that a linear recurrence of order r `):
print(``):
print(`a[n]=C1*a[n-1]+ ... +Cr*a[n-r]`):
print(``):
print(`is abbreviated as a pair of lists: [[a[1], ..., a[r]],[C1, ..., Cr]], where the first list consists of  the r initial values . `):
print(``):

for  x from 2 to X do
 for y0 from 1 to x^2-1 do
   if not y0 mod x=0 then
    gu:=PtoRsy(x,y0,r,k,Ord,4*Ord) :

   if gu=FAIL then
     kishalon:=kishalon union {[x,y0]}:
   else
   print(``):
   print(PISOT(x,y0	+k*x^2)=gu):
  print(``):
  fi:
 fi:
od:
od:

if kishalon={} then
 print(`Everything with x<=`, X):
 print(`has a recurrence`):
else
print(`The following pairs [a,b]`, PISOT(a,b+k*a^2,r),  `do not  posses recurrences of order<=`, Ord):
print(``):
print(kishalon):
print(``):
print(`and it is very possible that none exist. `):
print(``):
fi:
print(`This ends this fascinating article that took`,time()-zman0, `to generated. `):
end:


#MamarSyTtex(X,k,r,Ord)  Like MamarSyT(X,k,r,Ord) but typset with LaTeX
#(typically 1/2), a positive integer Ord, outputs all the recurrences for the infinite families
#PS(x,y0+k*x^2,r) in terms of the symbol k, valid for k>=1 (and if y0>x, probably also for k=0)  for all y0
#from 1 to x^2-1 that are NOT multiples of x, for all 2<=x<=X. Try:
# MamarSyT(5,k,1/2,6);
MamarSyTtex:=proc(X,k,r,Ord) local zman0,kishalon,PISOT,x,y0,gu,a,b:


zman0:=time():
kishalon:={}:
print(`A list of Succesful Recurrences or Order<=`, Ord,  `for the Infinite Families of Pisot Sequence with parameter`, 1/2):
print(`with initial conditions x, x^2*k+y, for x from 2 to `, X, `and y from 1 to x^2-1 that is not a multiple of x`):
print(`For GENERAL (symbolic!) k`):
print(``):
print(`By Shalosh B. Ekhad `):
print(``):


print(`The Pisot Sequence a(n)=PS(x,y)(n), is  defined by`, a(1)=x,a(2)=y, `and for n>1`, a(n+2)=trunc(r+a(n+1)^2/a(n))):
print(``):
print(`In this  gripping article, we will list PROVED recurrences for the Pisot Sequenes`):
print(``):
print(PISOT(x,k*x^2+y)):
print(``):
print(`for all (numeric) integers x between 2 and`, X, `and all numeric integers, y, that are not divisible by x from 1 to x^2-1.` ):
print(`We will use the convention that a linear recurrence of order r `):
print(``):
print(`a[n]=C1*a[n-1]+ ... +Cr*a[n-r]`):
print(``):
print(`is abbreviated as a pair of lists: [[a[1], ..., a[r]],[C1, ..., Cr]], where the first list consists of  the r initial values . `):
print(``):

for  x from 2 to X do
 for y0 from 1 to x^2-1 do
   if not y0 mod x=0 then
    gu:=PtoRsy(x,y0,r,k,Ord,4*Ord) :

   if gu=FAIL then
     kishalon:=kishalon union {[x,y0]}:
   else
   print(``):
   print(latex(PISOT(x,y0	+k*x^2)=gu)):
  print(``):
  fi:
 fi:
od:
od:

if kishalon={} then
 print(`Everything with x<=`, X):
 print(`has a recurrence`):
else
print(`The following pairs [a,b]`, PISOT(a,b+k*a^2,r),  `do not  posses recurrences of order<=`, Ord):
print(``):
print(kishalon):
print(``):
print(`and it is very possible that none exist. `):
print(``):
fi:
print(`This ends this fascinating article that took`,time()-zman0, `to generated. `):
end:
