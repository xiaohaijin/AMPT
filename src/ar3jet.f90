Subroutine ar3jet(s, x1, x3, jl)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /rndf77/nseed
  Save
  c = 1./3.
  If (k(jl,2)/=21 .And. k(jl+1,2)/=21) c = 8./27.
  exp1 = 3
  exp3 = 3
  If (k(jl,2)/=21) exp1 = 2
  If (k(jl+1,2)/=21) exp3 = 2
  a = 0.24**2/s
  yma = alog(.5/sqrt(a)+sqrt(.25/a-1))
  d = 4.*c*yma
  sm1 = p(jl, 5)**2/s
  sm3 = p(jl+1, 5)**2/s
  xt2m = (1.-2.*sqrt(sm1)+sm1-sm3)*(1.-2.*sqrt(sm3)-sm1+sm3)
  xt2m = min(.25, xt2m)
  ntry = 0
  1 If (ntry==5000) Then
    x1 = .5*(2.*sqrt(sm1)+1.+sm1-sm3)
    x3 = .5*(2.*sqrt(sm3)+1.-sm1+sm3)
    Return
  End If
  ntry = ntry + 1
  xt2 = a*(xt2m/a)**(ranart(nseed)**(1./d))
  ymax = alog(.5/sqrt(xt2)+sqrt(.25/xt2-1.))
  y = (2.*ranart(nseed)-1.)*ymax
  x1 = 1. - sqrt(xt2)*exp(y)
  x3 = 1. - sqrt(xt2)*exp(-y)
  x2 = 2. - x1 - x3
  neg = 0
  If (k(jl,2)/=21 .Or. k(jl+1,2)/=21) Then
    If ((1.-x1)*(1.-x2)*(1.-x3)-x2*sm1*(1.-x1)-x2*sm3*(1.-x3)<=0. .Or. x1<=2.*sqrt(sm1)-sm1+sm3 .Or. x3<=2.*sqrt(sm3)-sm3+sm1) neg = 1
    x1 = x1 + sm1 - sm3
    x3 = x3 - sm1 + sm3
  End If
  If (neg==1) Goto 1
  fg = 2.*ymax*c*(x1**exp1+x3**exp3)/d
  xt2m = xt2
  If (fg<ranart(nseed)) Goto 1
  Return
End Subroutine ar3jet
