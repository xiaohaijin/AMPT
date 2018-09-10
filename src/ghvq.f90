Function ghvq(y1, y2, amt2)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Dimension f(2, 7)
  Save
  xt = 2.0D0*dsqrt(amt2)/dble(hint1(1))
  x1 = 0.5D0*xt*(dexp(y1)+dexp(y2))
  x2 = 0.5D0*xt*(dexp(-y1)+dexp(-y2))
  ss = x1*x2*dble(hint1(1))**2
  af = 4.0D0
  If (ihpr2(18)/=0) af = 5.0D0
  dlam = dble(hipr1(15))
  aph = 12.0D0*3.1415926D0/(33.D0-2.D0*af)/dlog(amt2/dlam**2)
  Call parton(f, x1, x2, amt2)
  gqq = 4.D0*(dcosh(y1-y2)+dble(hipr1(7))**2/amt2)/(1.D0+dcosh(y1-y2))/9.D0*(f(1,1)*f(2,2)+f(1,2)*f(2,1)+f(1,3)*f(2,4)+f(1,4)*f(2,3)+f(1,5)*f(2,6)+f(1,6)*f(2,5))
  ggg = (8.D0*dcosh(y1-y2)-1.D0)*(dcosh(y1-y2)+2.D0*dble(hipr1(7))**2/amt2-2.D0*dble(hipr1(7))**4/amt2**2)/(1.D0+dcosh(y1-y2))/24.D0*f(1, 7)*f(2, 7)
  ghvq = (gqq+ggg)*dble(hipr1(23))*3.14159D0*aph**2/ss**2
  Return
End Function ghvq
