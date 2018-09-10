Function gphotn(y1, y2, pt2)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Dimension f(2, 7)
  Save
  xt = 2.D0*dsqrt(pt2)/dble(hint1(1))
  x1 = 0.5D0*xt*(dexp(y1)+dexp(y2))
  x2 = 0.5D0*xt*(dexp(-y1)+dexp(-y2))
  z = dsqrt(1.D0-xt**2/x1/x2)
  ss = x1*x2*dble(hint1(1))**2
  t = -(1.D0-z)/2.D0
  u = -(1.D0+z)/2.D0
  af = 3.D0
  dlam = dble(hipr1(15))
  aph = 12.D0*3.1415926D0/(33.D0-2.D0*af)/dlog(pt2/dlam**2)
  aphem = 1.D0/137.D0
  Call parton(f, x1, x2, pt2)
  g11 = -(u**2+1.D0)/u/3.D0*f(1, 7)*(4.D0*f(2,1)+4.D0*f(2,2)+f(2,3)+f(2,4)+f(2,5)+f(2,6))/9.D0
  g12 = -(t**2+1.D0)/t/3.D0*f(2, 7)*(4.D0*f(1,1)+4.D0*f(1,2)+f(1,3)+f(1,4)+f(1,5)+f(1,6))/9.D0
  g2 = 8.D0*(u**2+t**2)/u/t/9.D0*(4.D0*f(1,1)*f(2,2)+4.D0*f(1,2)*f(2,1)+f(1,3)*f(2,4)+f(1,4)*f(2,3)+f(1,5)*f(2,6)+f(1,6)*f(2,5))/9.D0
  gphotn = (g11+g12+g2)*dble(hipr1(23))*3.14159D0*aph*aphem/ss**2
  Return
End Function gphotn
