Function g(y1, y2, pt2)
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
  Call parton(f, x1, x2, pt2)
  g11 = ((f(1,1)+f(1,2))*(f(2,3)+f(2,4)+f(2,5)+f(2,6))+(f(1,3)+f(1,4))*(f(2,5)+f(2,6)))*subcr1(t, u)
  g12 = ((f(2,1)+f(2,2))*(f(1,3)+f(1,4)+f(1,5)+f(1,6))+(f(2,3)+f(2,4))*(f(1,5)+f(1,6)))*subcr1(u, t)
  g13 = (f(1,1)*f(2,1)+f(1,2)*f(2,2)+f(1,3)*f(2,3)+f(1,4)*f(2,4)+f(1,5)*f(2,5)+f(1,6)*f(2,6))*(subcr1(u,t)+subcr1(t,u)-8.D0/t/u/27.D0)
  g2 = (af-1)*(f(1,1)*f(2,2)+f(2,1)*f(1,2)+f(1,3)*f(2,4)+f(2,3)*f(1,4)+f(1,5)*f(2,6)+f(2,5)*f(1,6))*subcr2(t, u)
  g31 = (f(1,1)*f(2,2)+f(1,3)*f(2,4)+f(1,5)*f(2,6))*subcr3(t, u)
  g32 = (f(2,1)*f(1,2)+f(2,3)*f(1,4)+f(2,5)*f(1,6))*subcr3(u, t)
  g4 = (f(1,1)*f(2,2)+f(2,1)*f(1,2)+f(1,3)*f(2,4)+f(2,3)*f(1,4)+f(1,5)*f(2,6)+f(2,5)*f(1,6))*subcr4(t, u)
  g5 = af*f(1, 7)*f(2, 7)*subcr5(t, u)
  g61 = f(1, 7)*(f(2,1)+f(2,2)+f(2,3)+f(2,4)+f(2,5)+f(2,6))*subcr6(t, u)
  g62 = f(2, 7)*(f(1,1)+f(1,2)+f(1,3)+f(1,4)+f(1,5)+f(1,6))*subcr6(u, t)
  g7 = f(1, 7)*f(2, 7)*subcr7(t, u)
  g = (g11+g12+g13+g2+g31+g32+g4+g5+g61+g62+g7)*dble(hipr1(17))*3.14159D0*aph**2/ss**2
  Return
End Function g
