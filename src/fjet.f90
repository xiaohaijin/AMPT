Function fjet(x, wgt)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Dimension x(10)
  Save
  pt2 = dble(hint1(1)**2/4.0-hipr1(8)**2)*x(1) + dble(hipr1(8))**2
  xt = 2.0D0*dsqrt(pt2)/dble(hint1(1))
  ymx1 = dlog(1.0D0/xt+dsqrt(1.0D0/xt**2-1.0D0))
  y1 = 2.0D0*ymx1*x(2) - ymx1
  ymx2 = dlog(2.0D0/xt-dexp(y1))
  ymn2 = dlog(2.0D0/xt-dexp(-y1))
  y2 = (ymx2+ymn2)*x(3) - ymn2
  fjet = 2.0D0*ymx1*(ymx2+ymn2)*dble(hint1(1)**2/4.0-hipr1(8)**2)*g(y1, y2, pt2)/2.0D0
  Return
End Function fjet
