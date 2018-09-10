Function fjetrg(x, wgt)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100), ptmax, ptmin
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Dimension x(10)
  Save
  ptmin = abs(hipr1(10)) - 0.25
  ptmin = max(ptmin, hipr1(8))
  am2 = 0.D0
  If (ihpr2(3)==3) Then
    am2 = dble(hipr1(7)**2)
    ptmin = max(0.0, hipr1(10))
  End If
  ptmax = abs(hipr1(10)) + 0.25
  If (hipr1(10)<=0.0) ptmax = hint1(1)/2.0 - sngl(am2)
  If (ptmax<=ptmin) ptmax = ptmin + 0.25
  pt2 = dble(ptmax**2-ptmin**2)*x(1) + dble(ptmin)**2
  amt2 = pt2 + am2
  xt = 2.0D0*dsqrt(amt2)/dble(hint1(1))
  ymx1 = dlog(1.0D0/xt+dsqrt(1.0D0/xt**2-1.0D0))
  y1 = 2.0D0*ymx1*x(2) - ymx1
  ymx2 = dlog(2.0D0/xt-dexp(y1))
  ymn2 = dlog(2.0D0/xt-dexp(-y1))
  y2 = (ymx2+ymn2)*x(3) - ymn2
  If (ihpr2(3)==3) Then
    gtrig = 2.0D0*ghvq(y1, y2, amt2)
  Else If (ihpr2(3)==2) Then
    gtrig = 2.0D0*gphotn(y1, y2, pt2)
  Else
    gtrig = g(y1, y2, pt2)
  End If
  fjetrg = 2.0D0*ymx1*(ymx2+ymn2)*dble(ptmax**2-ptmin**2)*gtrig/2.0D0
  Return
End Function fjetrg
