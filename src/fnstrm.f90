Function fnstrm(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  fnstrm = 1.0/((1.0-x)**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)/(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
  Return
End Function fnstrm
