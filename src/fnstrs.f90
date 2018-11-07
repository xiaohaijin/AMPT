Function fnstrs(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  fnstrs = (1.0-x)**hipr1(47)/(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(48)
  Return
End Function fnstrs
