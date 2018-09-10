Function fnkc2(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  fnkc2 = x*exp(-2.0*x/hipr1(42))
  Return
End Function fnkc2
