Function omg0(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /besel/x4
  External bk
  Save
  x4 = hipr1(32)*sqrt(x)
  omg0 = hipr1(32)**2*gauss2(bk, x4, x4+20.0, 0.01)/96.0
  Return
End Function omg0
