Function ftotjt(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  omg = omg0(x)*hint1(11)/hipr1(31)/2.0
  ftotjt = 1.0 - exp(-2.0*omg)
  Return
End Function ftotjt
