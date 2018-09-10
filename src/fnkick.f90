Function fnkick(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  fnkick = 1.0/(x+hipr1(19)**2)/(x+hipr1(20)**2)/(1+exp((sqrt(x)-hipr1(20))/0.4))
  Return
End Function fnkick
