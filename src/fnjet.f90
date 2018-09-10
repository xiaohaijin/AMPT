Function fnjet(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /njet/n, ipcrs
  Save
  omg1 = omg0(x)*hint1(11)/hipr1(31)
  c0 = exp(n*alog(omg1)-sgmin(n+1))
  If (n==0) c0 = 1.0 - exp(-2.0*omg0(x)*hipr1(30)/hipr1(31)/2.0)
  fnjet = c0*exp(-omg1)
  Return
End Function fnjet
