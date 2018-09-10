Function subcr3(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr3 = 4.D0/9.D0*(t**2+u**2+(1.D0+u**2)/t**2-2.D0*u**2/3.D0/t)
  Return
End Function subcr3
