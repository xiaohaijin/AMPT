Function pyw1au(eps, ireim)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  acosh(x) = log(x+sqrt(x**2-1.))
  If (eps<0.) Then
    w1re = 2.*sqrt(1.-eps)*asinh(sqrt(-1./eps))
    w1im = 0.
  Else If (eps<1.) Then
    w1re = 2.*sqrt(1.-eps)*acosh(sqrt(1./eps))
    w1im = -paru(1)*sqrt(1.-eps)
  Else
    w1re = 2.*sqrt(eps-1.)*asin(sqrt(1./eps))
    w1im = 0.
  End If
  If (ireim==1) pyw1au = w1re
  If (ireim==2) pyw1au = w1im
  Return
End Function pyw1au
