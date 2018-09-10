Function pyw2au(eps, ireim)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  acosh(x) = log(x+sqrt(x**2-1.))
  If (eps<0.) Then
    w2re = 4.*(asinh(sqrt(-1./eps)))**2
    w2im = 0.
  Else If (eps<1.) Then
    w2re = 4.*(acosh(sqrt(1./eps)))**2 - paru(1)**2
    w2im = -4.*paru(1)*acosh(sqrt(1./eps))
  Else
    w2re = -4.*(asin(sqrt(1./eps)))**2
    w2im = 0.
  End If
  If (ireim==1) pyw2au = w2re
  If (ireim==2) pyw2au = w2im
  Return
End Function pyw2au
