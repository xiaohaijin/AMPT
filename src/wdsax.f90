Function wdsax(x)
  Common /wood/r, d, fnorm, w
  Save
  wdsax = fnorm*(1.+w*(x/r)**2)/(1+exp((x-r)/d))
  If (w<0.) Then
    If (x>=r/sqrt(abs(w))) wdsax = 0.
  End If
  Return
End Function wdsax
