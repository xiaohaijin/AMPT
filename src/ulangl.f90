Function ulangl(x, y)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  ulangl = 0.
  r = sqrt(x**2+y**2)
  If (r<1E-20) Return
  If (abs(x)/r<0.8) Then
    ulangl = sign(acos(x/r), y)
  Else
    ulangl = asin(y/r)
    If (x<0. .And. ulangl>=0.) Then
      ulangl = paru(1) - ulangl
    Else If (x<0.) Then
      ulangl = -paru(1) - ulangl
    End If
  End If
  Return
End Function ulangl
