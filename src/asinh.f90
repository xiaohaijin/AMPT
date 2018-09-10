Function asinh(x)
  Save
  If (x>0) Then
    asinh = alog(x+sqrt(x**2+1.))
  Else
    asinh = -alog(-x+sqrt(x**2+1.))
  End If
  Return
End Function asinh
