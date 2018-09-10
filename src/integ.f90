Integer Function integ(x)
  Implicit Double Precision (A-H, O-Z)
  Save
  If (x<0D0) Then
     integ = int(x-1D0)
  Else
     integ = int(x)
  End If
  Return
End Function integ
