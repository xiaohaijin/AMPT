Subroutine zpcou
  Implicit Double Precision (A-H, O-Z)
  Common /para5/iconfg, iordsc
  Save
  If (iconfg<=3) Then
     Call zpcou1
  Else
     Call zpcou2
  End If
  Return
End Subroutine zpcou
