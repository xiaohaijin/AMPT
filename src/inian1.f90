Subroutine inian1
  Implicit Double Precision (A-H, O-Z)
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Common /ana1/ts(12)
  Save
  If (ibstfg/=0) Then
     a = cosh(6D0)
     Do i = 1, 12
        ts(i) = ts(i)*a
     End Do
  End If
  Return
End Subroutine inian1
