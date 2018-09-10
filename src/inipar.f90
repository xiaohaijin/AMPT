Subroutine inipar
  Implicit Double Precision (A-H, O-Z)
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Common /para6/centy
  Save
  If (ibstfg/=0) Then
     centy = -6D0
  End If
  Return
End Subroutine inipar
