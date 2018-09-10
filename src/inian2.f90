Subroutine inian2
  Implicit Double Precision (A-H, O-Z)
  Common /para5/iconfg, iordsc
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  Save
  If (iconfg<=3) Then
     Do i = 1, 12
        det(i) = 0D0
        dn(i) = 0D0
        det1(i) = 0D0
        dn1(i) = 0D0
        det2(i) = 0D0
        dn2(i) = 0D0
     End Do
  End If
  Return
End Subroutine inian2
