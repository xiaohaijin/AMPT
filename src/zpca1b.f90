Subroutine zpca1b(rapi, et, ian)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para6/centy
  Common /ilist6/t, iopern, icolln
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  Save
  If (rapi>centy-0.5D0 .And. rapi<centy+0.5D0) Then
     det2(ian) = det2(ian) + et
     dn2(ian) = dn2(ian) + 1D0
     If (ian==10) Then
     End If
     If (ian==11) Then
     End If
     If (ian==12) Then
     End If
     If (rapi>centy-0.25D0 .And. rapi<centy+0.25D0) Then
        det1(ian) = det1(ian) + et
        dn1(ian) = dn1(ian) + 1D0
        If (rapi>centy-0.1D0 .And. rapi<centy+0.1D0) Then
           det(ian) = det(ian) + et
           dn(ian) = dn(ian) + 1D0
        End If
     End If
  End If
  Return
End Subroutine zpca1b
