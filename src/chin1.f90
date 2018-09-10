Subroutine chin1(l, i1, i2, i3, last0, t, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Save
  itest = 0
  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call chcell(l, i, j, k, last0, t, tmin, nc)
           Else If (itest==0) Then
              m1 = 11
              m2 = 11
              m3 = 11
              Call chcell(l, m1, m2, m3, last0, t, tmin, nc)
              itest = 1
           End If
        End Do
     End Do
  End Do
  Return
End Subroutine chin1
