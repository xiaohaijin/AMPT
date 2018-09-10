Subroutine chkin2(l, i1, i2, i3, t, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Save
  itest = 0
  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ib = j
           ic = k
           If (k>=1 .And. k<=10) Then
              If (i==0) ia = 10
              If (i==11) ia = 1
              If (j==0) ib = 10
              If (j==11) ib = 1
              Call chkcel(l, ia, ib, ic, t, tmin, nc)
           End If
        End Do
     End Do
  End Do
  Return
End Subroutine chkin2
