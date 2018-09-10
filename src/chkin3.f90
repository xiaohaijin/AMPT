Subroutine chkin3(l, i1, i2, i3, t, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Save
  itest = 0
  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i==0) Then
              ia = 10
           Else If (i==11) Then
              ia = 1
           Else
              ia = i
           End If
           If (j==0) Then
              ib = 10
           Else If (j==11) Then
              ib = 1
           Else
              ib = j
           End If
           If (k==0) Then
              ic = 10
           Else If (k==11) Then
              ic = 1
           Else
              ic = k
           End If
           Call chkcel(l, ia, ib, ic, t, tmin, nc)
        End Do
     End Do
  End Do
  Return
End Subroutine chkin3
