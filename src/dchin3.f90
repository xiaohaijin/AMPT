Subroutine dchin3(l, ii, i1, i2, i3, t)
  Implicit Double Precision (A-H, O-Z)
  Save
  If (ii==1) Then
     i = i1 - 2
     If (i<=0) i = i + 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  If (ii==2) Then
     i = i1 + 2
     If (i>=11) i = i - 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  If (ii==3) Then
     j = i2 - 2
     If (j<=0) j = j + 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  If (ii==4) Then
     j = i2 + 2
     If (j>=11) j = j - 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  If (ii==5) Then
     k = i3 - 2
     If (k<=0) k = k + 10
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  If (ii==6) Then
     k = i3 + 2
     If (k>=11) k = k - 10
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  Return
End Subroutine dchin3
