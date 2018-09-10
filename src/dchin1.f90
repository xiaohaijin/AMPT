Subroutine dchin1(l, ii, i1, i2, i3, t)
  Implicit Double Precision (A-H, O-Z)
  Save
  itest = 0
  If (ii==1) Then
     If (i1==1) Goto 100
     If (i1==2) Then
        If (i2>=2 .And. i2<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     i = i1 - 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
  If (ii==2) Then
     If (i1==10) Goto 100
     If (i1==9) Then
        If (i2>=2 .And. i2<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     i = i1 + 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
  If (ii==3) Then
     If (i2==1) Goto 100
     If (i2==2) Then
        If (i1>=2 .And. i1<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     j = i2 - 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
  If (ii==4) Then
     If (i2==10) Goto 100
     If (i2==9) Then
        If (i1>=2 .And. i1<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     j = i2 + 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
  If (ii==5) Then
     If (i3==1) Goto 100
     If (i3==2) Then
        If (i1>=2 .And. i1<=9 .And. i2>=2 .And. i2<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     k = i3 - 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
  If (ii==6) Then
     If (i3==10) Goto 100
     If (i3==9) Then
        If (i1>=2 .And. i1<=9 .And. i2>=2 .And. i2<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     k = i3 + 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If
100 Continue
  Return
End Subroutine dchin1
