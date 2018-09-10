Subroutine isco(i, j, allok, tm, t1, t2)
  Implicit Double Precision (A-H, O-Z)
  Common /para5/iconfg, iordsc
  Save
  Logical allok
  iorder = iordsc/10
  If (iconfg==1) Then
     If (iorder==1) Then
        Call isco1(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco2(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco3(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==2 .Or. iconfg==4) Then
     If (iorder==1) Then
        Call isco4(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco5(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco6(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==3) Then
     If (iorder==1) Then
        Call isco7(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco8(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco9(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==5) Then
     If (iorder==1) Then
        Call isco10(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco11(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco12(i, j, allok, tm, t1, t2)
     End If
  End If
  Return
End Subroutine isco
