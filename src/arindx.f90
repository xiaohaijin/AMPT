Subroutine arindx(n, m, arrin, indx)
  Dimension arrin(n), indx(n)
  Save
  Do j = 1, m
    indx(j) = j
  End Do
  l = m/2 + 1
  ir = m
  10 Continue
  If (l>1) Then
    l = l - 1
    indxt = indx(l)
    q = arrin(indxt)
  Else
    indxt = indx(ir)
    q = arrin(indxt)
    indx(ir) = indx(1)
    ir = ir - 1
    If (ir==1) Then
      indx(1) = indxt
      Return
    End If
  End If
  i = l
  j = l + l
  20 If (j<=ir) Then
    If (j<ir) Then
      If (arrin(indx(j))<arrin(indx(j+1))) j = j + 1
    End If
    If (q<arrin(indx(j))) Then
      indx(i) = indx(j)
      i = j
      j = j + j
    Else
      j = ir + 1
    End If
    Goto 20
  End If
  indx(i) = indxt
  Goto 10
End Subroutine arindx
