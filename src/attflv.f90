Subroutine attflv(id, idq, idqq)
  Common /rndf77/nseed
  Save
  If (abs(id)<100) Then
    nsign = 1
    idq = id/100
    idqq = -id/10 + idq*10
    If (abs(idq)==3) nsign = -1
    idq = nsign*idq
    idqq = nsign*idqq
    If (idq<0) Then
      id0 = idq
      idq = idqq
      idqq = id0
    End If
    Return
  End If
  idq = 2
  If (abs(id)==2112) idq = 1
  idqq = 2101
  x = ranart(nseed)
  If (x<=0.5) Goto 30
  If (x>0.666667) Goto 10
  idqq = 2103
  Goto 30
  10 idq = 1
  idqq = 2203
  If (abs(id)==2112) Then
    idq = 2
    idqq = 1103
  End If
  30 If (id<0) Then
    id00 = idqq
    idqq = -idq
    idq = -id00
  End If
  Return
End Subroutine attflv
