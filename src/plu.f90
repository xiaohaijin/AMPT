Function plu(i, j)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension psum(4)
  plu = 0.
  If (i<0 .Or. i>mstu(4) .Or. j<=0) Then
  Else If (i==0 .And. j<=4) Then
    Do i1 = 1, n
      If (k(i1,1)>0 .And. k(i1,1)<=10) plu = plu + p(i1, j)
    End Do
  Else If (i==0 .And. j==5) Then
    Do j1 = 1, 4
      psum(j1) = 0.
      Do i1 = 1, n
        If (k(i1,1)>0 .And. k(i1,1)<=10) psum(j1) = psum(j1) + p(i1, j1)
      End Do
    End Do
    plu = sqrt(max(0.,psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2))
  Else If (i==0 .And. j==6) Then
    Do i1 = 1, n
      If (k(i1,1)>0 .And. k(i1,1)<=10) plu = plu + luchge(k(i1,2))/3.
    End Do
  Else If (i==0) Then
  Else If (j<=5) Then
    plu = p(i, j)
  Else If (j<=12) Then
    If (j==6) plu = luchge(k(i,2))/3.
    If (j==7 .Or. j==8) plu = p(i, 1)**2 + p(i, 2)**2 + p(i, 3)**2
    If (j==9 .Or. j==10) plu = p(i, 1)**2 + p(i, 2)**2
    If (j==11 .Or. j==12) plu = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
    If (j==8 .Or. j==10 .Or. j==12) plu = sqrt(plu)
  Else If (j<=16) Then
    If (j<=14) plu = ulangl(p(i,3), sqrt(p(i,1)**2+p(i,2)**2))
    If (j>=15) plu = ulangl(p(i,1), p(i,2))
    If (j==14 .Or. j==16) plu = plu*180./paru(1)
  Else If (j<=19) Then
    pmr = 0.
    If (j==17) pmr = p(i, 5)
    If (j==18) pmr = ulmass(211)
    pr = max(1E-20, pmr**2+p(i,1)**2+p(i,2)**2)
    plu = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))
  Else If (j<=25) Then
    If (j==20) plu = 2.*sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)/paru(21)
    If (j==21) plu = 2.*p(i, 3)/paru(21)
    If (j==22) plu = 2.*sqrt(p(i,1)**2+p(i,2)**2)/paru(21)
    If (j==23) plu = 2.*p(i, 4)/paru(21)
    If (j==24) plu = (p(i,4)+p(i,3))/paru(21)
    If (j==25) plu = (p(i,4)-p(i,3))/paru(21)
  End If
  Return
End Function plu
