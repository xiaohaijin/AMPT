Subroutine hboost
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save
  Do i = 1, n
    dbeta = dble(p(i,3)/p(i,4))
    If (abs(dbeta)>=1.D0) Then
      db = dble(hint1(2))
      If (db>0.99999999D0) Then
        Write (6, *) '(HIBOOT:) boost vector too large'
        db = 0.99999999D0
      End If
      dga = 1D0/sqrt(1D0-db**2)
      dp3 = dble(p(i,3))
      dp4 = dble(p(i,4))
      p(i, 3) = sngl((dp3+db*dp4)*dga)
      p(i, 4) = sngl((dp4+db*dp3)*dga)
      Goto 100
    End If
    y = 0.5*sngl(dlog((1.D0+dbeta)/(1.D0-dbeta)))
    amt = sqrt(p(i,1)**2+p(i,2)**2+p(i,5)**2)
    p(i, 3) = amt*sinh(y+hint1(3))
    p(i, 4) = amt*cosh(y+hint1(3))
  100 End Do
  Return
End Subroutine hboost
