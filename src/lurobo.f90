Subroutine lurobo(the, phi, bex, bey, bez)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension rot(3, 3), pr(3), vr(3), dp(4), dv(4)
  imin = 1
  If (mstu(1)>0) imin = mstu(1)
  imax = n
  If (mstu(2)>0) imax = mstu(2)
  dbx = dble(bex)
  dby = dble(bey)
  dbz = dble(bez)
  Goto 100
  Entry ludbrb(imi, ima, the, phi, dbex, dbey, dbez)
  imin = imi
  If (imin<=0) imin = 1
  imax = ima
  If (imax<=0) imax = n
  dbx = dbex
  dby = dbey
  dbz = dbez
  100 If (imin>mstu(4) .Or. imax>mstu(4)) Then
    Call luerrm(11, '(LUROBO:) range outside LUJETS memory')
    Return
  End If
  If ((the**2+phi**2)>1E-20) Then
    rot(1, 1) = cos(the)*cos(phi)
    rot(1, 2) = -sin(phi)
    rot(1, 3) = sin(the)*cos(phi)
    rot(2, 1) = cos(the)*sin(phi)
    rot(2, 2) = cos(phi)
    rot(2, 3) = sin(the)*sin(phi)
    rot(3, 1) = -sin(the)
    rot(3, 2) = 0.
    rot(3, 3) = cos(the)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 130
      Do j = 1, 3
        pr(j) = p(i, j)
        vr(j) = v(i, j)
      End Do
      Do j = 1, 3
        p(i, j) = rot(j, 1)*pr(1) + rot(j, 2)*pr(2) + rot(j, 3)*pr(3)
        v(i, j) = rot(j, 1)*vr(1) + rot(j, 2)*vr(2) + rot(j, 3)*vr(3)
      End Do
    130 End Do
  End If
  If ((dbx**2+dby**2+dbz**2)>1D-20) Then
    db = sqrt(dbx**2+dby**2+dbz**2)
    If (db>0.99999999D0) Then
      Call luerrm(3, '(LUROBO:) boost vector too large')
      dbx = dbx*(0.99999999D0/db)
      dby = dby*(0.99999999D0/db)
      dbz = dbz*(0.99999999D0/db)
      db = 0.99999999D0
    End If
    dga = 1D0/sqrt(1D0-db**2)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 150
      Do j = 1, 4
        dp(j) = dble(p(i,j))
        dv(j) = dble(v(i,j))
      End Do
      dbp = dbx*dp(1) + dby*dp(2) + dbz*dp(3)
      dgabp = dga*(dga*dbp/(1D0+dga)+dp(4))
      p(i, 1) = sngl(dp(1)+dgabp*dbx)
      p(i, 2) = sngl(dp(2)+dgabp*dby)
      p(i, 3) = sngl(dp(3)+dgabp*dbz)
      p(i, 4) = sngl(dga*(dp(4)+dbp))
      dbv = dbx*dv(1) + dby*dv(2) + dbz*dv(3)
      dgabv = dga*(dga*dbv/(1D0+dga)+dv(4))
      v(i, 1) = sngl(dv(1)+dgabv*dbx)
      v(i, 2) = sngl(dv(2)+dgabv*dby)
      v(i, 3) = sngl(dv(3)+dgabv*dbz)
      v(i, 4) = sngl(dga*(dv(4)+dbv))
    150 End Do
  End If
  Return
End Subroutine lurobo
