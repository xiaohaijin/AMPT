Subroutine atrobo(the, phi, bex, bey, bez, imin, imax, ierror)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Dimension rot(3, 3), pv(3)
  Double Precision dp(4), dbex, dbey, dbez, dga, dga2, dbep, dgabep
  Save
  ierror = 0
  If (imin<=0 .Or. imax>n .Or. imin>imax) Return
  If (the**2+phi**2>1E-20) Then
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
      Do j = 1, 3
        pv(j) = p(i, j)
      End Do
      Do j = 1, 3
        p(i, j) = rot(j, 1)*pv(1) + rot(j, 2)*pv(2) + rot(j, 3)*pv(3)
      End Do
    End Do
  End If
  If (bex**2+bey**2+bez**2>1E-20) Then
    dbex = dble(bex)
    dbey = dble(bey)
    dbez = dble(bez)
    dga2 = 1D0 - dbex**2 - dbey**2 - dbez**2
    If (dga2<=0D0) Then
      ierror = 1
      Return
    End If
    dga = 1D0/dsqrt(dga2)
    Do i = imin, imax
      Do j = 1, 4
        dp(j) = dble(p(i,j))
      End Do
      dbep = dbex*dp(1) + dbey*dp(2) + dbez*dp(3)
      dgabep = dga*(dga*dbep/(1D0+dga)+dp(4))
      p(i, 1) = sngl(dp(1)+dgabep*dbex)
      p(i, 2) = sngl(dp(2)+dgabep*dbey)
      p(i, 3) = sngl(dp(3)+dgabep*dbez)
      p(i, 4) = sngl(dga*(dp(4)+dbep))
    End Do
  End If
  Return
End Subroutine atrobo
