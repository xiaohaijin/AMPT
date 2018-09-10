Subroutine pystfe(kf, x, q2, xpq)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Dimension xpq(-6:6), xfdflm(9)
  Character chdflm(9)*5, header*40
  Data chdflm/'UPVAL', 'DOVAL', 'GLUON', 'QBAR ', 'UBAR ', 'SBAR ', 'CBAR ', 'BBAR ', 'TBAR '/
  Data header/'Tung evolution package has been invoked'/
  Data init/0/
  If (mstp(51)>=11 .And. mstp(51)<=13 .And. mstp(52)<=1) Then
    xdflm = max(0.51E-4, x)
    q2dflm = max(10., min(1E8,q2))
    If (mstp(52)==0) q2dflm = 10.
    Do j = 1, 9
      If (mstp(52)==1 .And. j==9) Then
        q2dflm = q2dflm*(40./pmas(6,1))**2
        q2dflm = max(10., min(1E8,q2))
      End If
      xfdflm(j) = 0.
    End Do
    If (x<0.51E-4 .And. abs(parp(51)-1.)>0.01) Then
      cxs = (0.51E-4/x)**(parp(51)-1.)
      Do j = 1, 7
        xfdflm(j) = xfdflm(j)*cxs
      End Do
    End If
    xpq(0) = xfdflm(3)
    xpq(1) = xfdflm(2) + xfdflm(5)
    xpq(2) = xfdflm(1) + xfdflm(5)
    xpq(3) = xfdflm(6)
    xpq(4) = xfdflm(7)
    xpq(5) = xfdflm(8)
    xpq(6) = xfdflm(9)
    xpq(-1) = xfdflm(5)
    xpq(-2) = xfdflm(5)
    xpq(-3) = xfdflm(6)
    xpq(-4) = xfdflm(7)
    xpq(-5) = xfdflm(8)
    xpq(-6) = xfdflm(9)
  Else
    If (init==0) Then
      i1 = 0
      If (mstp(52)==4) i1 = 1
      ihdrn = 1
      nu = mstp(53)
      i2 = mstp(51)
      If (mstp(51)>=11) i2 = mstp(51) - 3
      i3 = 0
      If (mstp(52)==3) i3 = 1
      alam = 0.75*parp(1)
      tpms = pmas(6, 1)
      qini = parp(52)
      qmax = parp(53)
      xmin = parp(54)
      init = 1
    End If
    q = sqrt(q2)
    Do i = -6, 6
      fixq = 0.
      xpq(i) = x*fixq
    End Do
    xps = xpq(1)
    xpq(1) = xpq(2)
    xpq(2) = xps
    xps = xpq(-1)
    xpq(-1) = xpq(-2)
    xpq(-2) = xps
  End If
  Return
End Subroutine pystfe
