Subroutine hijels(psc1, psc2)
  Implicit Double Precision (D)
  Dimension psc1(5), psc2(5)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /rndf77/nseed
  Save
  cc = 1.0 - hint1(12)/hint1(13)
  rr = (1.0-cc)*hint1(13)/hint1(12)/(1.0-hipr1(33)) - 1.0
  bb = 0.5*(3.0+rr+sqrt(9.0+10.0*rr+rr**2))
  ep = sqrt((psc1(1)-psc2(1))**2+(psc1(2)-psc2(2))**2+(psc1(3)-psc2(3))**2)
  If (ep<=0.1) Return
  els0 = 98.0/ep + 52.0*(1.0+rr)**2
  pcm1 = psc1(1) + psc2(1)
  pcm2 = psc1(2) + psc2(2)
  pcm3 = psc1(3) + psc2(3)
  ecm = psc1(4) + psc2(4)
  am1 = psc1(5)**2
  am2 = psc2(5)**2
  amm = ecm**2 - pcm1**2 - pcm2**2 - pcm3**2
  If (amm<=psc1(5)+psc2(5)) Return
  pmax = (amm**2+am1**2+am2**2-2.0*amm*am1-2.0*amm*am2-2.0*am1*am2)/4.0/amm
  pmax = abs(pmax)
  20 tt = ranart(nseed)*min(pmax, 1.5)
  els = 98.0*exp(-2.8*tt)/ep + 52.0*exp(-9.2*tt)*(1.0+rr*exp(-4.6*(bb-1.0)*tt))**2
  If (ranart(nseed)>els/els0) Goto 20
  phi = 2.0*hipr1(40)*ranart(nseed)
  dbx = dble(pcm1/ecm)
  dby = dble(pcm2/ecm)
  dbz = dble(pcm3/ecm)
  db = dsqrt(dbx**2+dby**2+dbz**2)
  If (db>0.99999999D0) Then
    dbx = dbx*(0.99999999D0/db)
    dby = dby*(0.99999999D0/db)
    dbz = dbz*(0.99999999D0/db)
    db = 0.99999999D0
    Write (6, *) ' (HIJELS) boost vector too large'
  End If
  dga = 1D0/sqrt(1D0-db**2)
  dp1 = dble(sqrt(tt)*sin(phi))
  dp2 = dble(sqrt(tt)*cos(phi))
  dp3 = dble(sqrt(pmax-tt))
  dp4 = dble(sqrt(pmax+am1))
  dbp = dbx*dp1 + dby*dp2 + dbz*dp3
  dgabp = dga*(dga*dbp/(1D0+dga)+dp4)
  psc1(1) = sngl(dp1+dgabp*dbx)
  psc1(2) = sngl(dp2+dgabp*dby)
  psc1(3) = sngl(dp3+dgabp*dbz)
  psc1(4) = sngl(dga*(dp4+dbp))
  dp1 = -dble(sqrt(tt)*sin(phi))
  dp2 = -dble(sqrt(tt)*cos(phi))
  dp3 = -dble(sqrt(pmax-tt))
  dp4 = dble(sqrt(pmax+am2))
  dbp = dbx*dp1 + dby*dp2 + dbz*dp3
  dgabp = dga*(dga*dbp/(1D0+dga)+dp4)
  psc2(1) = sngl(dp1+dgabp*dbx)
  psc2(2) = sngl(dp2+dgabp*dby)
  psc2(3) = sngl(dp3+dgabp*dbz)
  psc2(4) = sngl(dga*(dp4+dbp))
  Return
End Subroutine hijels
