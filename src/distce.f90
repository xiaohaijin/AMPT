Subroutine distce(i1, i2, deltar, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
  Parameter (maxstr=150001)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /ee/id(maxstr), lb(maxstr)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Save
  ic = 0
  x1 = r(1, i1)
  y1 = r(2, i1)
  z1 = r(3, i1)
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  x2 = r(1, i2)
  y2 = r(2, i2)
  z2 = r(3, i2)
  px2 = p(1, i2)
  py2 = p(2, i2)
  pz2 = p(3, i2)
  em1 = e(i1)
  em2 = e(i2)
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  rsqare = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
  If (rsqare>deltar**2) Goto 400
  e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
  s = srt*srt
  If (s<ec) Goto 400
  p1beta = px1*betax + py1*betay + pz1*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)-e1)
  prcm = sqrt(px1cm**2+py1cm**2+pz1cm**2)
  If (prcm<=0.00001) Goto 400
  drbeta = betax*(x1-x2) + betay*(y1-y2) + betaz*(z1-z2)
  transf = gamma*gamma*drbeta/(gamma+1)
  dxcm = betax*transf + x1 - x2
  dycm = betay*transf + y1 - y2
  dzcm = betaz*transf + z1 - z2
  drcm = sqrt(dxcm**2+dycm**2+dzcm**2)
  dzz = (px1cm*dxcm+py1cm*dycm+pz1cm*dzcm)/prcm
  If ((drcm**2-dzz**2)<=0.) Then
     bbb = 0.
  Else
     bbb = sqrt(drcm**2-dzz**2)
  End If
  If (bbb>ds) Goto 400
  relvel = prcm*(1.0/e1+1.0/e2)
  ddd = relvel*dt*0.5
  If (abs(ddd)<abs(dzz)) Goto 400
  ic = 1
  Goto 500
400 ic = -1
500 Continue
  Return
End Subroutine distce
