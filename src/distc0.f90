Subroutine distc0(drmax, deltr0, dt, ifirst, px1cm, py1cm, pz1cm, x1, y1, z1, px1, py1, pz1, em1, x2, y2, z2, px2, py2, pz2, em2)
  Common /bg/betax, betay, betaz, gamma
  Save
  ifirst = -1
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
  p1beta = px1*betax + py1*betay + pz1*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)-e1)
  prcm = sqrt(px1cm**2+py1cm**2+pz1cm**2)
  If (prcm<=0.00001) Return
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
  If (bbb>drmax) Return
  relvel = prcm*(1.0/e1+1.0/e2)
  ddd = relvel*dt*0.5
  If (abs(ddd)<abs(dzz)) Return
  ifirst = 1
  Return
End Subroutine distc0
