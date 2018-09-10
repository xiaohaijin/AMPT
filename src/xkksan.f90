Subroutine xkksan(i1, i2, srt, sigks1, sigks2, sigks3, sigks4, sigk, prkk)
  Parameter (aka=0.498, pimass=0.140, rhom=0.770, aks=0.895, omegam=0.7819, etam=0.5473)
  Parameter (maxstr=150001)
  Common /cc/e(maxstr)
  Save
  s = srt**2
  sigks1 = 1.E-08
  sigks2 = 1.E-08
  sigks3 = 1.E-08
  sigks4 = 1.E-08
  xpion0 = prkk
  xpion0 = xpion0/2
  pi2 = (s-(e(i1)+e(i2))**2)*(s-(e(i1)-e(i2))**2)
  sigk = 1.E-08
  If (pi2<=0.0) Return
  xm1 = pimass
  xm2 = rhom
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pi2>0.0 .And. pf2>0.0) Then
    sigks1 = 27.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = pimass
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pi2>0.0 .And. pf2>0.0) Then
    sigks2 = 9.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = rhom
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    sigks3 = 9.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = omegam
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    sigks4 = 3.0/4.0*pf2/pi2*xpion0
  End If
  sigk = sigks1 + sigks2 + sigks3 + sigks4
  Return
End Subroutine xkksan
