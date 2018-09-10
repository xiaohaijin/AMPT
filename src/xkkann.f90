Subroutine xkkann(srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigk, rrkk)
  Parameter (maxstr=150001, maxx=20, maxz=24)
  Parameter (aka=0.498, pimass=0.140, rhom=0.770, omegam=0.7819, etam=0.5473, aphi=1.02)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Save
  s = srt**2
  sigk = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  xsk7 = 0.0
  xsk8 = 0.0
  xsk9 = 0.0
  xsk10 = 0.0
  xsk11 = 0.0
  xpion0 = pipik(srt)
  xpion0 = 2.0*xpion0
  pi2 = s*(s-4.0*aka**2)
  If (pi2<=0.0) Return
  xm1 = pimass
  xm2 = pimass
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk1 = 9.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = pimass
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk4 = 3.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = etam
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk10 = 1.0/4.0*pf2/pi2*xpion0
  End If
  xpion0 = rrkk
  xm1 = rhom
  xm2 = rhom
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk5 = 81.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = rhom
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk6 = 27.0/4.0*pf2/pi2*xpion0
  End If
  xm1 = omegam
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk8 = 9.0/4.0*pf2/pi2*xpion0
  End If
  fwdp = 1.68*(aphi**2-4.*aka**2)**1.5/6./aphi/aphi
  scheck = srt**2 - 4.0*aka**2
  If (scheck<=0) Then
    Write (99, *) 'scheck47: ', scheck
    Stop
  End If
  pkaon = 0.5*sqrt(scheck)
  xsk11 = 30.*3.14159*0.1973**2*(aphi*fwdp)**2/((srt**2-aphi**2)**2+(aphi*fwdp)**2)/pkaon**2
  sigk = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6 + xsk7 + xsk8 + xsk9 + xsk10 + xsk11
  Return
End Subroutine xkkann
