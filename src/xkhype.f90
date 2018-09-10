Subroutine xkhype(i1, i2, srt, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, aphi=1.02, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (pimass=0.140, ameta=0.5473, aka=0.498, aml=1.116, ams=1.193, am1440=1.44, am1535=1.535)
  Common /ee/id(maxstr), lb(maxstr)
  Save
  s = srt**2
  sigk = 1.E-08
  xky1 = 0.0
  xky2 = 0.0
  xky3 = 0.0
  xky4 = 0.0
  xky5 = 0.0
  xky6 = 0.0
  xky7 = 0.0
  xky8 = 0.0
  xky9 = 0.0
  xky10 = 0.0
  xky11 = 0.0
  xky12 = 0.0
  xky13 = 0.0
  xky14 = 0.0
  xky15 = 0.0
  xky16 = 0.0
  xky17 = 0.0
  lb1 = lb(i1)
  lb2 = lb(i2)
  If (iabs(lb1)==14 .Or. iabs(lb2)==14) Then
    xkaon0 = pnlka(srt)
    xkaon0 = 2.0*xkaon0
    pi2 = (s-(aml+aka)**2)*(s-(aml-aka)**2)
  Else
    xkaon0 = pnska(srt)
    xkaon0 = 2.0*xkaon0
    pi2 = (s-(ams+aka)**2)*(s-(ams-aka)**2)
  End If
  If (pi2<=0.0) Return
  xm1 = pimass
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky1 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = pimass
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky2 = 12.0*pf2/pi2*xkaon0
  End If
  xm1 = pimass
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky3 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = pimass
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky4 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = amrho
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky5 = 9.0*pf2/pi2*xkaon0
  End If
  xm1 = amrho
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky6 = 36.0*pf2/pi2*xkaon0
  End If
  xm1 = amrho
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky7 = 9.0*pf2/pi2*xkaon0
  End If
  xm1 = amrho
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky8 = 9.0*pf2/pi2*xkaon0
  End If
  xm1 = amomga
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky9 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = amomga
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky10 = 12.0*pf2/pi2*xkaon0
  End If
  xm1 = amomga
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky11 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = amomga
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky12 = 3.0*pf2/pi2*xkaon0
  End If
  xm1 = ameta
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky13 = 1.0*pf2/pi2*xkaon0
  End If
  xm1 = ameta
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky14 = 4.0*pf2/pi2*xkaon0
  End If
  xm1 = ameta
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky15 = 1.0*pf2/pi2*xkaon0
  End If
  xm1 = ameta
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky16 = 1.0*pf2/pi2*xkaon0
  End If
  If (lb1==14 .Or. lb2==14) Then
    If (srt>(aphi+amn)) Then
      srrt = srt - (aphi+amn)
      sig = 1.715/((srrt+3.508)**2-12.138)
      xm1 = amn
      xm2 = aphi
      pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
      xky17 = 3.0*pf2/pi2*sig/10.
    End If
  End If
  If ((iabs(lb1)>=15 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=15 .And. iabs(lb2)<=17)) Then
    ddf = 3.0
    xky1 = xky1/ddf
    xky2 = xky2/ddf
    xky3 = xky3/ddf
    xky4 = xky4/ddf
    xky5 = xky5/ddf
    xky6 = xky6/ddf
    xky7 = xky7/ddf
    xky8 = xky8/ddf
    xky9 = xky9/ddf
    xky10 = xky10/ddf
    xky11 = xky11/ddf
    xky12 = xky12/ddf
    xky13 = xky13/ddf
    xky14 = xky14/ddf
    xky15 = xky15/ddf
    xky16 = xky16/ddf
  End If
  sigk = xky1 + xky2 + xky3 + xky4 + xky5 + xky6 + xky7 + xky8 + xky9 + xky10 + xky11 + xky12 + xky13 + xky14 + xky15 + xky16 + xky17
  Return
End Subroutine xkhype
