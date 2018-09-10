Subroutine rmasdd(srt, am10, am20, dmin1, dmin2, iseed, ic, dm1, dm2)
  Common /rndf77/nseed
  Save
  amn = 0.94
  amp = 0.14
  dmax1 = srt - dmin2
  5 ntry1 = 0
  ntry2 = 0
  ntry = 0
  ictrl = 0
  10 dm1 = ranart(nseed)*(dmax1-dmin1) + dmin1
  ntry1 = ntry1 + 1
  If (ictrl==0) dmax2 = srt - dm1
  20 dm2 = ranart(nseed)*(dmax2-dmin2) + dmin2
  ntry2 = ntry2 + 1
  q2 = ((srt**2-dm1**2-dm2**2)**2-4.*dm1**2*dm2**2)
  If (q2<=0) Then
    dmax2 = dm2 - 0.01
    ictrl = 1
    Goto 20
  End If
  If (dmax1<am10) Then
    If (ic==1) fm1 = fmassd(dmax1)
    If (ic==2) fm1 = fmassn(dmax1)
    If (ic==3) fm1 = fmassd(dmax1)
    If (ic==4) fm1 = fmassd(dmax1)
  Else
    If (ic==1) fm1 = fmassd(am10)
    If (ic==2) fm1 = fmassn(am10)
    If (ic==3) fm1 = fmassd(am10)
    If (ic==4) fm1 = fmassd(am10)
  End If
  If (dmax2<am20) Then
    If (ic==1) fm2 = fmassd(dmax2)
    If (ic==2) fm2 = fmassn(dmax2)
    If (ic==3) fm2 = fmassn(dmax2)
    If (ic==4) fm2 = fmassr(dmax2)
  Else
    If (ic==1) fm2 = fmassd(am20)
    If (ic==2) fm2 = fmassn(am20)
    If (ic==3) fm2 = fmassn(am20)
    If (ic==4) fm2 = fmassr(am20)
  End If
  If (fm1==0.) fm1 = 1.E-04
  If (fm2==0.) fm2 = 1.E-04
  prob0 = fm1*fm2
  If (ic==1) prob = fmassd(dm1)*fmassd(dm2)
  If (ic==2) prob = fmassn(dm1)*fmassn(dm2)
  If (ic==3) prob = fmassd(dm1)*fmassn(dm2)
  If (ic==4) prob = fmassd(dm1)*fmassr(dm2)
  If (prob<=1.E-06) prob = 1.E-06
  fff = prob/prob0
  ntry = ntry + 1
  If (ranart(nseed)>fff .And. ntry<=20) Goto 10
  If ((abs(am10-0.77)<=0.01 .And. dm1>1.07) .Or. (abs(am10-1.232)<=0.01 .And. dm1>1.47) .Or. (abs(am10-1.44)<=0.01 .And. dm1>2.14)) Goto 5
  If ((abs(am20-0.77)<=0.01 .And. dm2>1.07) .Or. (abs(am20-1.232)<=0.01 .And. dm2>1.47) .Or. (abs(am20-1.44)<=0.01 .And. dm2>2.14)) Goto 5
  Return
End Subroutine rmasdd
