Subroutine crkphi(px, py, pz, ec, srt, iblock, emm1, emm2, lbp1, lbp2, i1, i2, ikk, icase, rrkk, prkk)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, aphi=1.02, am0=1.232, amns=1.52, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, acas=1.3213)
  Parameter (aks=0.895, aomega=0.7819, arho=0.77)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  lb1 = lb(i1)
  lb2 = lb(i2)
  icase = 0
  If (srt<(aphi+ap1)) Then
    sig1 = 0.
    sig2 = 0.
    sig3 = 0.
  Else
    If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
      dnr = 4.
      ikk = 2
    Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
      dnr = 12.
      ikk = 1
    Else
      dnr = 36.
      ikk = 0
    End If
    sig1 = 0.
    sig2 = 0.
    sig3 = 0.
    srri = e(i1) + e(i2)
    srr1 = aphi + ap1
    srr2 = aphi + aomega
    srr3 = aphi + arho
    pii = (srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
    srrt = srt - amax1(srri, srr1)
    If (srrt<0.3 .And. srrt>0.01) Then
      sig = 1.69/(srrt**0.141-0.407)
    Else
      sig = 3.74 + 0.008*srrt**1.9
    End If
    sig1 = sig*(9./dnr)*(srt**2-(aphi+ap1)**2)*(srt**2-(aphi-ap1)**2)/pii
    If (srt>aphi+aomega) Then
      srrt = srt - amax1(srri, srr2)
      If (srrt<0.3 .And. srrt>0.01) Then
        sig = 1.69/(srrt**0.141-0.407)
      Else
        sig = 3.74 + 0.008*srrt**1.9
      End If
      sig2 = sig*(9./dnr)*(srt**2-(aphi+aomega)**2)*(srt**2-(aphi-aomega)**2)/pii
    End If
    If (srt>aphi+arho) Then
      srrt = srt - amax1(srri, srr3)
      If (srrt<0.3 .And. srrt>0.01) Then
        sig = 1.69/(srrt**0.141-0.407)
      Else
        sig = 3.74 + 0.008*srrt**1.9
      End If
      sig3 = sig*(27./dnr)*(srt**2-(aphi+arho)**2)*(srt**2-(aphi-arho)**2)/pii
    End If
  End If
  rrkk0 = rrkk
  prkk0 = prkk
  sigm = 0.
  If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
    Call xkkann(srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigm, rrkk0)
  Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
    Call xkksan(i1, i2, srt, sigks1, sigks2, sigks3, sigks4, sigm, prkk0)
  Else
  End If
  sigm0 = sigm
  sigks = sig1 + sig2 + sig3 + sigm
  dskn = sqrt(sigks/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  icase = 1
  ranx = ranart(nseed)
  lbp1 = 29
  emm1 = aphi
  If (ranx<=sig1/sigks) Then
    lbp2 = 3 + int(3*ranart(nseed))
    emm2 = ap1
  Else If (ranx<=(sig1+sig2)/sigks) Then
    lbp2 = 28
    emm2 = aomega
  Else If (ranx<=(sig1+sig2+sig3)/sigks) Then
    lbp2 = 25 + int(3*ranart(nseed))
    emm2 = arho
  Else
    If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
      Call crkkpi(i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigm0, iblock, lbp1, lbp2, emm1, emm2)
    Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
      Call crkspi(i1, i2, sigks1, sigks2, sigks3, sigks4, sigm0, iblock, lbp1, lbp2, emm1, emm2)
    Else
    End If
  End If
  px0 = px
  py0 = py
  pz0 = pz
  pr2 = (srt**2-emm1**2-emm2**2)**2 - 4.0*(emm1*emm2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Return
End Subroutine crkphi
