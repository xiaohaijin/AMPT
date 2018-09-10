Subroutine crksph(px, py, pz, ec, srt, emm1, emm2, lbp1, lbp2, i1, i2, ikkg, ikkl, iblock, icase, srhoks)
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
  sigela = 10.
  sigkm = 0.
  If ((lb1>=25 .And. lb1<=28) .Or. (lb2>=25 .And. lb2<=28)) Then
    If (iabs(lb1)==30 .Or. iabs(lb2)==30) Then
      sigkm = srhoks
    Else If ((lb1==23 .Or. lb1==21 .Or. lb2==23 .Or. lb2==21) .And. srt>(ap2+aks)) Then
      sigkm = srhoks
    End If
  End If
  If (srt<(aphi+aka)) Then
    sig11 = 0.
    sig22 = 0.
  Else
    If ((iabs(lb1)==30 .And. (lb2>=3 .And. lb2<=5)) .Or. (iabs(lb2)==30 .And. (lb1>=3 .And. lb1<=5))) Then
      dnr = 18.
      ikkl = 0
      iblock = 225
      sig1 = 2047.042
      sig2 = 1496.692
    Else If ((lb1==23 .Or. lb1==21 .And. (lb2>=25 .And. lb2<=27)) .Or. (lb2==23 .Or. lb2==21 .And. (lb1>=25 .And. lb1<=27))) Then
      dnr = 18.
      ikkl = 1
      iblock = 224
      sig1 = 526.702
      sig2 = 1313.960
    Else If ((iabs(lb1)==30 .And. (lb2>=25 .And. lb2<=27)) .Or. (iabs(lb2)==30 .And. (lb1>=25 .And. lb1<=27))) Then
      dnr = 54.
      ikkl = 0
      iblock = 225
      sig1 = 1371.257
      sig2 = 6999.840
    Else If (((lb1==23 .Or. lb1==21) .And. lb2==28) .Or. ((lb2==23 .Or. lb2==21) .And. lb1==28)) Then
      dnr = 6.
      ikkl = 1
      iblock = 224
      sig1 = 355.429
      sig2 = 440.558
    Else
      dnr = 18.
      ikkl = 0
      iblock = 225
      sig1 = 482.292
      sig2 = 1698.903
    End If
    sig11 = 0.
    sig22 = 0.
    scheck = (srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
    If (scheck<=0) Then
      Write (99, *) 'scheck42: ', scheck
      Stop
    End If
    pii = sqrt(scheck)
    scheck = (srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2)
    If (scheck<0) Then
      Write (99, *) 'scheck43: ', scheck
      scheck = 0.
    End If
    pff = sqrt(scheck)
    sig11 = sig1*pff/pii*6./dnr/32./pi/srt**2
    If (srt>aphi+aks) Then
      pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))
      sig22 = sig2*pff/pii*18./dnr/32./pi/srt**2
    End If
  End If
  sigks = sig11 + sig22 + sigela + sigkm
  dskn = sqrt(sigks/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  icase = 1
  ranx = ranart(nseed)
  If (ranx<=(sigela/sigks)) Then
    lbp1 = lb1
    emm1 = e(i1)
    lbp2 = lb2
    emm2 = e(i2)
    iblock = 111
  Else If (ranx<=((sigela+sigkm)/sigks)) Then
    lbp1 = 3 + int(3*ranart(nseed))
    emm1 = 0.14
    If (lb1==23 .Or. lb2==23) Then
      lbp2 = 30
      emm2 = aks
    Else If (lb1==21 .Or. lb2==21) Then
      lbp2 = -30
      emm2 = aks
    Else If (lb1==30 .Or. lb2==30) Then
      lbp2 = 23
      emm2 = aka
    Else
      lbp2 = 21
      emm2 = aka
    End If
    iblock = 112
  Else If (ranx<=((sigela+sigkm+sig11)/sigks)) Then
    lbp2 = 23
    emm2 = aka
    ikkg = 1
    If (lb1==21 .Or. lb2==21 .Or. lb1==-30 .Or. lb2==-30) Then
      lbp2 = 21
      iblock = iblock - 100
    End If
    lbp1 = 29
    emm1 = aphi
  Else
    lbp2 = 30
    emm2 = aks
    ikkg = 0
    iblock = iblock + 2
    If (lb1==21 .Or. lb2==21 .Or. lb1==-30 .Or. lb2==-30) Then
      lbp2 = -30
      iblock = iblock - 100
    End If
    lbp1 = 29
    emm1 = aphi
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
End Subroutine crksph
