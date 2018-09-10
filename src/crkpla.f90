Subroutine crkpla(px, py, pz, ec, srt, spika, emm1, emm2, lbp1, lbp2, i1, i2, icase, srhoks)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, amrho=0.769, amomga=0.782, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895, ala=1.1157, asa=1.1974, aphi=1.02)
  Parameter (am1440=1.44, am1535=1.535)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  emm1 = 0.
  emm2 = 0.
  lbp1 = 0
  lbp2 = 0
  xkp0 = spika
  xkp1 = 0.
  xkp2 = 0.
  xkp3 = 0.
  xkp4 = 0.
  xkp5 = 0.
  xkp6 = 0.
  xkp7 = 0.
  xkp8 = 0.
  xkp9 = 0.
  xkp10 = 0.
  sigm = 15.
  pdd = (srt**2-(aka+ap1)**2)*(srt**2-(aka-ap1)**2)
  If (srt<(ala+amn)) Goto 70
  xkp1 = sigm*(4./3.)*(srt**2-(ala+amn)**2)*(srt**2-(ala-amn)**2)/pdd
  If (srt>(ala+am0)) Then
    xkp2 = sigm*(16./3.)*(srt**2-(ala+am0)**2)*(srt**2-(ala-am0)**2)/pdd
  End If
  If (srt>(ala+am1440)) Then
    xkp3 = sigm*(4./3.)*(srt**2-(ala+am1440)**2)*(srt**2-(ala-am1440)**2)/pdd
  End If
  If (srt>(ala+am1535)) Then
    xkp4 = sigm*(4./3.)*(srt**2-(ala+am1535)**2)*(srt**2-(ala-am1535)**2)/pdd
  End If
  If (srt>(asa+amn)) Then
    xkp5 = sigm*4.*(srt**2-(asa+amn)**2)*(srt**2-(asa-amn)**2)/pdd
  End If
  If (srt>(asa+am0)) Then
    xkp6 = sigm*16.*(srt**2-(asa+am0)**2)*(srt**2-(asa-am0)**2)/pdd
  End If
  If (srt>(asa+am1440)) Then
    xkp7 = sigm*4.*(srt**2-(asa+am1440)**2)*(srt**2-(asa-am1440)**2)/pdd
  End If
  If (srt>(asa+am1535)) Then
    xkp8 = sigm*4.*(srt**2-(asa+am1535)**2)*(srt**2-(asa-am1535)**2)/pdd
  End If
  70 Continue
  sig1 = 195.639
  sig2 = 372.378
  If (srt>aphi+aka) Then
    pff = sqrt((srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2))
    scheck = pdd
    If (scheck<=0) Then
      Write (99, *) 'scheck40: ', scheck
      Stop
    End If
    xkp9 = sig1*pff/sqrt(pdd)*1./32./pi/srt**2
    If (srt>aphi+aks) Then
      pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))
      scheck = pdd
      If (scheck<=0) Then
        Write (99, *) 'scheck41: ', scheck
        Stop
      End If
      xkp10 = sig2*pff/sqrt(pdd)*3./32./pi/srt**2
    End If
  End If
  sigpik = 0.
  If (srt>(amrho+aks)) Then
    sigpik = srhoks*9.*(srt**2-(0.77-aks)**2)*(srt**2-(0.77+aks)**2)/4/srt**2/(px**2+py**2+pz**2)
    If (srt>(amomga+aks)) sigpik = sigpik*12./9.
  End If
  sigkp = xkp0 + xkp1 + xkp2 + xkp3 + xkp4 + xkp5 + xkp6 + xkp7 + xkp8 + xkp9 + xkp10 + sigpik
  icase = 0
  dskn = sqrt(sigkp/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  randu = ranart(nseed)*sigkp
  xkp1 = xkp0 + xkp1
  xkp2 = xkp1 + xkp2
  xkp3 = xkp2 + xkp3
  xkp4 = xkp3 + xkp4
  xkp5 = xkp4 + xkp5
  xkp6 = xkp5 + xkp6
  xkp7 = xkp6 + xkp7
  xkp8 = xkp7 + xkp8
  xkp9 = xkp8 + xkp9
  xkp10 = xkp9 + xkp10
  If (randu<=xkp0) Then
    icase = 1
    Return
  Else
    icase = 2
    If (randu<=xkp1) Then
      lbp1 = -14
      lbp2 = 1 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = amn
      Goto 60
    Else If (randu<=xkp2) Then
      lbp1 = -14
      lbp2 = 6 + int(4*ranart(nseed))
      emm1 = ala
      emm2 = am0
      Goto 60
    Else If (randu<=xkp3) Then
      lbp1 = -14
      lbp2 = 10 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = am1440
      Goto 60
    Else If (randu<=xkp4) Then
      lbp1 = -14
      lbp2 = 12 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = am1535
      Goto 60
    Else If (randu<=xkp5) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 1 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = amn
      Goto 60
    Else If (randu<=xkp6) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 6 + int(4*ranart(nseed))
      emm1 = asa
      emm2 = am0
      Goto 60
    Else If (randu<xkp7) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 10 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = am1440
      Goto 60
    Else If (randu<xkp8) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 12 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = am1535
      Goto 60
    Else If (randu<xkp9) Then
      icase = 3
      lbp1 = 29
      lbp2 = 23
      emm1 = aphi
      emm2 = aka
      If (lb(i1)==21 .Or. lb(i2)==21) Then
        lbp2 = 21
        icase = -3
      End If
      Goto 60
    Else If (randu<xkp10) Then
      icase = 4
      lbp1 = 29
      lbp2 = 30
      emm1 = aphi
      emm2 = aks
      If (lb(i1)==21 .Or. lb(i2)==21) Then
        lbp2 = -30
        icase = -4
      End If
      Goto 60
    Else
      icase = 5
      lbp1 = 25 + int(3*ranart(nseed))
      lbp2 = 30
      emm1 = amrho
      emm2 = aks
      If (srt>(amomga+aks) .And. ranart(nseed)<0.25) Then
        lbp1 = 28
        emm1 = amomga
      End If
      If (lb(i1)==21 .Or. lb(i2)==21) Then
        lbp2 = -30
        icase = -5
      End If
    End If
  End If
  60 If (icase==2 .And. (lb(i1)==21 .Or. lb(i2)==21)) Then
    lbp1 = -lbp1
    lbp2 = -lbp2
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
End Subroutine crkpla
