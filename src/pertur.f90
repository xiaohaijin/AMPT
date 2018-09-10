Subroutine pertur(px, py, pz, srt, irun, i1, i2, nt, kp, icont)
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, aks=0.895)
  Parameter (acas=1.3213, aome=1.6724, amrho=0.769, amomga=0.782)
  Parameter (aeta=0.548, adiomg=3.2288)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /hh/proper(maxstr)
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /input/nstar, ndirct, dir
  Common /nn/nnn
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /pe/propi(maxstr, maxr)
  Common /rr/massr(0:maxr)
  Common /bg/betax, betay, betaz, gamma
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  px0 = px
  py0 = py
  pz0 = pz
  lb1 = lb(i1)
  em1 = e(i1)
  x1 = r(1, i1)
  y1 = r(2, i1)
  z1 = r(3, i1)
  prob1 = proper(i1)
  lb2 = lb(i2)
  em2 = e(i2)
  x2 = r(1, i2)
  y2 = r(2, i2)
  z2 = r(3, i2)
  prob2 = proper(i2)
  icont = 1
  icsbel = -1
  If ((lb1==21 .Or. lb1==23 .Or. iabs(lb1)==30) .And. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 60
  If ((lb2==21 .Or. lb2==23 .Or. iabs(lb2)==30) .And. (iabs(lb1)>=14 .And. iabs(lb1)<=17)) Goto 60
  If ((lb1==21 .Or. lb1==23 .Or. iabs(lb1)==30) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) Goto 70
  If ((lb2==21 .Or. lb2==23 .Or. iabs(lb2)==30) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41)) Goto 70
  If ((((lb1>=3 .And. lb1<=5) .Or. lb1==0) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) .Or. (((lb2>=3 .And. lb2<=5) .Or. lb2==0) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41))) Goto 90
  If (((lb1>=3 .And. lb1<=5) .And. iabs(lb2)==45) .Or. ((lb2>=3 .And. lb2<=5) .And. iabs(lb1)==45)) Goto 110
  60 If (iabs(lb1)>=14 .And. iabs(lb1)<=17) Then
    asap = e(i1)
    akap = e(i2)
    idp = i1
  Else
    asap = e(i2)
    akap = e(i1)
    idp = i2
  End If
  app = 0.138
  If (srt<(acas+app)) Return
  srrt = srt - (acas+app) + (amn+akap)
  pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  pii = sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))
  pff = sqrt((srt**2-(asap+app)**2)*(srt**2-(asap-app)**2))
  cmat = sigca*pii/pff
  sigpi = cmat*sqrt((srt**2-(acas+app)**2)*(srt**2-(acas-app)**2))/sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
  sigeta = 0.
  If (srt>(acas+aeta)) Then
    srrt = srt - (acas+aeta) + (amn+akap)
    pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
    sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
    cmat = sigca*pii/pff
    sigeta = cmat*sqrt((srt**2-(acas+aeta)**2)*(srt**2-(acas-aeta)**2))/sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
  End If
  sigca = sigpi + sigeta
  sigpe = 0.
  sig = amax1(sigpe, sigca)
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  brpp = sigca/sig
  If ((lb1>=14 .And. lb1<=17) .Or. (lb2>=14 .And. lb2<=17)) Then
    lbpp1 = 40 + int(2*ranart(nseed))
  Else
    lbpp1 = -40 - int(2*ranart(nseed))
  End If
  empp1 = acas
  If (ranart(nseed)<sigpi/sigca) Then
    lbpp2 = 3 + int(3*ranart(nseed))
    empp2 = 0.138
  Else
    lbpp2 = 0
    empp2 = aeta
  End If
  If (ranart(nseed)<brpp) Then
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
    proper(i1) = brpp
    lb(i2) = lbpp2
    e(i2) = empp2
    proper(i2) = 1.
  End If
  Goto 700
  70 If (iabs(lb1)==40 .Or. iabs(lb1)==41) Then
    acap = e(i1)
    akap = e(i2)
    idp = i1
  Else
    acap = e(i2)
    akap = e(i1)
    idp = i2
  End If
  app = 0.138
  ames = 0.138
  If (srt<(aome+ames)) Return
  srrt = srt - (aome+ames) + (amn+akap)
  pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
  sigomm = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigomm*sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))/sqrt((srt**2-(asa+app)**2)*(srt**2-(asa-app)**2))
  sigom = cmat*sqrt((srt**2-(aome+ames)**2)*(srt**2-(aome-ames)**2))/sqrt((srt**2-(acap+akap)**2)*(srt**2-(acap-akap)**2))
  sigpe = 0.
  sig = amax1(sigpe, sigom)
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  brpp = sigom/sig
  If ((lb1>=40 .And. lb1<=41) .Or. (lb2>=40 .And. lb2<=41)) Then
    lbpp1 = 45
  Else
    lbpp1 = -45
  End If
  empp1 = aome
  lbpp2 = 3 + int(3*ranart(nseed))
  empp2 = ames
  xrand = ranart(nseed)
  If (xrand<(proper(idp)*brpp)) Then
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
    proper(i1) = proper(idp)*brpp
    lb(i2) = lbpp2
    e(i2) = empp2
    proper(i2) = 1.
  Else If (xrand<brpp) Then
    e(idp) = 0.
  End If
  Goto 700
  90 If (iabs(lb1)==40 .Or. iabs(lb1)==41) Then
    acap = e(i1)
    app = e(i2)
    idp = i1
    idn = i2
  Else
    acap = e(i2)
    app = e(i1)
    idp = i2
    idn = i1
  End If
  akal = aka
  alas = ala
  If (srt<=(alas+aka)) Return
  srrt = srt - (acap+app) + (amn+aka)
  pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
  sigca = cmat*sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
  dfr = 1./3.
  If (lb(idn)==0) dfr = 1.
  sigcal = sigca*dfr*(srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/(srt**2-(acap-app)**2)
  alas = asa
  If (srt<=(alas+aka)) Then
    sigcas = 0.
  Else
    srrt = srt - (acap+app) + (amn+aka)
    pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
    sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
    cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
    sigca = cmat*sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
    dfr = 1.
    If (lb(idn)==0) dfr = 3.
    sigcas = sigca*dfr*(srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/(srt**2-(acap-app)**2)
  End If
  sig = sigcal + sigcas
  brpp = 1.
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Then
    ds = sqrt(20.0/31.4)
    dsr = ds + 0.1
    Call distce(i1, i2, dsr, ds, dt, ec, srt, icsbel, px, py, pz)
    If (icsbel==-1) Return
    empp1 = em1
    empp2 = em2
    Goto 700
  End If
  If (sigcal/sig>ranart(nseed)) Then
    If (lb1==40 .Or. lb1==41 .Or. lb2==40 .Or. lb2==41) Then
      lbpp1 = 21
      lbpp2 = 14
    Else
      lbpp1 = 23
      lbpp2 = -14
    End If
    alas = ala
  Else
    If (lb1==40 .Or. lb1==41 .Or. lb2==40 .Or. lb2==41) Then
      lbpp1 = 21
      lbpp2 = 15 + int(3*ranart(nseed))
    Else
      lbpp1 = 23
      lbpp2 = -15 - int(3*ranart(nseed))
    End If
    alas = asa
  End If
  empp1 = aka
  empp2 = alas
  If (ranart(nseed)<proper(idp)) Then
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
    proper(i1) = 1.
    lb(i2) = lbpp2
    e(i2) = empp2
    proper(i2) = 1.
    Goto 700
  Else
    e(idp) = 0.
  End If
  Return
  110 If (lb1==45 .Or. lb1==-45) Then
    aomp = e(i1)
    app = e(i2)
    idp = i1
    idn = i2
  Else
    aomp = e(i2)
    app = e(i1)
    idp = i2
    idn = i1
  End If
  akal = aka
  If (srt<=(acas+aka)) Return
  srrt = srt - (aome+app) + (amn+aka)
  pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(asa+0.138)**2)*(srt**2-(asa-0.138)**2))
  sigom = cmat*sqrt((srt**2-(aomp+app)**2)*(srt**2-(aomp-app)**2))/sqrt((srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2))
  dfr = 2./3.
  sigom = sigom*dfr*(srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2)/(srt**2-(aomp+app)**2)/(srt**2-(aomp-app)**2)
  brpp = 1.
  ds = sqrt(sigom/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Then
    ds = sqrt(20.0/31.4)
    dsr = ds + 0.1
    Call distce(i1, i2, dsr, ds, dt, ec, srt, icsbel, px, py, pz)
    If (icsbel==-1) Return
    empp1 = em1
    empp2 = em2
    Goto 700
  End If
  If (lb1==45 .Or. lb2==45) Then
    lbpp1 = 40 + int(2*ranart(nseed))
    lbpp2 = 21
  Else
    lbpp1 = -40 - int(2*ranart(nseed))
    lbpp2 = 23
  End If
  empp1 = acas
  empp2 = aka
  If (ranart(nseed)<proper(idp)) Then
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
    proper(i1) = proper(idp)
    lb(i2) = lbpp2
    e(i2) = empp2
    proper(i2) = 1.
  Else
    e(idp) = 0.
  End If
  Goto 700
  700 Continue
  pr2 = (srt**2-empp1**2-empp2**2)**2 - 4.0*(empp1*empp2)**2
  If (pr2<=0.) pr2 = 0.00000001
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Call rotate(px0, py0, pz0, px, py, pz)
  If (icont==0) Return
  e1cm = sqrt(empp1**2+px**2+py**2+pz**2)
  p1beta = px*betax + py*betay + pz*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  ppt11 = betax*transf + px
  ppt12 = betay*transf + py
  ppt13 = betaz*transf + pz
  If (icsbel/=-1) Then
    p(1, i1) = ppt11
    p(2, i1) = ppt12
    p(3, i1) = ppt13
    e2cm = sqrt(empp2**2+px**2+py**2+pz**2)
    transf = gamma*(-gamma*p1beta/(gamma+1)+e2cm)
    ppt21 = betax*transf - px
    ppt22 = betay*transf - py
    ppt23 = betaz*transf - pz
    p(1, i2) = ppt21
    p(2, i2) = ppt22
    p(3, i2) = ppt23
    Return
  End If
  xpt = x1
  ypt = y1
  zpt = z1
  nnn = nnn + 1
  propi(nnn, irun) = proper(idp)*brpp
  lpion(nnn, irun) = lbpp1
  epion(nnn, irun) = empp1
  rpion(1, nnn, irun) = xpt
  rpion(2, nnn, irun) = ypt
  rpion(3, nnn, irun) = zpt
  ppion(1, nnn, irun) = ppt11
  ppion(2, nnn, irun) = ppt12
  ppion(3, nnn, irun) = ppt13
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  Return
End Subroutine pertur
