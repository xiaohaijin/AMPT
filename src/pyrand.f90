Subroutine pyrand
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  mint(17) = 0
  mint(18) = 0
  vint(143) = 1.
  vint(144) = 1.
  If (msub(95)==1 .Or. mint(82)>=2) Call pymult(2)
  isub = 0
  100 mint(51) = 0
  If (mint(82)==1 .And. (isub<=90 .Or. isub>96)) Then
    rsub = xsec(0, 1)*rlu(0)
    Do i = 1, 200
      If (msub(i)/=1) Goto 110
      isub = i
      rsub = rsub - xsec(i, 1)
      If (rsub<=0.) Goto 120
    110 End Do
    120 If (isub==95) isub = 96
  Else If (mint(82)>=2 .And. isub==0) Then
    rsub = vint(131)*rlu(0)
    isub = 96
    If (rsub>vint(106)) isub = 93
    If (rsub>vint(106)+vint(104)) isub = 92
    If (rsub>vint(106)+vint(104)+vint(103)) isub = 91
  End If
  If (mint(82)==1) ngen(0, 1) = ngen(0, 1) + 1
  If (mint(82)==1) ngen(isub, 1) = ngen(isub, 1) + 1
  mint(1) = isub
  mint(72) = 0
  kfr1 = 0
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    kfr1 = kfpr(isub, 1)
  Else If (isub>=71 .And. isub<=77) Then
    kfr1 = 25
  End If
  If (kfr1/=0) Then
    taur1 = pmas(kfr1, 1)**2/vint(2)
    gamr1 = pmas(kfr1, 1)*pmas(kfr1, 2)/vint(2)
    mint(72) = 1
    mint(73) = kfr1
    vint(73) = taur1
    vint(74) = gamr1
  End If
  If (isub==141) Then
    kfr2 = 23
    taur2 = pmas(kfr2, 1)**2/vint(2)
    gamr2 = pmas(kfr2, 1)*pmas(kfr2, 2)/vint(2)
    mint(72) = 2
    mint(74) = kfr2
    vint(75) = taur2
    vint(76) = gamr2
  End If
  vint(63) = 0.
  vint(64) = 0.
  mint(71) = 0
  vint(71) = ckin(3)
  If (mint(82)>=2) vint(71) = 0.
  If (iset(isub)==2 .Or. iset(isub)==4) Then
    Do i = 1, 2
      If (kfpr(isub,i)==0) Then
      Else If (mstp(42)<=0) Then
        vint(62+i) = pmas(kfpr(isub,i), 1)**2
      Else
        vint(62+i) = ulmass(kfpr(isub,i))**2
      End If
    End Do
    If (min(vint(63),vint(64))<ckin(6)**2) mint(71) = 1
    If (mint(71)==1) vint(71) = max(ckin(3), ckin(5))
  End If
  If (iset(isub)==0) Then
    is = int(1.5+rlu(0))
    vint(63) = vint(3)**2
    vint(64) = vint(4)**2
    If (isub==92 .Or. isub==93) vint(62+is) = parp(111)**2
    If (isub==93) vint(65-is) = parp(111)**2
    sh = vint(2)
    sqm1 = vint(3)**2
    sqm2 = vint(4)**2
    sqm3 = vint(63)
    sqm4 = vint(64)
    sqla12 = (sh-sqm1-sqm2)**2 - 4.*sqm1*sqm2
    sqla34 = (sh-sqm3-sqm4)**2 - 4.*sqm3*sqm4
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1-sqm2)*(sqm3-sqm4)/sh - sh
    thter2 = sqrt(max(0.,sqla12))*sqrt(max(0.,sqla34))/sh
    thl = 0.5*(thter1-thter2)
    thu = 0.5*(thter1+thter2)
    thm = min(max(thl,parp(101)), thu)
    jtmax = 0
    If (isub==92 .Or. isub==93) jtmax = isub - 91
    Do jt = 1, jtmax
      mint(13+3*jt-is*(2*jt-3)) = 1
      sqmmin = vint(59+3*jt-is*(2*jt-3))
      sqmi = vint(8-3*jt+is*(2*jt-3))**2
      sqmj = vint(3*jt-1-is*(2*jt-3))**2
      sqmf = vint(68-3*jt+is*(2*jt-3))
      squa = 0.5*sh/sqmi*((1.+(sqmi-sqmj)/sh)*thm+sqmi-sqmf-sqmj**2/sh+(sqmi+sqmj)*sqmf/sh+(sqmi-sqmj)**2/sh**2*sqmf)
      quar = sh/sqmi*(thm*(thm+sh-sqmi-sqmj-sqmf*(1.-(sqmi-sqmj)/sh))+sqmi*sqmj-sqmj*sqmf*(1.+(sqmi-sqmj-sqmf)/sh))
      sqmmax = squa + sqrt(max(0.,squa**2-quar))
      If (abs(quar/squa**2)<1.E-06) sqmmax = 0.5*quar/squa
      sqmmax = min(sqmmax, (vint(1)-sqrt(sqmf))**2)
      vint(59+3*jt-is*(2*jt-3)) = sqmmin*(sqmmax/sqmmin)**rlu(0)
    End Do
    sqm3 = vint(63)
    sqm4 = vint(64)
    sqla34 = (sh-sqm3-sqm4)**2 - 4.*sqm3*sqm4
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1-sqm2)*(sqm3-sqm4)/sh - sh
    thter2 = sqrt(max(0.,sqla12))*sqrt(max(0.,sqla34))/sh
    thl = 0.5*(thter1-thter2)
    thu = 0.5*(thter1+thter2)
    b = vint(121)
    c = vint(122)
    If (isub==92 .Or. isub==93) Then
      b = 0.5*b
      c = 0.5*c
    End If
    thm = min(max(thl,parp(101)), thu)
    expth = 0.
    tharg = b*(thm-thu)
    If (tharg>-20.) expth = exp(tharg)
    150 th = thu + log(expth+(1.-expth)*rlu(0))/b
    th = max(thm, min(thu,th))
    ratlog = min((b+c*(th+thm))*(th-thm), (b+c*(th+thu))*(th-thu))
    If (ratlog<log(rlu(0))) Goto 150
    vint(21) = 1.
    vint(22) = 0.
    vint(23) = min(1., max(-1.,(2.*th-thter1)/thter2))
  Else If (iset(isub)>=1 .And. iset(isub)<=4) Then
    Call pyklim(1)
    If (mint(51)/=0) Goto 100
    rtau = rlu(0)
    mtau = 1
    If (rtau>coef(isub,1)) mtau = 2
    If (rtau>coef(isub,1)+coef(isub,2)) mtau = 3
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)) mtau = 4
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)) mtau = 5
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)+coef(isub,5)) mtau = 6
    Call pykmap(1, mtau, rlu(0))
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      Call pyklim(4)
      If (mint(51)/=0) Goto 100
      rtaup = rlu(0)
      mtaup = 1
      If (rtaup>coef(isub,15)) mtaup = 2
      Call pykmap(4, mtaup, rlu(0))
    End If
    Call pyklim(2)
    If (mint(51)/=0) Goto 100
    ryst = rlu(0)
    myst = 1
    If (ryst>coef(isub,7)) myst = 2
    If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
    Call pykmap(2, myst, rlu(0))
    Call pyklim(3)
    If (mint(51)/=0) Goto 100
    If (iset(isub)==2 .Or. iset(isub)==4) Then
      rcth = rlu(0)
      mcth = 1
      If (rcth>coef(isub,10)) mcth = 2
      If (rcth>coef(isub,10)+coef(isub,11)) mcth = 3
      If (rcth>coef(isub,10)+coef(isub,11)+coef(isub,12)) mcth = 4
      If (rcth>coef(isub,10)+coef(isub,11)+coef(isub,12)+coef(isub,13)) mcth = 5
      Call pykmap(3, mcth, rlu(0))
    End If
  Else If (iset(isub)==5) Then
    Call pymult(3)
    isub = mint(1)
  End If
  vint(24) = paru(2)*rlu(0)
  mint(51) = 0
  If (isub<=90 .Or. isub>100) Call pyklim(0)
  If (mint(51)/=0) Goto 100
  If (mint(82)==1 .And. mstp(141)>=1) Then
    mcut = 0
    If (msub(91)+msub(92)+msub(93)+msub(94)+msub(95)==0) Call pykcut(mcut)
    If (mcut/=0) Goto 100
  End If
  Call pysigh(nchn, sigs)
  If (mint(82)==1 .And. isub<=90 .Or. isub>=96) Then
    xsec(isub, 2) = xsec(isub, 2) + sigs
  Else If (mint(82)==1) Then
    xsec(isub, 2) = xsec(isub, 2) + xsec(isub, 1)
  End If
  If (mint(43)==4 .And. mstp(82)>=3) Then
    vint(153) = sigs
    Call pymult(4)
  End If
  viol = sigs/xsec(isub, 1)
  If (viol<rlu(0)) Goto 100
  If (mstp(123)<=0) Then
    If (viol>1.) Then
      Write (mstu(11), 1000) viol, ngen(0, 3) + 1
      Write (mstu(11), 1100) isub, vint(21), vint(22), vint(23), vint(26)
      Stop
    End If
  Else If (mstp(123)==1) Then
    If (viol>vint(108)) Then
      vint(108) = viol
    End If
  Else If (viol>vint(108)) Then
    vint(108) = viol
    If (viol>1.) Then
      xdif = xsec(isub, 1)*(viol-1.)
      xsec(isub, 1) = xsec(isub, 1) + xdif
      If (msub(isub)==1 .And. (isub<=90 .Or. isub>96)) xsec(0, 1) = xsec(0, 1) + xdif
      vint(108) = 1.
    End If
  End If
  vint(148) = 1.
  If (mint(43)==4 .And. (isub<=90 .Or. isub>=96) .And. mstp(82)>=3) Then
    Call pymult(5)
    If (vint(150)<rlu(0)) Goto 100
  End If
  If (mint(82)==1 .And. msub(95)==1) Then
    If (isub<=90 .Or. isub>=95) ngen(95, 1) = ngen(95, 1) + 1
    If (isub<=90 .Or. isub>=96) ngen(96, 2) = ngen(96, 2) + 1
  End If
  If (isub<=90 .Or. isub>=96) mint(31) = mint(31) + 1
  rsigs = sigs*rlu(0)
  qt2 = vint(48)
  rqqbar = parp(87)*(1.-(qt2/(qt2+(parp(88)*parp(82))**2))**2)
  If (isub/=95 .And. (isub/=96 .Or. mstp(82)<=1 .Or. rlu(0)>rqqbar)) Then
    Do ichn = 1, nchn
      kfl1 = isig(ichn, 1)
      kfl2 = isig(ichn, 2)
      mint(2) = isig(ichn, 3)
      rsigs = rsigs - sigh(ichn)
      If (rsigs<=0.) Goto 210
    End Do
  Else If (isub==96) Then
    Call pyspli(mint(11), 21, kfl1, kfldum)
    Call pyspli(mint(12), 21, kfl2, kfldum)
    mint(1) = 11
    mint(2) = 1
    If (kfl1==kfl2 .And. rlu(0)<0.5) mint(2) = 2
  Else
    kfl1 = 21
    kfl2 = 21
    rsigs = 6.*rlu(0)
    mint(2) = 1
    If (rsigs>1.) mint(2) = 2
    If (rsigs>2.) mint(2) = 3
  End If
  210 If (mint(2)>10) Then
    mint(1) = mint(2)/10
    mint(2) = mod(mint(2), 10)
  End If
  mint(15) = kfl1
  mint(16) = kfl2
  mint(13) = mint(15)
  mint(14) = mint(16)
  vint(141) = vint(41)
  vint(142) = vint(42)
  Return
  1000 Format (1X, 'Error: maximum violated by', 1P, E11.3, 1X, 'in event', 1X, I7, '.'/1X, 'Execution stopped!')
  1100 Format (1X, 'ISUB = ', I3, '; Point of violation:'/1X, 'tau=', 1P, E11.3, ', y* =', E11.3, ', cthe = ', 0P, F11.7, ', tau'' =', 1P, E11.3)
End Subroutine pyrand
