Subroutine pythia
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
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
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  mint(7) = 0
  mint(8) = 0
  novl = 1
  If (mstp(131)/=0) Call pyovly(2)
  If (mstp(131)/=0) novl = mint(81)
  mint(83) = 0
  mint(84) = mstp(126)
  mstu(70) = 0
  Do iovl = 1, novl
    If (mint(84)+100>=mstu(4)) Then
      Call luerrm(11, '(PYTHIA:) no more space in LUJETS for overlayed events')
      If (mstu(21)>=1) Goto 200
    End If
    mint(82) = iovl
    100 Continue
    If (iovl==1) ngen(0, 2) = ngen(0, 2) + 1
    mint(31) = 0
    mint(51) = 0
    Call pyrand
    isub = mint(1)
    If (iovl==1) Then
      ngen(isub, 2) = ngen(isub, 2) + 1
      Do j = 1, 200
        msti(j) = 0
        pari(j) = 0.
      End Do
      msti(1) = mint(1)
      msti(2) = mint(2)
      msti(11) = mint(11)
      msti(12) = mint(12)
      msti(15) = mint(15)
      msti(16) = mint(16)
      msti(17) = mint(17)
      msti(18) = mint(18)
      pari(11) = vint(1)
      pari(12) = vint(2)
      If (isub/=95) Then
        Do j = 13, 22
          pari(j) = vint(30+j)
        End Do
        pari(33) = vint(41)
        pari(34) = vint(42)
        pari(35) = pari(33) - pari(34)
        pari(36) = vint(21)
        pari(37) = vint(22)
        pari(38) = vint(26)
        pari(41) = vint(23)
      End If
    End If
    If (mstp(111)==-1) Goto 160
    If (isub<=90 .Or. isub>=95) Then
      Call pyscat
      If (mint(51)==1) Goto 100
      ipu1 = mint(84) + 1
      ipu2 = mint(84) + 2
      If (mstp(61)>=1 .And. mint(43)/=1 .And. isub/=95) Call pysspa(ipu1, ipu2)
      nsav1 = n
      If (mstp(81)>=1 .And. mint(43)==4 .And. isub/=95) Call pymult(6)
      mint(1) = isub
      nsav2 = n
      Call pyremn(ipu1, ipu2)
      If (mint(51)==1) Goto 100
      nsav3 = n
      ipu3 = mint(84) + 3
      ipu4 = mint(84) + 4
      If (mstp(71)>=1 .And. isub/=95 .And. k(ipu3,1)>0 .And. k(ipu3,1)<=10 .And. k(ipu4,1)>0 .And. k(ipu4,1)<=10) Then
        qmax = sqrt(parp(71)*vint(52))
        If (isub==5) qmax = sqrt(pmas(23,1)**2)
        If (isub==8) qmax = sqrt(pmas(24,1)**2)
        Call lushow(ipu3, ipu4, qmax)
      End If
      If (iovl==1) Then
        pari(65) = 2.*pari(17)
        Do i = mstp(126) + 1, n
          If (k(i,1)<=0 .Or. k(i,1)>10) Goto 130
          pt = sqrt(p(i,1)**2+p(i,2)**2)
          pari(69) = pari(69) + pt
          If (i<=nsav1 .Or. i>nsav3) pari(66) = pari(66) + pt
          If (i>nsav1 .And. i<=nsav2) pari(68) = pari(68) + pt
        130 End Do
        pari(67) = pari(68)
        pari(71) = vint(151)
        pari(72) = vint(152)
        pari(73) = vint(151)
        pari(74) = vint(152)
      End If
      If (mstp(41)>=1 .And. isub/=95) Call pyresd
    Else
      Call pydiff
      If (iovl==1) Then
        pari(65) = 2.*pari(17)
        pari(66) = pari(65)
        pari(69) = pari(65)
      End If
    End If
    If (mstp(113)>=1) Then
      Do i = mint(83) + 1, n
        If (k(i,1)>0 .And. k(i,1)<=10) p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      End Do
    End If
    mstu(28) = 0
    Call luprep(mint(84)+1)
    If (mstp(112)==1 .And. mstu(28)==3) Goto 100
    If (mstp(125)==0 .Or. mstp(125)==1) Then
      Do i = mint(84) + 1, n
        If (k(i,2)/=94) Goto 150
        k(i+1, 3) = mod(k(i+1,4)/mstu(5), mstu(5))
        k(i+2, 3) = mod(k(i+2,4)/mstu(5), mstu(5))
      150 End Do
      Call luedit(12)
      Call luedit(14)
      If (mstp(125)==0) Call luedit(15)
      If (mstp(125)==0) mint(4) = 0
    End If
    If (iovl==1 .And. mstp(125)<=0) Then
      mstu(70) = 1
      mstu(71) = n
    Else If (iovl==1) Then
      mstu(70) = 3
      mstu(71) = 2
      mstu(72) = mint(4)
      mstu(73) = n
    End If
    If (mstp(111)>=1) Call luexec
    If (mstp(125)==0 .Or. mstp(125)==1) Call luedit(14)
    160 If (iovl==1) Then
      If (mstp(111)/=-1) ngen(isub, 3) = ngen(isub, 3) + 1
      ngen(0, 3) = ngen(0, 3) + 1
      xsec(0, 3) = 0.
      Do i = 1, 200
        If (i==96) Then
          xsec(i, 3) = 0.
        Else If (msub(95)==1 .And. (i==11 .Or. i==12 .Or. i==13 .Or. i==28 .Or. i==53 .Or. i==68)) Then
          xsec(i, 3) = xsec(96, 2)*ngen(i, 3)/max(1., float(ngen(96,1))*float(ngen(96,2)))
        Else If (ngen(i,1)==0) Then
          xsec(i, 3) = 0.
        Else If (ngen(i,2)==0) Then
          xsec(i, 3) = xsec(i, 2)*ngen(0, 3)/(float(ngen(i,1))*float(ngen(0,2)))
        Else
          xsec(i, 3) = xsec(i, 2)*ngen(i, 3)/(float(ngen(i,1))*float(ngen(i,2)))
        End If
        xsec(0, 3) = xsec(0, 3) + xsec(i, 3)
      End Do
      If (msub(95)==1) Then
        ngens = ngen(91, 3) + ngen(92, 3) + ngen(93, 3) + ngen(94, 3) + ngen(95, 3)
        xsecs = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94, 3) + xsec(95, 3)
        xmaxs = xsec(95, 1)
        If (msub(91)==1) xmaxs = xmaxs + xsec(91, 1)
        If (msub(92)==1) xmaxs = xmaxs + xsec(92, 1)
        If (msub(93)==1) xmaxs = xmaxs + xsec(93, 1)
        If (msub(94)==1) xmaxs = xmaxs + xsec(94, 1)
        fac = 1.
        If (ngens<ngen(0,3)) fac = (xmaxs-xsecs)/(xsec(0,3)-xsecs)
        xsec(11, 3) = fac*xsec(11, 3)
        xsec(12, 3) = fac*xsec(12, 3)
        xsec(13, 3) = fac*xsec(13, 3)
        xsec(28, 3) = fac*xsec(28, 3)
        xsec(53, 3) = fac*xsec(53, 3)
        xsec(68, 3) = fac*xsec(68, 3)
        xsec(0, 3) = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94, 3) + xsec(95, 1)
      End If
      mint(5) = mint(5) + 1
      msti(3) = mint(3)
      msti(4) = mint(4)
      msti(5) = mint(5)
      msti(6) = mint(6)
      msti(7) = mint(7)
      msti(8) = mint(8)
      msti(13) = mint(13)
      msti(14) = mint(14)
      msti(21) = mint(21)
      msti(22) = mint(22)
      msti(23) = mint(23)
      msti(24) = mint(24)
      msti(25) = mint(25)
      msti(26) = mint(26)
      msti(31) = mint(31)
      pari(1) = xsec(0, 3)
      pari(2) = xsec(0, 3)/mint(5)
      pari(31) = vint(141)
      pari(32) = vint(142)
      If (isub/=95 .And. mint(7)*mint(8)/=0) Then
        pari(42) = 2.*vint(47)/vint(1)
        Do is = 7, 8
          pari(36+is) = p(mint(is), 3)/vint(1)
          pari(38+is) = p(mint(is), 4)/vint(1)
          i = mint(is)
          pr = max(1E-20, p(i,5)**2+p(i,1)**2+p(i,2)**2)
          pari(40+is) = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))
          pr = max(1E-20, p(i,1)**2+p(i,2)**2)
          pari(42+is) = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))
          pari(44+is) = p(i, 3)/sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
          pari(46+is) = ulangl(p(i,3), sqrt(p(i,1)**2+p(i,2)**2))
          pari(48+is) = ulangl(p(i,1), p(i,2))
        End Do
      End If
      pari(61) = vint(148)
      If (iset(isub)==1 .Or. iset(isub)==3) Then
        mstu(161) = mint(21)
        mstu(162) = 0
      Else
        mstu(161) = mint(21)
        mstu(162) = mint(22)
      End If
    End If
    msti(41) = iovl
    If (iovl>=2 .And. iovl<=10) msti(40+iovl) = isub
    If (mstu(70)<10) Then
      mstu(70) = mstu(70) + 1
      mstu(70+mstu(70)) = n
    End If
    mint(83) = n
    mint(84) = n + mstp(126)
  End Do
  If (mstp(131)==1 .And. mstp(133)>=1) Then
    pari(91) = vint(132)
    pari(92) = vint(133)
    pari(93) = vint(134)
    If (mstp(133)==2) pari(93) = pari(93)*xsec(0, 3)/vint(131)
  End If
  200 Call pyfram(mstp(124))
  Return
End Subroutine pythia
