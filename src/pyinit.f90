Subroutine pyinit(frame, beam, target, win)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
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
  Character *(*) frame, beam, target
  Character chfram*8, chbeam*8, chtarg*8, chmo(12)*3, chlh(2)*6
  Data chmo/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/, chlh/'lepton', 'hadron'/
  Write (mstu(11), *) 'In PYINIT: BEAM,TARGET= ', beam, target
  Call lulist(0)
  chfram = frame // ' '
  chbeam = beam // ' '
  chtarg = target // ' '
  Call pyinki(chfram, chbeam, chtarg, win)
  If (msel/=0) Then
    Do i = 1, 200
      msub(i) = 0
    End Do
  End If
  If (mint(43)==1 .And. (msel==1 .Or. msel==2)) Then
    If (mint(11)+mint(12)==0) msub(1) = 1
    If (mint(11)+mint(12)/=0) msub(2) = 1
  Else If (msel==1) Then
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    If (mstp(82)<=1 .And. ckin(3)<parp(81)) msub(95) = 1
    If (mstp(82)>=2 .And. ckin(3)<parp(82)) msub(95) = 1
  Else If (msel==2) Then
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    msub(91) = 1
    msub(92) = 1
    msub(93) = 1
    msub(95) = 1
  Else If (msel>=4 .And. msel<=8) Then
    msub(81) = 1
    msub(82) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    mdme(mdcy(21,2)+msel-1, 1) = 1
  Else If (msel==10) Then
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
  Else If (msel==11) Then
    msub(1) = 1
  Else If (msel==12) Then
    msub(2) = 1
  Else If (msel==13) Then
    msub(15) = 1
    msub(30) = 1
  Else If (msel==14) Then
    msub(16) = 1
    msub(31) = 1
  Else If (msel==15) Then
    msub(19) = 1
    msub(20) = 1
    msub(22) = 1
    msub(23) = 1
    msub(25) = 1
  Else If (msel==16) Then
    msub(3) = 1
    msub(5) = 1
    msub(8) = 1
    msub(102) = 1
  Else If (msel==17) Then
    msub(24) = 1
    msub(26) = 1
  Else If (msel==21) Then
    msub(141) = 1
  Else If (msel==22) Then
    msub(142) = 1
  Else If (msel==23) Then
    msub(143) = 1
  End If
  mint(44) = 0
  Do isub = 1, 200
    If (mint(43)<4 .And. isub>=91 .And. isub<=96 .And. msub(isub)==1) Then
      Write (mstu(11), 1200) isub, chlh(mint(41)), chlh(mint(42))
      Stop
    Else If (msub(isub)==1 .And. iset(isub)==-1) Then
      Write (mstu(11), 1300) isub
      Stop
    Else If (msub(isub)==1 .And. iset(isub)<=-2) Then
      Write (mstu(11), 1400) isub
      Stop
    Else If (msub(isub)==1) Then
      mint(44) = mint(44) + 1
    End If
  End Do
  If (mint(44)==0) Then
    Write (mstu(11), 1500)
    Stop
  End If
  mint(45) = mint(44) - msub(91) - msub(92) - msub(93) - msub(94)
  mstp(1) = min(4, mstp(1))
  mstu(114) = min(mstu(114), 2*mstp(1))
  mstp(54) = min(mstp(54), 2*mstp(1))
  Do i = -20, 20
    vint(180+i) = 0.
    ia = iabs(i)
    If (ia>=1 .And. ia<=2*mstp(1)) Then
      Do j = 1, mstp(1)
        ib = 2*j - 1 + mod(ia, 2)
        ipm = (5-isign(1,i))/2
        idc = j + mdcy(ia, 2) + 2
        If (mdme(idc,1)==1 .Or. mdme(idc,1)==ipm) vint(180+i) = vint(180+i) + vckm((ia+1)/2, (ib+1)/2)
      End Do
    Else If (ia>=11 .And. ia<=10+2*mstp(1)) Then
      vint(180+i) = 1.
    End If
  End Do
  mstu(111) = mstp(2)
  If (mstp(3)>=1) Then
    alam = parp(1)
    If (mstp(51)==1) alam = 0.2
    If (mstp(51)==2) alam = 0.29
    If (mstp(51)==3) alam = 0.2
    If (mstp(51)==4) alam = 0.4
    If (mstp(51)==11) alam = 0.16
    If (mstp(51)==12) alam = 0.26
    If (mstp(51)==13) alam = 0.36
    parp(1) = alam
    parp(61) = alam
    paru(112) = alam
    parj(81) = alam
  End If
  Call pyinre
  Do i = 0, 200
    Do j = 1, 3
      ngen(i, j) = 0
      xsec(i, j) = 0.
    End Do
  End Do
  vint(108) = 0.
  If (mint(43)==4) Call pyxtot
  If (mstp(121)<=0) Call pymaxi
  If (mstp(131)/=0) Call pyovly(1)
  If (mint(43)==4 .And. (mint(45)/=0 .Or. mstp(131)/=0) .And. mstp(82)>=2) Call pymult(1)
  Return
  1200 Format (1X, 'Error: process number ', I3, ' not meaningful for ', A6, '-', A6, ' interactions.'/1X, 'Execution stopped!')
  1300 Format (1X, 'Error: requested subprocess', I4, ' not implemented.'/1X, 'Execution stopped!')
  1400 Format (1X, 'Error: requested subprocess', I4, ' not existing.'/1X, 'Execution stopped!')
  1500 Format (1X, 'Error: no subprocess switched on.'/1X, 'Execution stopped.')
End Subroutine pyinit
