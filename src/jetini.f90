Subroutine jetini(jp, jt, itrig)
  Character beam*16, targ*16
  Dimension xsec0(8, 0:200), coef0(8, 200, 20), ini(8), mint44(8), mint45(8)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Common /pyint1/mint(400), vint(400)
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save
  Data ini/8*0/, ilast/ -1/
  ihnt2(11) = jp
  ihnt2(12) = jt
  If (ihnt2(5)/=0 .And. ihnt2(6)/=0) Then
    itype = 1
  Else If (ihnt2(5)/=0 .And. ihnt2(6)==0) Then
    itype = 1
    If (nft(jt,4)==2112) itype = 2
  Else If (ihnt2(5)==0 .And. ihnt2(6)/=0) Then
    itype = 1
    If (nfp(jp,4)==2112) itype = 2
  Else
    If (nfp(jp,4)==2212 .And. nft(jt,4)==2212) Then
      itype = 1
    Else If (nfp(jp,4)==2212 .And. nft(jt,4)==2112) Then
      itype = 2
    Else If (nfp(jp,4)==2112 .And. nft(jt,4)==2212) Then
      itype = 3
    Else
      itype = 4
    End If
  End If
  If (itrig/=0) Goto 160
  If (itrig==ilast) Goto 150
  mstp(2) = 2
  mstp(33) = 1
  parp(31) = hipr1(17)
  mstp(51) = 3
  mstp(61) = 1
  mstp(71) = 1
  If (ihpr2(2)==0 .Or. ihpr2(2)==2) mstp(61) = 0
  If (ihpr2(2)==0 .Or. ihpr2(2)==1) mstp(71) = 0
  mstp(81) = 0
  mstp(82) = 1
  mstp(111) = 0
  If (ihpr2(10)==0) mstp(122) = 0
  parp(81) = hipr1(8)
  ckin(5) = hipr1(8)
  ckin(3) = hipr1(8)
  ckin(4) = hipr1(9)
  If (hipr1(9)<=hipr1(8)) ckin(4) = -1.0
  ckin(9) = -10.0
  ckin(10) = 10.0
  msel = 0
  Do isub = 1, 200
    msub(isub) = 0
  End Do
  msub(11) = 1
  msub(12) = 1
  msub(13) = 1
  msub(28) = 1
  msub(53) = 1
  msub(68) = 1
  msub(81) = 1
  msub(82) = 1
  Do j = 1, min(8, mdcy(21,3))
    mdme(mdcy(21,2)+j-1, 1) = 0
  End Do
  isel = 4
  If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
  mdme(mdcy(21,2)+isel-1, 1) = 1
  msub(14) = 1
  msub(18) = 1
  msub(29) = 1
  150 If (ini(itype)/=0) Goto 800
  Goto 400
  160 itype = 4 + itype
  If (itrig==ilast) Goto 260
  parp(81) = abs(hipr1(10)) - 0.25
  ckin(5) = abs(hipr1(10)) - 0.25
  ckin(3) = abs(hipr1(10)) - 0.25
  ckin(4) = abs(hipr1(10)) + 0.25
  If (hipr1(10)<hipr1(8)) ckin(4) = -1.0
  msel = 0
  Do isub = 1, 200
    msub(isub) = 0
  End Do
  If (ihpr2(3)==1) Then
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    msub(81) = 1
    msub(82) = 1
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    isel = 4
    If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
    mdme(mdcy(21,2)+isel-1, 1) = 1
  Else If (ihpr2(3)==2) Then
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
  Else If (ihpr2(3)==3) Then
    ckin(3) = max(0.0, hipr1(10))
    ckin(5) = hipr1(8)
    parp(81) = hipr1(8)
    msub(81) = 1
    msub(82) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    isel = 4
    If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
    mdme(mdcy(21,2)+isel-1, 1) = 1
  End If
  260 If (ini(itype)/=0) Goto 800
  400 ini(itype) = 1
  If (ihpr2(10)==0) mstp(122) = 0
  If (nfp(jp,4)==2212) Then
    beam = 'P'
  Else If (nfp(jp,4)==-2212) Then
    beam = 'P~'
  Else If (nfp(jp,4)==2112) Then
    beam = 'N'
  Else If (nfp(jp,4)==-2112) Then
    beam = 'N~'
  Else If (nfp(jp,4)==211) Then
    beam = 'PI+'
  Else If (nfp(jp,4)==-211) Then
    beam = 'PI-'
  Else If (nfp(jp,4)==321) Then
    beam = 'PI+'
  Else If (nfp(jp,4)==-321) Then
    beam = 'PI-'
  Else
    Write (6, *) 'unavailable beam type', nfp(jp, 4)
  End If
  If (nft(jt,4)==2212) Then
    targ = 'P'
  Else If (nft(jt,4)==-2212) Then
    targ = 'P~'
  Else If (nft(jt,4)==2112) Then
    targ = 'N'
  Else If (nft(jt,4)==-2112) Then
    targ = 'N~'
  Else If (nft(jt,4)==211) Then
    targ = 'PI+'
  Else If (nft(jt,4)==-211) Then
    targ = 'PI-'
  Else If (nft(jt,4)==321) Then
    targ = 'PI+'
  Else If (nft(jt,4)==-321) Then
    targ = 'PI-'
  Else
    Write (6, *) 'unavailable target type', nft(jt, 4)
  End If
  ihnt2(16) = 1
  Call pyinit('CMS', beam, targ, hint1(1))
  mint4 = mint(44)
  mint5 = mint(45)
  mint44(itype) = mint(44)
  mint45(itype) = mint(45)
  atxs(0) = xsec(0, 1)
  xsec0(itype, 0) = xsec(0, 1)
  Do i = 1, 200
    atxs(i) = xsec(i, 1)
    xsec0(itype, i) = xsec(i, 1)
    Do j = 1, 20
      atco(i, j) = coef(i, j)
      coef0(itype, i, j) = coef(i, j)
    End Do
  End Do
  ihnt2(16) = 0
  Return
  800 mint(44) = mint44(itype)
  mint(45) = mint45(itype)
  mint4 = mint(44)
  mint5 = mint(45)
  xsec(0, 1) = xsec0(itype, 0)
  atxs(0) = xsec(0, 1)
  Do i = 1, 200
    xsec(i, 1) = xsec0(itype, i)
    atxs(i) = xsec(i, 1)
    Do j = 1, 20
      coef(i, j) = coef0(itype, i, j)
      atco(i, j) = coef(i, j)
    End Do
  End Do
  ilast = itrig
  mint(11) = nfp(jp, 4)
  mint(12) = nft(jt, 4)
  Return
End Subroutine jetini
