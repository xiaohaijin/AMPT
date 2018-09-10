Subroutine hijini
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
  Common /rndf77/nseed
  Save
  nsg = 0
  ndr = 0
  ipp = 2212
  ipt = 2212
  If (ihnt2(5)/=0) ipp = ihnt2(5)
  If (ihnt2(6)/=0) ipt = ihnt2(6)
  Do i = 1, ihnt2(1)
    pp(i, 1) = 0.0
    pp(i, 2) = 0.0
    pp(i, 3) = sqrt(hint1(1)**2/4.0-hint1(8)**2)
    pp(i, 4) = hint1(1)/2
    pp(i, 5) = hint1(8)
    pp(i, 6) = 0.0
    pp(i, 7) = 0.0
    pp(i, 8) = 0.0
    pp(i, 9) = 0.0
    pp(i, 10) = 0.0
    pp(i, 11) = 0.0
    pp(i, 12) = 0.0
    nfp(i, 3) = ipp
    nfp(i, 4) = ipp
    nfp(i, 5) = 0
    nfp(i, 6) = 0
    nfp(i, 7) = 0
    nfp(i, 8) = 0
    nfp(i, 9) = 0
    nfp(i, 10) = 0
    nfp(i, 11) = 0
    npj(i) = 0
    If (i>abs(ihnt2(2))) nfp(i, 3) = 2112
    If (i>abs(ihnt2(2))) nfp(i, 4) = 2112
    Call attflv(nfp(i,3), idq, idqq)
    nfp(i, 1) = idq
    nfp(i, 2) = idqq
    nfp(i, 15) = -1
    If (abs(idq)>1000 .Or. (abs(idq*idqq)<100 .And. ranart(nseed)<0.5)) nfp(i, 15) = 1
    pp(i, 14) = ulmass(idq)
    pp(i, 15) = ulmass(idqq)
  End Do
  Do i = 1, ihnt2(3)
    pt(i, 1) = 0.0
    pt(i, 2) = 0.0
    pt(i, 3) = -sqrt(hint1(1)**2/4.0-hint1(9)**2)
    pt(i, 4) = hint1(1)/2.0
    pt(i, 5) = hint1(9)
    pt(i, 6) = 0.0
    pt(i, 7) = 0.0
    pt(i, 8) = 0.0
    pt(i, 9) = 0.0
    pt(i, 10) = 0.0
    pt(i, 11) = 0.0
    pt(i, 12) = 0.0
    nft(i, 3) = ipt
    nft(i, 4) = ipt
    nft(i, 5) = 0
    nft(i, 6) = 0
    nft(i, 7) = 0
    nft(i, 8) = 0
    nft(i, 9) = 0
    nft(i, 10) = 0
    nft(i, 11) = 0
    ntj(i) = 0
    If (i>abs(ihnt2(4))) nft(i, 3) = 2112
    If (i>abs(ihnt2(4))) nft(i, 4) = 2112
    Call attflv(nft(i,3), idq, idqq)
    nft(i, 1) = idq
    nft(i, 2) = idqq
    nft(i, 15) = 1
    If (abs(idq)>1000 .Or. (abs(idq*idqq)<100 .And. ranart(nseed)<0.5)) nft(i, 15) = -1
    pt(i, 14) = ulmass(idq)
    pt(i, 15) = ulmass(idqq)
  End Do
  Return
End Subroutine hijini
