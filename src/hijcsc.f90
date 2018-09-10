Subroutine hijcsc(jp, jt)
  Dimension psc1(5), psc2(5)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /rndf77/nseed
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Save
  If (jp==0 .Or. jt==0) Goto 25
  Do i = 1, 5
    psc1(i) = pp(jp, i)
    psc2(i) = pt(jt, i)
  End Do
  Call hijels(psc1, psc2)
  dpp1 = psc1(1) - pp(jp, 1)
  dpp2 = psc1(2) - pp(jp, 2)
  dpt1 = psc2(1) - pt(jt, 1)
  dpt2 = psc2(2) - pt(jt, 2)
  pp(jp, 6) = pp(jp, 6) + dpp1/2.0
  pp(jp, 7) = pp(jp, 7) + dpp2/2.0
  pp(jp, 8) = pp(jp, 8) + dpp1/2.0
  pp(jp, 9) = pp(jp, 9) + dpp2/2.0
  pt(jt, 6) = pt(jt, 6) + dpt1/2.0
  pt(jt, 7) = pt(jt, 7) + dpt2/2.0
  pt(jt, 8) = pt(jt, 8) + dpt1/2.0
  pt(jt, 9) = pt(jt, 9) + dpt2/2.0
  Do i = 1, 4
    pp(jp, i) = psc1(i)
    pt(jt, i) = psc2(i)
  End Do
  nfp(jp, 5) = max(1, nfp(jp,5))
  nft(jt, 5) = max(1, nft(jt,5))
  Return
  25 If (jp==0) Goto 45
  pabs = sqrt(pp(jp,1)**2+pp(jp,2)**2+pp(jp,3)**2)
  bx = pp(jp, 1)/pabs
  by = pp(jp, 2)/pabs
  bz = pp(jp, 3)/pabs
  Do i = 1, ihnt2(1)
    If (i==jp) Goto 40
    dx = yp(1, i) - yp(1, jp)
    dy = yp(2, i) - yp(2, jp)
    dz = yp(3, i) - yp(3, jp)
    dis = dx*bx + dy*by + dz*bz
    If (dis<=0) Goto 40
    bb = dx**2 + dy**2 + dz**2 - dis**2
    r2 = bb*hipr1(40)/hipr1(31)/0.1
    gs = 1.0 - exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(r2))**2
    gs0 = 1.0 - exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(0.0))**2
    If (ranart(nseed)>gs/gs0) Goto 40
    Do k = 1, 5
      psc1(k) = pp(jp, k)
      psc2(k) = pp(i, k)
    End Do
    Call hijels(psc1, psc2)
    dpp1 = psc1(1) - pp(jp, 1)
    dpp2 = psc1(2) - pp(jp, 2)
    dpt1 = psc2(1) - pp(i, 1)
    dpt2 = psc2(2) - pp(i, 2)
    pp(jp, 6) = pp(jp, 6) + dpp1/2.0
    pp(jp, 7) = pp(jp, 7) + dpp2/2.0
    pp(jp, 8) = pp(jp, 8) + dpp1/2.0
    pp(jp, 9) = pp(jp, 9) + dpp2/2.0
    pp(i, 6) = pp(i, 6) + dpt1/2.0
    pp(i, 7) = pp(i, 7) + dpt2/2.0
    pp(i, 8) = pp(i, 8) + dpt1/2.0
    pp(i, 9) = pp(i, 9) + dpt2/2.0
    Do k = 1, 5
      pp(jp, k) = psc1(k)
      pp(i, k) = psc2(k)
    End Do
    nfp(i, 5) = max(1, nfp(i,5))
    Goto 45
  40 End Do
  45 If (jt==0) Goto 80
  pabs = sqrt(pt(jt,1)**2+pt(jt,2)**2+pt(jt,3)**2)
  bx = pt(jt, 1)/pabs
  by = pt(jt, 2)/pabs
  bz = pt(jt, 3)/pabs
  Do i = 1, ihnt2(3)
    If (i==jt) Goto 70
    dx = yt(1, i) - yt(1, jt)
    dy = yt(2, i) - yt(2, jt)
    dz = yt(3, i) - yt(3, jt)
    dis = dx*bx + dy*by + dz*bz
    If (dis<=0) Goto 70
    bb = dx**2 + dy**2 + dz**2 - dis**2
    r2 = bb*hipr1(40)/hipr1(31)/0.1
    gs = (1.0-exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(r2)))**2
    gs0 = (1.0-exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(0.0)))**2
    If (ranart(nseed)>gs/gs0) Goto 70
    Do k = 1, 5
      psc1(k) = pt(jt, k)
      psc2(k) = pt(i, k)
    End Do
    Call hijels(psc1, psc2)
    dpp1 = psc1(1) - pt(jt, 1)
    dpp2 = psc1(2) - pt(jt, 2)
    dpt1 = psc2(1) - pt(i, 1)
    dpt2 = psc2(2) - pt(i, 2)
    pt(jt, 6) = pt(jt, 6) + dpp1/2.0
    pt(jt, 7) = pt(jt, 7) + dpp2/2.0
    pt(jt, 8) = pt(jt, 8) + dpp1/2.0
    pt(jt, 9) = pt(jt, 9) + dpp2/2.0
    pt(i, 6) = pt(i, 6) + dpt1/2.0
    pt(i, 7) = pt(i, 7) + dpt2/2.0
    pt(i, 8) = pt(i, 8) + dpt1/2.0
    pt(i, 9) = pt(i, 9) + dpt2/2.0
    Do k = 1, 5
      pt(jt, k) = psc1(k)
      pt(i, k) = psc2(k)
    End Do
    nft(i, 5) = max(1, nft(i,5))
    Goto 80
  70 End Do
  80 Return
End Subroutine hijcsc
