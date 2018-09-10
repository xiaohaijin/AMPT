Subroutine quench(jpjt, ntp)
  Parameter (maxstr=150001)
  Dimension rdp(300), lqp(300), rdt(300), lqt(300)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /rndf77/nseed
  Save
  bb = hint1(19)
  phi = hint1(20)
  bbx = bb*cos(phi)
  bby = bb*sin(phi)
  If (ntp==2) Goto 400
  If (ntp==3) Goto 2000
  If (nfp(jpjt,7)/=1) Return
  jp = jpjt
  Do i = 1, npj(jp)
    ptjet0 = sqrt(pjpx(jp,i)**2+pjpy(jp,i)**2)
    If (ptjet0<=hipr1(11)) Goto 290
    ptot = sqrt(ptjet0*ptjet0+pjpz(jp,i)**2)
    If (ptot<hipr1(8)) Goto 290
    phip = ulangl(pjpx(jp,i), pjpy(jp,i))
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3 .Or. i2==jp) Goto 100
      dx = yp(1, i2) - yp(1, jp)
      dy = yp(2, i2) - yp(2, jp)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phip)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>=hipr1(40)/2.0) Goto 100
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 100
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    100 End Do
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 110
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      110 End Do
    End Do
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3) Goto 120
      dx = yt(1, i2) - yp(1, jp) - bbx
      dy = yt(2, i2) - yp(2, jp) - bby
      phi = ulangl(dx, dy)
      dphi = abs(phi-phip)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 120
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 120
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    120 End Do
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 130
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      130 End Do
    End Do
    mp = 0
    mt = 0
    r0 = 0.0
    nq = 0
    dp = 0.0
    ptot = sqrt(pjpx(jp,i)**2+pjpy(jp,i)**2+pjpz(jp,i)**2)
    v1 = pjpx(jp, i)/ptot
    v2 = pjpy(jp, i)/ptot
    v3 = pjpz(jp, i)/ptot
    200 rn = ranart(nseed)
    210 If (mt>=kt .And. mp>=kp) Goto 290
    If (mt>=kt) Goto 220
    If (mp>=kp) Goto 240
    If (rdp(mp+1)>rdt(mt+1)) Goto 240
    220 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 210
    dp = drr*hipr1(14)
    If (kfpj(jp,i)/=21) dp = 0.5*dp
    If (dp<=0.2) Goto 210
    If (ptot<=0.4) Goto 290
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (kfpj(jp,i)/=21) Then
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pjpm(jp,i)**2+ptot**2) - sqrt(pjpm(jp,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 210
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0
    Goto 260
    240 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 210
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 210
    If (ptot<=0.4) Goto 290
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (kfpj(jp,i)/=21) Then
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pjpm(jp,i)**2+ptot**2) - sqrt(pjpm(jp,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 210
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0
    260 pjpx(jp, i) = (ptot-dp)*v1
    pjpy(jp, i) = (ptot-dp)*v2
    pjpz(jp, i) = (ptot-dp)*v3
    pjpe(jp, i) = pjpe(jp, i) - de
    ptot = ptot - dp
    nq = nq + 1
    Goto 200
  290 End Do
  Return
  400 If (nft(jpjt,7)/=1) Return
  jt = jpjt
  Do i = 1, ntj(jt)
    ptjet0 = sqrt(pjtx(jt,i)**2+pjty(jt,i)**2)
    If (ptjet0<=hipr1(11)) Goto 690
    ptot = sqrt(ptjet0*ptjet0+pjtz(jt,i)**2)
    If (ptot<hipr1(8)) Goto 690
    phit = ulangl(pjtx(jt,i), pjty(jt,i))
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3) Goto 500
      dx = yp(1, i2) + bbx - yt(1, jt)
      dy = yp(2, i2) + bby - yt(2, jt)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phit)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 500
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 500
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    500 End Do
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 510
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      510 End Do
    End Do
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3 .Or. i2==jt) Goto 520
      dx = yt(1, i2) - yt(1, jt)
      dy = yt(2, i2) - yt(2, jt)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phit)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 520
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 520
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    520 End Do
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 530
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      530 End Do
    End Do
    mp = 0
    mt = 0
    nq = 0
    dp = 0.0
    r0 = 0.0
    ptot = sqrt(pjtx(jt,i)**2+pjty(jt,i)**2+pjtz(jt,i)**2)
    v1 = pjtx(jt, i)/ptot
    v2 = pjty(jt, i)/ptot
    v3 = pjtz(jt, i)/ptot
    600 rn = ranart(nseed)
    610 If (mt>=kt .And. mp>=kp) Goto 690
    If (mt>=kt) Goto 620
    If (mp>=kp) Goto 640
    If (rdp(mp+1)>rdt(mt+1)) Goto 640
    620 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 610
    dp = drr*hipr1(14)
    If (kftj(jt,i)/=21) dp = 0.5*dp
    If (dp<=0.2) Goto 610
    If (ptot<=0.4) Goto 690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (kftj(jt,i)/=21) Then
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pjtm(jt,i)**2+ptot**2) - sqrt(pjtm(jt,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 610
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0
    Goto 660
    640 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 610
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 610
    If (ptot<=0.4) Goto 690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (kftj(jt,i)/=21) Then
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pjtm(jt,i)**2+ptot**2) - sqrt(pjtm(jt,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 610
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0
    660 pjtx(jt, i) = (ptot-dp)*v1
    pjty(jt, i) = (ptot-dp)*v2
    pjtz(jt, i) = (ptot-dp)*v3
    pjte(jt, i) = pjte(jt, i) - de
    ptot = ptot - dp
    nq = nq + 1
    Goto 600
  690 End Do
  Return
  2000 isg = jpjt
  If (iasg(isg,3)/=1) Return
  jp = iasg(isg, 1)
  jt = iasg(isg, 2)
  xj = (yp(1,jp)+bbx+yt(1,jt))/2.0
  yj = (yp(2,jp)+bby+yt(2,jt))/2.0
  Do i = 1, njsg(isg)
    ptjet0 = sqrt(pxsg(isg,i)**2+pysg(isg,i)**2)
    If (ptjet0<=hipr1(11) .Or. pesg(isg,i)<hipr1(1)) Goto 2690
    ptot = sqrt(ptjet0*ptjet0+pzsg(isg,i)**2)
    If (ptot<max(hipr1(1),hipr1(8))) Goto 2690
    phiq = ulangl(pxsg(isg,i), pysg(isg,i))
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3 .Or. i2==jp) Goto 2500
      dx = yp(1, i2) + bbx - xj
      dy = yp(2, i2) + bby - yj
      phi = ulangl(dx, dy)
      dphi = abs(phi-phiq)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 2500
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 2500
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    2500 End Do
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 2510
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      2510 End Do
    End Do
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3 .Or. i2==jt) Goto 2520
      dx = yt(1, i2) - xj
      dy = yt(2, i2) - yj
      phi = ulangl(dx, dy)
      dphi = abs(phi-phiq)
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 2520
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 2520
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    2520 End Do
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 2530
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      2530 End Do
    End Do
    mp = 0
    mt = 0
    nq = 0
    dp = 0.0
    r0 = 0.0
    ptot = sqrt(pxsg(isg,i)**2+pysg(isg,i)**2+pzsg(isg,i)**2)
    v1 = pxsg(isg, i)/ptot
    v2 = pysg(isg, i)/ptot
    v3 = pzsg(isg, i)/ptot
    2600 rn = ranart(nseed)
    2610 If (mt>=kt .And. mp>=kp) Goto 2690
    If (mt>=kt) Goto 2620
    If (mp>=kp) Goto 2640
    If (rdp(mp+1)>rdt(mt+1)) Goto 2640
    2620 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 2610
    dp = drr*hipr1(14)/2.0
    If (dp<=0.2) Goto 2610
    If (ptot<=0.4) Goto 2690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (k2sg(isg,i)/=21) Then
      If (ptot<dp+hipr1(1)) Goto 2690
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pmsg(isg,i)**2+ptot**2) - sqrt(pmsg(isg,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 2610
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0
    Goto 2660
    2640 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 2610
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 2610
    If (ptot<=0.4) Goto 2690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
    If (k2sg(isg,i)/=21) Then
      If (ptot<dp+hipr1(1)) Goto 2690
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pmsg(isg,i)**2+ptot**2) - sqrt(pmsg(isg,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 2610
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0
    2660 pxsg(isg, i) = (ptot-dp)*v1
    pysg(isg, i) = (ptot-dp)*v2
    pzsg(isg, i) = (ptot-dp)*v3
    pesg(isg, i) = pesg(isg, i) - de
    ptot = ptot - dp
    nq = nq + 1
    Goto 2600
  2690 End Do
  Return
End Subroutine quench
