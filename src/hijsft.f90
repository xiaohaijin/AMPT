Subroutine hijsft(jp, jt, jout, ierror)
  Parameter (maxstr=150001)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /rndf77/nseed
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /dpmcm1/jjp, jjt, amp, amt, apx0, atx0, ampn, amtn, amp0, amt0, nfdp, nfdt, wp, wm, sw, xremp, xremt, dpkc1, dpkc2, pp11, pp12, pt11, pt12, ptp2, ptt2
  Common /dpmcm2/ndpm, kdpm(20, 2), pdpm1(20, 5), pdpm2(20, 5)
  Save
  ierror = 0
  jjp = jp
  jjt = jt
  ndpm = 0
  If (jp>ihnt2(1) .Or. jt>ihnt2(3)) Return
  epp = pp(jp, 4) + pp(jp, 3)
  epm = pp(jp, 4) - pp(jp, 3)
  etp = pt(jt, 4) + pt(jt, 3)
  etm = pt(jt, 4) - pt(jt, 3)
  wp = epp + etp
  wm = epm + etm
  sw = wp*wm
  If (wp<0.0 .Or. wm<0.0) Goto 1000
  If (jout==0) Then
    If (epp<0.0) Goto 1000
    If (epm<0.0) Goto 1000
    If (etp<0.0) Goto 1000
    If (etm<0.0) Goto 1000
    If (epp/(epm+0.01)<=etp/(etm+0.01)) Return
  End If
  ihnt2(11) = jp
  ihnt2(12) = jt
  miss = 0
  pkc1 = 0.0
  pkc2 = 0.0
  pkc11 = 0.0
  pkc12 = 0.0
  pkc21 = 0.0
  pkc22 = 0.0
  dpkc11 = 0.0
  dpkc12 = 0.0
  dpkc21 = 0.0
  dpkc22 = 0.0
  If (nfp(jp,10)==1 .Or. nft(jt,10)==1) Then
    If (nfp(jp,10)==1) Then
      phi1 = ulangl(pp(jp,10), pp(jp,11))
      ppjet = sqrt(pp(jp,10)**2+pp(jp,11)**2)
      pkc1 = ppjet
      pkc11 = pp(jp, 10)
      pkc12 = pp(jp, 11)
    End If
    If (nft(jt,10)==1) Then
      phi2 = ulangl(pt(jt,10), pt(jt,11))
      ptjet = sqrt(pt(jt,10)**2+pt(jt,11)**2)
      pkc2 = ptjet
      pkc21 = pt(jt, 10)
      pkc22 = pt(jt, 11)
    End If
    If (ihpr2(4)>0 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
      If (nfp(jp,10)==0) Then
        phi = -phi2
      Else If (nft(jt,10)==0) Then
        phi = phi1
      Else
        phi = (phi1+phi2-hipr1(40))/2.0
      End If
      bx = hint1(19)*cos(hint1(20))
      by = hint1(19)*sin(hint1(20))
      xp0 = yp(1, jp)
      yp0 = yp(2, jp)
      xt0 = yt(1, jt) + bx
      yt0 = yt(2, jt) + by
      r1 = max(1.2*ihnt2(1)**0.3333333, sqrt(xp0**2+yp0**2))
      r2 = max(1.2*ihnt2(3)**0.3333333, sqrt((xt0-bx)**2+(yt0-by)**2))
      If (abs(cos(phi))<1.0E-5) Then
        dd1 = r1
        dd2 = r1
        dd3 = abs(by+sqrt(r2**2-(xp0-bx)**2)-yp0)
        dd4 = abs(by-sqrt(r2**2-(xp0-bx)**2)-yp0)
        Goto 5
      End If
      bb = 2.0*sin(phi)*(cos(phi)*yp0-sin(phi)*xp0)
      cc = (yp0**2-r1**2)*cos(phi)**2 + xp0*sin(phi)*(xp0*sin(phi)-2.0*yp0*cos(phi))
      dd = bb**2 - 4.0*cc
      If (dd<0.0) Goto 10
      xx1 = (-bb+sqrt(dd))/2.0
      xx2 = (-bb-sqrt(dd))/2.0
      dd1 = abs((xx1-xp0)/cos(phi))
      dd2 = abs((xx2-xp0)/cos(phi))
      bb = 2.0*sin(phi)*(cos(phi)*(yt0-by)-sin(phi)*xt0) - 2.0*bx
      cc = (bx**2+(yt0-by)**2-r2**2)*cos(phi)**2 + xt0*sin(phi)*(xt0*sin(phi)-2.0*cos(phi)*(yt0-by)) - 2.0*bx*sin(phi)*(cos(phi)*(yt0-by)-sin(phi)*xt0)
      dd = bb**2 - 4.0*cc
      If (dd<0.0) Goto 10
      xx1 = (-bb+sqrt(dd))/2.0
      xx2 = (-bb-sqrt(dd))/2.0
      dd3 = abs((xx1-xt0)/cos(phi))
      dd4 = abs((xx2-xt0)/cos(phi))
      5 dd1 = min(dd1, dd3)
      dd2 = min(dd2, dd4)
      If (dd1<hipr1(13)) dd1 = 0.0
      If (dd2<hipr1(13)) dd2 = 0.0
      If (nfp(jp,10)==1 .And. ppjet>hipr1(11)) Then
        dp1 = dd1*hipr1(14)/2.0
        dp1 = min(dp1, ppjet-hipr1(11))
        pkc1 = ppjet - dp1
        dpx1 = cos(phi1)*dp1
        dpy1 = sin(phi1)*dp1
        pkc11 = pp(jp, 10) - dpx1
        pkc12 = pp(jp, 11) - dpy1
        If (dp1>0.0) Then
          cthep = pp(jp, 12)/sqrt(pp(jp,12)**2+ppjet**2)
          dpz1 = dp1*cthep/sqrt(1.0-cthep**2)
          dpe1 = sqrt(dpx1**2+dpy1**2+dpz1**2)
          eppprm = pp(jp, 4) + pp(jp, 3) - dpe1 - dpz1
          epmprm = pp(jp, 4) - pp(jp, 3) - dpe1 + dpz1
          If (eppprm<=0.0 .Or. epmprm<=0.0) Goto 15
          epp = eppprm
          epm = epmprm
          pp(jp, 10) = pkc11
          pp(jp, 11) = pkc12
          npj(jp) = npj(jp) + 1
          kfpj(jp, npj(jp)) = 21
          pjpx(jp, npj(jp)) = dpx1
          pjpy(jp, npj(jp)) = dpy1
          pjpz(jp, npj(jp)) = dpz1
          pjpe(jp, npj(jp)) = dpe1
          pjpm(jp, npj(jp)) = 0.0
          pp(jp, 3) = pp(jp, 3) - dpz1
          pp(jp, 4) = pp(jp, 4) - dpe1
        End If
      End If
      15 If (nft(jt,10)==1 .And. ptjet>hipr1(11)) Then
        dp2 = dd2*hipr1(14)/2.0
        dp2 = min(dp2, ptjet-hipr1(11))
        pkc2 = ptjet - dp2
        dpx2 = cos(phi2)*dp2
        dpy2 = sin(phi2)*dp2
        pkc21 = pt(jt, 10) - dpx2
        pkc22 = pt(jt, 11) - dpy2
        If (dp2>0.0) Then
          cthet = pt(jt, 12)/sqrt(pt(jt,12)**2+ptjet**2)
          dpz2 = dp2*cthet/sqrt(1.0-cthet**2)
          dpe2 = sqrt(dpx2**2+dpy2**2+dpz2**2)
          etpprm = pt(jt, 4) + pt(jt, 3) - dpe2 - dpz2
          etmprm = pt(jt, 4) - pt(jt, 3) - dpe2 + dpz2
          If (etpprm<=0.0 .Or. etmprm<=0.0) Goto 16
          etp = etpprm
          etm = etmprm
          pt(jt, 10) = pkc21
          pt(jt, 11) = pkc22
          ntj(jt) = ntj(jt) + 1
          kftj(jt, ntj(jt)) = 21
          pjtx(jt, ntj(jt)) = dpx2
          pjty(jt, ntj(jt)) = dpy2
          pjtz(jt, ntj(jt)) = dpz2
          pjte(jt, ntj(jt)) = dpe2
          pjtm(jt, ntj(jt)) = 0.0
          pt(jt, 3) = pt(jt, 3) - dpz2
          pt(jt, 4) = pt(jt, 4) - dpe2
        End If
      End If
      16 dpkc11 = -(pp(jp,10)-pkc11)/2.0
      dpkc12 = -(pp(jp,11)-pkc12)/2.0
      dpkc21 = -(pt(jt,10)-pkc21)/2.0
      dpkc22 = -(pt(jt,11)-pkc22)/2.0
      wp = epp + etp
      wm = epm + etm
      sw = wp*wm
    End If
  End If
  10 ptp02 = pp(jp, 1)**2 + pp(jp, 2)**2
  ptt02 = pt(jt, 1)**2 + pt(jt, 2)**2
  amq = max(pp(jp,14)+pp(jp,15), pt(jt,14)+pt(jt,15))
  amx = hipr1(1) + amq
  amp0 = amx
  dpm0 = amx
  nfdp = 0
  If (nfp(jp,5)<=2 .And. nfp(jp,3)/=0) Then
    amp0 = ulmass(nfp(jp,3))
    nfdp = nfp(jp, 3) + 2*nfp(jp, 3)/abs(nfp(jp,3))
    dpm0 = ulmass(nfdp)
    If (dpm0<=0.0) Then
      nfdp = nfdp - 2*nfdp/abs(nfdp)
      dpm0 = ulmass(nfdp)
    End If
  End If
  amt0 = amx
  dtm0 = amx
  nfdt = 0
  If (nft(jt,5)<=2 .And. nft(jt,3)/=0) Then
    amt0 = ulmass(nft(jt,3))
    nfdt = nft(jt, 3) + 2*nft(jt, 3)/abs(nft(jt,3))
    dtm0 = ulmass(nfdt)
    If (dtm0<=0.0) Then
      nfdt = nfdt - 2*nfdt/abs(nfdt)
      dtm0 = ulmass(nfdt)
    End If
  End If
  ampn = sqrt(amp0**2+ptp02)
  amtn = sqrt(amt0**2+ptt02)
  snn = (ampn+amtn)**2 + 0.001
  If (sw<snn+0.001) Goto 4000
  swptn = 4.0*(max(amp0,amt0)**2+max(ptp02,ptt02))
  swptd = 4.0*(max(dpm0,dtm0)**2+max(ptp02,ptt02))
  swptx = 4.0*(amx**2+max(ptp02,ptt02))
  If (sw<=swptn) Then
    pkcmx = 0.0
  Else If (sw>swptn .And. sw<=swptd .And. npj(jp)==0 .And. ntj(jt)==0) Then
    pkcmx = sqrt(sw/4.0-max(amp0,amt0)**2) - sqrt(max(ptp02,ptt02))
  Else If (sw>swptd .And. sw<=swptx .And. npj(jp)==0 .And. ntj(jt)==0) Then
    pkcmx = sqrt(sw/4.0-max(dpm0,dtm0)**2) - sqrt(max(ptp02,ptt02))
  Else If (sw>swptx) Then
    pkcmx = sqrt(sw/4.0-amx**2) - sqrt(max(ptp02,ptt02))
  End If
  If (nfp(jp,10)==1 .Or. nft(jt,10)==1) Then
    If (pkc1>pkcmx) Then
      pkc1 = pkcmx
      pkc11 = pkc1*cos(phi1)
      pkc12 = pkc1*sin(phi1)
      dpkc11 = -(pp(jp,10)-pkc11)/2.0
      dpkc12 = -(pp(jp,11)-pkc12)/2.0
    End If
    If (pkc2>pkcmx) Then
      pkc2 = pkcmx
      pkc21 = pkc2*cos(phi2)
      pkc22 = pkc2*sin(phi2)
      dpkc21 = -(pt(jt,10)-pkc21)/2.0
      dpkc22 = -(pt(jt,11)-pkc22)/2.0
    End If
    dpkc1 = dpkc11 + dpkc21
    dpkc2 = dpkc12 + dpkc22
    nfp(jp, 10) = -nfp(jp, 10)
    nft(jt, 10) = -nft(jt, 10)
    Goto 40
  End If
  isng = 0
  If (ihpr2(13)/=0 .And. ranart(nseed)<=hidat(4)) isng = 1
  If ((nfp(jp,5)==3 .Or. nft(jt,5)==3) .Or. (npj(jp)/=0 .Or. nfp(jp,10)/=0) .Or. (ntj(jt)/=0 .Or. nft(jt,10)/=0)) isng = 0
  If (ihpr2(5)==0) Then
    pkc = hipr1(2)*sqrt(-alog(1.0-ranart(nseed)*(1.0-exp(-pkcmx**2/hipr1(2)**2))))
    Goto 30
  End If
  xminhi = 0.0
  xmaxhi = pkcmx**2
  pkc = hirnd2(3, xminhi, xmaxhi)
  pkc = sqrt(pkc)
  If (pkc>hipr1(20)) pkc = hipr1(2)*sqrt(-alog(exp(-hipr1(20)**2/hipr1(2)**2)-ranart(nseed)*(exp(-hipr1(20)**2/hipr1(2)**2)-exp(-pkcmx**2/hipr1(2)**2))))
  If (isng==1) pkc = 0.65*sqrt(-alog(1.0-ranart(nseed)*(1.0-exp(-pkcmx**2/0.65**2))))
  30 phi0 = 2.0*hipr1(40)*ranart(nseed)
  pkc11 = pkc*sin(phi0)
  pkc12 = pkc*cos(phi0)
  pkc21 = -pkc11
  pkc22 = -pkc12
  dpkc1 = 0.0
  dpkc2 = 0.0
  40 pp11 = pp(jp, 1) + pkc11 - dpkc1
  pp12 = pp(jp, 2) + pkc12 - dpkc2
  pt11 = pt(jt, 1) + pkc21 - dpkc1
  pt12 = pt(jt, 2) + pkc22 - dpkc2
  ptp2 = pp11**2 + pp12**2
  ptt2 = pt11**2 + pt12**2
  ampn = sqrt(amp0**2+ptp2)
  amtn = sqrt(amt0**2+ptt2)
  snn = (ampn+amtn)**2 + 0.001
  wp = epp + etp
  wm = epm + etm
  sw = wp*wm
  If (sw<snn) Then
    miss = miss + 1
    If (miss<=100) Then
      pkc = 0.0
      Goto 30
    End If
    If (ihpr2(10)/=0) Write (6, *) 'Error occured in Pt kick section of HIJSFT'
    Goto 4000
  End If
  ampd = sqrt(dpm0**2+ptp2)
  amtd = sqrt(dtm0**2+ptt2)
  ampx = sqrt(amx**2+ptp2)
  amtx = sqrt(amx**2+ptt2)
  dpn = ampn**2/sw
  dtn = amtn**2/sw
  dpd = ampd**2/sw
  dtd = amtd**2/sw
  dpx = ampx**2/sw
  dtx = amtx**2/sw
  spntd = (ampn+amtd)**2
  spntx = (ampn+amtx)**2
  spdtn = (ampd+amtn)**2
  spxtn = (ampx+amtn)**2
  spdtx = (ampd+amtx)**2
  spxtd = (ampx+amtd)**2
  sdd = (ampd+amtd)**2
  sxx = (ampx+amtx)**2
  If (sw>sxx+0.001) Then
    If (isng==0) Then
      d1 = dpx
      d2 = dtx
      nfp3 = 0
      nft3 = 0
      Goto 400
    Else
      If ((nfp(jp,5)==3 .And. nft(jt,5)==3) .Or. (npj(jp)/=0 .Or. nfp(jp,10)/=0) .Or. (ntj(jt)/=0 .Or. nft(jt,10)/=0)) Then
        d1 = dpx
        d2 = dtx
        nfp3 = 0
        nft3 = 0
        Goto 400
      End If
      If (ranart(nseed)>0.5 .Or. (nft(jt,5)>2 .Or. ntj(jt)/=0 .Or. nft(jt,10)/=0)) Then
        d1 = dpn
        d2 = dtx
        nfp3 = nfp(jp, 3)
        nft3 = 0
        Goto 220
      Else
        d1 = dpx
        d2 = dtn
        nfp3 = 0
        nft3 = nft(jt, 3)
        Goto 240
      End If
    End If
  Else If (sw>max(spdtx,spxtd)+0.001 .And. sw<=sxx+0.001) Then
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpd
      d2 = dtx
      nfp3 = nfdp
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtd
      nfp3 = 0
      nft3 = nfdt
      Goto 240
    End If
    Goto 4000
  Else If (sw>min(spdtx,spxtd)+0.001 .And. sw<=max(spdtx,spxtd)+0.001) Then
    If (spdtx<=spxtd .And. npj(jp)==0 .And. nfp(jp,5)<=2) Then
      d1 = dpd
      d2 = dtx
      nfp3 = nfdp
      nft3 = 0
      Goto 220
    Else If (spdtx>spxtd .And. ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtd
      nfp3 = 0
      nft3 = nfdt
      Goto 240
    End If
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw>max(spntx,spxtn)+0.001 .And. sw<=min(spdtx,spxtd)+0.001) Then
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw>min(spntx,spxtn)+0.001 .And. sw<=max(spntx,spxtn)+0.001) Then
    If (spntx<=spxtn .And. npj(jp)==0 .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (spntx>spxtn .And. ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw<=min(spntx,spxtn)+0.001 .And. (npj(jp)/=0 .Or. ntj(jt)/=0)) Then
    Goto 4000
  Else If (sw<=min(spntx,spxtn)+0.001 .And. nfp(jp,5)>2 .And. nft(jt,5)>2) Then
    Goto 4000
  Else If (sw>sdd+0.001 .And. sw<=min(spntx,spxtn)+0.001) Then
    d1 = dpd
    d2 = dtd
    nfp3 = nfdp
    nft3 = nfdt
    Goto 100
  Else If (sw>max(spntd,spdtn)+0.001 .And. sw<=sdd+0.001) Then
    If (ranart(nseed)>0.5) Then
      d1 = dpd
      d2 = dtn
      nfp3 = nfdp
      nft3 = nft(jt, 3)
      Goto 100
    Else
      d1 = dpn
      d2 = dtd
      nfp3 = nfp(jp, 3)
      nft3 = nfdt
      Goto 100
    End If
  Else If (sw>min(spntd,spdtn)+0.001 .And. sw<=max(spntd,spdtn)+0.001) Then
    If (spntd>spdtn) Then
      d1 = dpd
      d2 = dtn
      nfp3 = nfdp
      nft3 = nft(jt, 3)
      Goto 100
    Else
      d1 = dpn
      d2 = dtd
      nfp3 = nfp(jp, 3)
      nft3 = nfdt
      Goto 100
    End If
  Else If (sw<=min(spntd,spdtn)+0.001) Then
    d1 = dpn
    d2 = dtn
    nfp3 = nfp(jp, 3)
    nft3 = nft(jt, 3)
    Goto 100
  End If
  Write (6, *) ' Error in HIJSFT: There is no path to here'
  Return
  100 nfp5 = max(2, nfp(jp,5))
  nft5 = max(2, nft(jt,5))
  bb1 = 1.0 + d1 - d2
  bb2 = 1.0 + d2 - d1
  If (bb1**2<4.0*d1 .Or. bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  If (ranart(nseed)<0.5) Then
    x1 = (bb1-sqrt(bb1**2-4.0*d1))/2.0
    x2 = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  Else
    x1 = (bb1+sqrt(bb1**2-4.0*d1))/2.0
    x2 = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  End If
  ihnt2(13) = 2
  Goto 600
  220 nfp5 = max(2, nfp(jp,5))
  nft5 = 3
  If (nfp3==0) nfp5 = 3
  bb2 = 1.0 + d2 - d1
  If (bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  xmax = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  miss4 = 0
  222 x2 = hirnd2(6, xmin, xmax)
  x1 = d1/(1.0-x2)
  If (x2*(1.0-x1)<(d2+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 222
    Goto 5000
  End If
  ihnt2(13) = 2
  Goto 600
  240 nfp5 = 3
  nft5 = max(2, nft(jt,5))
  If (nft3==0) nft5 = 3
  bb1 = 1.0 + d1 - d2
  If (bb1**2<4.0*d1) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin = (bb1-sqrt(bb1**2-4.0*d1))/2.0
  xmax = (bb1+sqrt(bb1**2-4.0*d1))/2.0
  miss4 = 0
  242 x1 = hirnd2(6, xmin, xmax)
  x2 = d2/(1.0-x1)
  If (x1*(1.0-x2)<(d1+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 242
    Goto 5000
  End If
  ihnt2(13) = 2
  Goto 600
  400 nfp5 = 3
  nft5 = 3
  bb1 = 1.0 + d1 - d2
  bb2 = 1.0 + d2 - d1
  If (bb1**2<4.0*d1 .Or. bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin1 = (bb1-sqrt(bb1**2-4.0*d1))/2.0
  xmax1 = (bb1+sqrt(bb1**2-4.0*d1))/2.0
  xmin2 = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  xmax2 = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  miss4 = 0
  410 x1 = hirnd2(4, xmin1, xmax1)
  x2 = hirnd2(4, xmin2, xmax2)
  If (nfp(jp,5)==3 .Or. nft(jt,5)==3) Then
    x1 = hirnd2(6, xmin1, xmax1)
    x2 = hirnd2(6, xmin2, xmax2)
  End If
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000 .Or. abs(nfp(jp,1)*nfp(jp,2))<100) Then
    x1 = hirnd2(5, xmin1, xmax1)
  End If
  If (abs(nft(jt,1)*nft(jt,2))>1000000 .Or. abs(nft(jt,1)*nft(jt,2))<100) Then
    x2 = hirnd2(5, xmin2, xmax2)
  End If
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000) x1 = 1.0 - x1
  xxp = x1*(1.0-x2)
  xxt = x2*(1.0-x1)
  If (xxp<(d1+1.E-4/sw) .Or. xxt<(d2+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 410
    Goto 5000
  End If
  ihnt2(13) = 3
  600 Continue
  If (x1*(1.0-x2)<(ampn**2-1.E-4)/sw .Or. x2*(1.0-x1)<(amtn**2-1.E-4)/sw) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 2000
    pkc = 0.0
    Goto 30
  End If
  epp = (1.0-x2)*wp
  epm = x1*wm
  etp = x2*wp
  etm = (1.0-x1)*wm
  pp(jp, 3) = (epp-epm)/2.0
  pp(jp, 4) = (epp+epm)/2.0
  If (epp*epm-ptp2<0.0) Goto 6000
  pp(jp, 5) = sqrt(epp*epm-ptp2)
  nfp(jp, 3) = nfp3
  nfp(jp, 5) = nfp5
  pt(jt, 3) = (etp-etm)/2.0
  pt(jt, 4) = (etp+etm)/2.0
  If (etp*etm-ptt2<0.0) Goto 6000
  pt(jt, 5) = sqrt(etp*etm-ptt2)
  nft(jt, 3) = nft3
  nft(jt, 5) = nft5
  pp(jp, 1) = pp11 - pkc11
  pp(jp, 2) = pp12 - pkc12
  kcdip = 1
  kcdit = 1
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000 .Or. abs(nfp(jp,1)*nfp(jp,2))<100) Then
    kcdip = 0
  End If
  If (abs(nft(jt,1)*nft(jt,2))>1000000 .Or. abs(nft(jt,1)*nft(jt,2))<100) Then
    kcdit = 0
  End If
  If ((kcdip==0 .And. ranart(nseed)<0.5) .Or. (kcdip/=0 .And. ranart(nseed)<0.5/(1.0+(pkc11**2+pkc12**2)/hipr1(22)**2))) Then
    pp(jp, 6) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 6)
    pp(jp, 7) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 7)
    pp(jp, 8) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 8) + pkc11
    pp(jp, 9) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 9) + pkc12
  Else
    pp(jp, 8) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 8)
    pp(jp, 9) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 9)
    pp(jp, 6) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 6) + pkc11
    pp(jp, 7) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 7) + pkc12
  End If
  pp(jp, 1) = pp(jp, 6) + pp(jp, 8)
  pp(jp, 2) = pp(jp, 7) + pp(jp, 9)
  pt(jt, 1) = pt11 - pkc21
  pt(jt, 2) = pt12 - pkc22
  If ((kcdit==0 .And. ranart(nseed)<0.5) .Or. (kcdit/=0 .And. ranart(nseed)<0.5/(1.0+(pkc21**2+pkc22**2)/hipr1(22)**2))) Then
    pt(jt, 6) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 6)
    pt(jt, 7) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 7)
    pt(jt, 8) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 8) + pkc21
    pt(jt, 9) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 9) + pkc22
  Else
    pt(jt, 8) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 8)
    pt(jt, 9) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 9)
    pt(jt, 6) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 6) + pkc21
    pt(jt, 7) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 7) + pkc22
  End If
  pt(jt, 1) = pt(jt, 6) + pt(jt, 8)
  pt(jt, 2) = pt(jt, 7) + pt(jt, 9)
  If (npj(jp)/=0) nfp(jp, 5) = 3
  If (ntj(jt)/=0) nft(jt, 5) = 3
  If (epp/(epm+0.0001)<etp/(etm+0.0001) .And. abs(nfp(jp,1)*nfp(jp,2))<1000000) Then
    Do jsb = 1, 15
      psb = pp(jp, jsb)
      pp(jp, jsb) = pt(jt, jsb)
      pt(jt, jsb) = psb
      nsb = nfp(jp, jsb)
      nfp(jp, jsb) = nft(jt, jsb)
      nft(jt, jsb) = nsb
    End Do
  End If
  Return
  1000 ierror = 1
  If (ihpr2(10)==0) Return
  Write (6, *) '     Fatal HIJSFT start error,abandon this event'
  Write (6, *) '     PROJ E+,E-,W+', epp, epm, wp
  Write (6, *) '     TARG E+,E-,W-', etp, etm, wm
  Write (6, *) '     W+*W-, (APN+ATN)^2', sw, snn
  Return
  2000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (2)energy partition fail,'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     MP1,MPN', x1*(1.0-x2)*sw, ampn**2
  Write (6, *) '     MT2,MTN', x2*(1.0-x1)*sw, amtn**2
  Return
  3000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (3)something is wrong with the pt kick, '
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     D1=', d1, ' D2=', d2, ' SW=', sw
  Write (6, *) '     HISTORY NFP5=', nfp(jp, 5), ' NFT5=', nft(jt, 5)
  Write (6, *) '     THIS COLLISON NFP5=', nfp5, ' NFT5=', nft5
  Write (6, *) '     # OF JET IN PROJ', npj(jp), ' IN TARG', ntj(jt)
  Return
  4000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (4)unable to choose process, but not harmful'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     PTP=', sqrt(ptp2), ' PTT=', sqrt(ptt2), ' SW=', sw
  Write (6, *) '     AMCUT=', amx, ' JP=', jp, ' JT=', jt
  Write (6, *) '     HISTORY NFP5=', nfp(jp, 5), ' NFT5=', nft(jt, 5)
  Return
  5000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     energy partition failed(5),for limited try'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     NFP5=', nfp5, ' NFT5=', nft5
  Write (6, *) '     D1', d1, ' X1(1-X2)', x1*(1.0-x2)
  Write (6, *) '     D2', d2, ' X2(1-X1)', x2*(1.0-x1)
  Return
  6000 pkc = 0.0
  miss = miss + 1
  If (miss<100) Goto 30
  ierror = 1
  If (ihpr2(10)==0) Return
  Write (6, *) ' ERROR OCCURED, HIJSFT NOT PERFORMED'
  Write (6, *) ' Abort this event'
  Write (6, *) 'MTP,PTP2', epp*epm, ptp2, '  MTT,PTT2', etp*etm, ptt2
  Return
End Subroutine hijsft
