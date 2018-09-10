Subroutine hijing(frame, bmin0, bmax0)
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Parameter (maxidl=4001)
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Double Precision ataui, zt1, zt2, zt3
  Double Precision xnprod, etprod, xnfrz, etfrz, dnprod, detpro, dnfrz, detfrz
  Double Precision vxp0, vyp0, vzp0, xstrg0, ystrg0, xstrg, ystrg
  Character frame*8
  Dimension scip(300, 300), rnip(300, 300), sjip(300, 300), jtp(3), ipcol(90000), itcol(90000)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjglbr/nelt, ninthj, nelp, ninp
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
  Common /xydr/rtdr(maxstr, 2)
  Common /rndf77/nseed
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /para1/mul
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  Common /srec1/nsp, nst, nsi
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
  Common /frzout/xnprod(30), etprod(30), xnfrz(30), etfrz(30), dnprod(30), detpro(30), dnfrz(30), detfrz(30)
  Common /anim/nevent, isoft, isflag, izpc
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Common /noprec/nnozpc, itypn(maxidl), gxn(maxidl), gyn(maxidl), gzn(maxidl), ftn(maxidl), pxn(maxidl), pyn(maxidl), pzn(maxidl), een(maxidl), xmn(maxidl)
  Common /lastt/itimeh, bimp
  Common /arevt/iaevt, iarun, miss
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /phihj/iphirp, phirp
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  Save
  bmax = min(bmax0, hipr1(34)+hipr1(35))
  bmin = min(bmin0, bmax)
  If (ihnt2(1)<=1 .And. ihnt2(3)<=1) Then
    bmin = 0.0
    bmax = 2.5*sqrt(hipr1(31)*0.1/hipr1(40))
  End If
  yp(1, 1) = 0.0
  yp(2, 1) = 0.0
  yp(3, 1) = 0.0
  If (ihnt2(1)<=1) Goto 14
  Do kp = 1, ihnt2(1)
    5 r = hirnd(1)
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
    yp(1, kp) = r*sx*cos(phi)
    yp(2, kp) = r*sx*sin(phi)
    yp(3, kp) = r*cx
    If (hipr1(29)==0.0) Goto 10
    Do kp2 = 1, kp - 1
      dnbp1 = (yp(1,kp)-yp(1,kp2))**2
      dnbp2 = (yp(2,kp)-yp(2,kp2))**2
      dnbp3 = (yp(3,kp)-yp(3,kp2))**2
      dnbp = dnbp1 + dnbp2 + dnbp3
      If (dnbp<hipr1(29)*hipr1(29)) Goto 5
    End Do
  10 End Do
  If (ihnt2(1)==2) Then
    rnd1 = max(ranart(nseed), 1.0E-20)
    rnd2 = max(ranart(nseed), 1.0E-20)
    rnd3 = max(ranart(nseed), 1.0E-20)
    r = -(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0+4.38*0.85*log(rnd3)/(4.38+0.85))
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
    r = r/2.
    yp(1, 1) = r*sx*cos(phi)
    yp(2, 1) = r*sx*sin(phi)
    yp(3, 1) = r*cx
    yp(1, 2) = -yp(1, 1)
    yp(2, 2) = -yp(2, 1)
    yp(3, 2) = -yp(3, 1)
  End If
  Do i = 1, ihnt2(1) - 1
    Do j = i + 1, ihnt2(1)
      If (yp(3,i)>yp(3,j)) Goto 12
      y1 = yp(1, i)
      y2 = yp(2, i)
      y3 = yp(3, i)
      yp(1, i) = yp(1, j)
      yp(2, i) = yp(2, j)
      yp(3, i) = yp(3, j)
      yp(1, j) = y1
      yp(2, j) = y2
      yp(3, j) = y3
    12 End Do
  End Do
  14 yt(1, 1) = 0.0
  yt(2, 1) = 0.0
  yt(3, 1) = 0.0
  If (ihnt2(3)<=1) Goto 24
  Do kt = 1, ihnt2(3)
    15 r = hirnd(2)
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
    yt(1, kt) = r*sx*cos(phi)
    yt(2, kt) = r*sx*sin(phi)
    yt(3, kt) = r*cx
    If (hipr1(29)==0.0) Goto 20
    Do kt2 = 1, kt - 1
      dnbt1 = (yt(1,kt)-yt(1,kt2))**2
      dnbt2 = (yt(2,kt)-yt(2,kt2))**2
      dnbt3 = (yt(3,kt)-yt(3,kt2))**2
      dnbt = dnbt1 + dnbt2 + dnbt3
      If (dnbt<hipr1(29)*hipr1(29)) Goto 15
    End Do
  20 End Do
  If (ihnt2(3)==2) Then
    rnd1 = max(ranart(nseed), 1.0E-20)
    rnd2 = max(ranart(nseed), 1.0E-20)
    rnd3 = max(ranart(nseed), 1.0E-20)
    r = -(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0+4.38*0.85*log(rnd3)/(4.38+0.85))
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
    r = r/2.
    yt(1, 1) = r*sx*cos(phi)
    yt(2, 1) = r*sx*sin(phi)
    yt(3, 1) = r*cx
    yt(1, 2) = -yt(1, 1)
    yt(2, 2) = -yt(2, 1)
    yt(3, 2) = -yt(3, 1)
  End If
  Do i = 1, ihnt2(3) - 1
    Do j = i + 1, ihnt2(3)
      If (yt(3,i)<yt(3,j)) Goto 22
      y1 = yt(1, i)
      y2 = yt(2, i)
      y3 = yt(3, i)
      yt(1, i) = yt(1, j)
      yt(2, i) = yt(2, j)
      yt(3, i) = yt(3, j)
      yt(1, j) = y1
      yt(2, j) = y2
      yt(3, j) = y3
    22 End Do
  End Do
  24 miss = -1
  50 miss = miss + 1
  If (miss>maxmiss) Then
    Write (6, *) 'infinite loop happened in  HIJING'
    Stop
  End If
  itest = 0
  natt = 0
  jatt = 0
  eatt = 0.0
  Call hijini
  nlop = 0
  60 nt = 0
  np = 0
  n0 = 0
  n01 = 0
  n10 = 0
  n11 = 0
  nelt = 0
  ninthj = 0
  nelp = 0
  ninp = 0
  nsg = 0
  ncolt = 0
  bb = sqrt(bmin**2+ranart(nseed)*(bmax**2-bmin**2))
  phi = 0.
  If (iphirp==1) phi = 2.0*hipr1(40)*ranart(nseed)
  phirp = phi
  bbx = bb*cos(phi)
  bby = bb*sin(phi)
  hint1(19) = bb
  hint1(20) = phi
  Do jp = 1, ihnt2(1)
    Do jt = 1, ihnt2(3)
      scip(jp, jt) = -1.0
      b2 = (yp(1,jp)+bbx-yt(1,jt))**2 + (yp(2,jp)+bby-yt(2,jt))**2
      r2 = b2*hipr1(40)/hipr1(31)/0.1
      rrb1 = min((yp(1,jp)**2+yp(2,jp)**2)/1.2**2/real(ihnt2(1))**0.6666667, 1.0)
      rrb2 = min((yt(1,jt)**2+yt(2,jt)**2)/1.2**2/real(ihnt2(3))**0.6666667, 1.0)
      aphx1 = hipr1(6)*4.0/3.0*(ihnt2(1)**0.3333333-1.0)*sqrt(1.0-rrb1)
      aphx2 = hipr1(6)*4.0/3.0*(ihnt2(3)**0.3333333-1.0)*sqrt(1.0-rrb2)
      hint1(18) = hint1(14) - aphx1*hint1(15) - aphx2*hint1(16) + aphx1*aphx2*hint1(17)
      If (ihpr2(14)==0 .Or. (ihnt2(1)==1 .And. ihnt2(3)==1)) Then
        gs = 1.0 - exp(-(hipr1(30)+hint1(18))*romg(r2)/hipr1(31))
        rantot = ranart(nseed)
        If (rantot>gs) Goto 70
        Goto 65
      End If
      gstot0 = 2.0*(1.0-exp(-(hipr1(30)+hint1(18))/hipr1(31)/2.0*romg(0.0)))
      r2 = r2/gstot0
      gs = 1.0 - exp(-(hipr1(30)+hint1(18))/hipr1(31)*romg(r2))
      gstot = 2.0*(1.0-sqrt(1.0-gs))
      rantot = ranart(nseed)*gstot0
      If (rantot>gstot) Goto 70
      If (rantot>gs) Then
        Call hijcsc(jp, jt)
        Goto 70
      End If
      65 scip(jp, jt) = r2
      rnip(jp, jt) = rantot
      sjip(jp, jt) = hint1(18)
      ncolt = ncolt + 1
      ipcol(ncolt) = jp
      itcol(ncolt) = jt
    70 End Do
  End Do
  bimp = bb
  Write (6, *) '#impact parameter,nlop,ncolt=', bimp, nlop, ncolt
  If (ncolt==0) Then
    nlop = nlop + 1
    If (nlop<=20 .Or. (ihnt2(1)==1 .And. ihnt2(3)==1)) Goto 60
    Return
  End If
  If (ihpr2(3)/=0) Then
    nhard = 1 + int(ranart(nseed)*(ncolt-1)+0.5)
    nhard = min(nhard, ncolt)
    jphard = ipcol(nhard)
    jthard = itcol(nhard)
  End If
  If (ihpr2(9)==1) Then
    nmini = 1 + int(ranart(nseed)*(ncolt-1)+0.5)
    nmini = min(nmini, ncolt)
    jpmini = ipcol(nmini)
    jtmini = itcol(nmini)
  End If
  Do jp = 1, ihnt2(1)
    Do jt = 1, ihnt2(3)
      If (scip(jp,jt)==-1.0) Goto 200
      nfp(jp, 11) = nfp(jp, 11) + 1
      nft(jt, 11) = nft(jt, 11) + 1
      If (nfp(jp,5)<=1 .And. nft(jt,5)>1) Then
        np = np + 1
        n01 = n01 + 1
      Else If (nfp(jp,5)>1 .And. nft(jt,5)<=1) Then
        nt = nt + 1
        n10 = n10 + 1
      Else If (nfp(jp,5)<=1 .And. nft(jt,5)<=1) Then
        np = np + 1
        nt = nt + 1
        n0 = n0 + 1
      Else If (nfp(jp,5)>1 .And. nft(jt,5)>1) Then
        n11 = n11 + 1
      End If
      jout = 0
      nfp(jp, 10) = 0
      nft(jt, 10) = 0
      If (ihpr2(8)==0 .And. ihpr2(3)==0) Goto 160
      If (nfp(jp,6)<0 .Or. nft(jt,6)<0) Goto 160
      r2 = scip(jp, jt)
      hint1(18) = sjip(jp, jt)
      tt = romg(r2)*hint1(18)/hipr1(31)
      tts = hipr1(30)*romg(r2)/hipr1(31)
      njet = 0
      If (ihpr2(3)/=0 .And. jp==jphard .And. jt==jthard) Then
        Call jetini(jp, jt, 1)
        Call hijhrd(jp, jt, 0, jflg, 0)
        hint1(26) = hint1(47)
        hint1(27) = hint1(48)
        hint1(28) = hint1(49)
        hint1(29) = hint1(50)
        hint1(36) = hint1(67)
        hint1(37) = hint1(68)
        hint1(38) = hint1(69)
        hint1(39) = hint1(70)
        If (abs(hint1(46))>hipr1(11) .And. jflg==2) nfp(jp, 7) = 1
        If (abs(hint1(56))>hipr1(11) .And. jflg==2) nft(jt, 7) = 1
        If (max(abs(hint1(46)),abs(hint1(56)))>hipr1(11) .And. jflg>=3) iasg(nsg, 3) = 1
        ihnt2(9) = ihnt2(14)
        ihnt2(10) = ihnt2(15)
        Do i05 = 1, 5
          hint1(20+i05) = hint1(40+i05)
          hint1(30+i05) = hint1(50+i05)
        End Do
        jout = 1
        If (ihpr2(8)==0) Goto 160
        rrb1 = min((yp(1,jp)**2+yp(2,jp)**2)/1.2**2/real(ihnt2(1))**0.6666667, 1.0)
        rrb2 = min((yt(1,jt)**2+yt(2,jt)**2)/1.2**2/real(ihnt2(3))**0.6666667, 1.0)
        aphx1 = hipr1(6)*4.0/3.0*(ihnt2(1)**0.3333333-1.0)*sqrt(1.0-rrb1)
        aphx2 = hipr1(6)*4.0/3.0*(ihnt2(3)**0.3333333-1.0)*sqrt(1.0-rrb2)
        hint1(65) = hint1(61) - aphx1*hint1(62) - aphx2*hint1(63) + aphx1*aphx2*hint1(64)
        ttrig = romg(r2)*hint1(65)/hipr1(31)
        njet = -1
        xr1 = -alog(exp(-ttrig)+ranart(nseed)*(1.0-exp(-ttrig)))
        106 njet = njet + 1
        xr1 = xr1 - alog(max(ranart(nseed),1.0E-20))
        If (xr1<ttrig) Goto 106
        xr = 0.0
        107 njet = njet + 1
        xr = xr - alog(max(ranart(nseed),1.0E-20))
        If (xr<tt-ttrig) Goto 107
        njet = njet - 1
        Goto 112
      End If
      If (ihpr2(9)==1 .And. jp==jpmini .And. jt==jtmini) Goto 110
      If (ihpr2(8)>0 .And. rnip(jp,jt)<=exp(-tt)*(1.0-exp(-tts))) Goto 160
      110 xr = -alog(exp(-tt)+ranart(nseed)*(1.0-exp(-tt)))
      111 njet = njet + 1
      xr = xr - alog(max(ranart(nseed),1.0E-20))
      If (xr<tt) Goto 111
      112 njet = min(njet, ihpr2(8))
      If (ihpr2(8)<0) njet = abs(ihpr2(8))
      Do ijet = 1, njet
        Call jetini(jp, jt, 0)
        Call hijhrd(jp, jt, jout, jflg, 1)
        If (jflg==0) Goto 160
        If (jflg<0) Then
          If (ihpr2(10)/=0) Write (6, *) 'error occured in HIJHRD'
          Goto 50
        End If
        jout = jout + 1
        If (abs(hint1(46))>hipr1(11) .And. jflg==2) nfp(jp, 7) = 1
        If (abs(hint1(56))>hipr1(11) .And. jflg==2) nft(jt, 7) = 1
        If (max(abs(hint1(46)),abs(hint1(56)))>hipr1(11) .And. jflg>=3) iasg(nsg, 3) = 1
      End Do
      160 Continue
      Call hijsft(jp, jt, jout, ierror)
      If (ierror/=0) Then
        If (ihpr2(10)/=0) Write (6, *) 'error occured in HIJSFT'
        Goto 50
      End If
      jatt = jatt + jout
    200 End Do
  End Do
  Call minijet_out(bb, phirp)
  If (pttrig>0 .And. ntrig==0) Goto 50
  Do jp = 1, ihnt2(1)
    If (nfp(jp,5)>2) Then
      ninp = ninp + 1
    Else If (nfp(jp,5)==2 .Or. nfp(jp,5)==1) Then
      nelp = nelp + 1
    End If
  End Do
  Do jt = 1, ihnt2(3)
    If (nft(jt,5)>2) Then
      ninthj = ninthj + 1
    Else If (nft(jt,5)==2 .Or. nft(jt,5)==1) Then
      nelt = nelt + 1
    End If
  End Do
  If ((ihpr2(8)/=0 .Or. ihpr2(3)/=0) .And. ihpr2(4)>0 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
    Do i = 1, ihnt2(1)
      If (nfp(i,7)==1) Call quench(i, 1)
    End Do
    Do i = 1, ihnt2(3)
      If (nft(i,7)==1) Call quench(i, 2)
    End Do
    Do isg = 1, nsg
      If (iasg(isg,3)==1) Call quench(isg, 3)
    End Do
  End If
  If (isoft==1) Then
    isflag = 1
    nsp = ihnt2(1)
    nst = ihnt2(3)
    nsi = nsg
    istr = 0
    npar = 0
    Do i = 1, ihnt2(1)
      istr = istr + 1
      Do j = 1, npj(i)
        If (kfpj(i,j)==21) Then
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = kfpj(i, j)
          gx0(npar) = dble(yp(1,i)+0.5*bb*cos(phirp))
          gy0(npar) = dble(yp(2,i)+0.5*bb*sin(phirp))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pjpx(i,j))
          py0(npar) = dble(pjpy(i,j))
          pz0(npar) = dble(pjpz(i,j))
          xmass0(npar) = dble(pjpm(i,j))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
        End If
      End Do
    End Do
    Do i = 1, ihnt2(3)
      istr = istr + 1
      Do j = 1, ntj(i)
        If (kftj(i,j)==21) Then
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = kftj(i, j)
          gx0(npar) = dble(yt(1,i)-0.5*bb*cos(phirp))
          gy0(npar) = dble(yt(2,i)-0.5*bb*sin(phirp))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pjtx(i,j))
          py0(npar) = dble(pjty(i,j))
          pz0(npar) = dble(pjtz(i,j))
          xmass0(npar) = dble(pjtm(i,j))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
        End If
      End Do
    End Do
    Do i = 1, nsg
      istr = istr + 1
      Do j = 1, njsg(i)
        If (k2sg(i,j)==21) Then
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = k2sg(i, j)
          gx0(npar) = 0.5D0*dble(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
          gy0(npar) = 0.5D0*dble(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pxsg(i,j))
          py0(npar) = dble(pysg(i,j))
          pz0(npar) = dble(pzsg(i,j))
          xmass0(npar) = dble(pmsg(i,j))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
        End If
      End Do
    End Do
    mul = npar
    Call hjana1
    If (ioscar==3) Write (95, *) iaevt, mul
    Call zpcmn
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    Do i = 1, mul
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If
    End Do
    itest = itest + 1
    Do i = 1, mul
      If (lstrg1(i)<=nsp) Then
        nstrg = lstrg1(i)
        npart = lpart1(i)
        kfpj(nstrg, npart) = ityp5(i)
        pjpx(nstrg, npart) = sngl(px5(i))
        pjpy(nstrg, npart) = sngl(py5(i))
        pjpz(nstrg, npart) = sngl(pz5(i))
        pjpe(nstrg, npart) = sngl(e5(i))
        pjpm(nstrg, npart) = sngl(xmass5(i))
      Else If (lstrg1(i)<=nsp+nst) Then
        nstrg = lstrg1(i) - nsp
        npart = lpart1(i)
        kftj(nstrg, npart) = ityp5(i)
        pjtx(nstrg, npart) = sngl(px5(i))
        pjty(nstrg, npart) = sngl(py5(i))
        pjtz(nstrg, npart) = sngl(pz5(i))
        pjte(nstrg, npart) = sngl(e5(i))
        pjtm(nstrg, npart) = sngl(xmass5(i))
      Else
        nstrg = lstrg1(i) - nsp - nst
        npart = lpart1(i)
        k2sg(nstrg, npart) = ityp5(i)
        pxsg(nstrg, npart) = sngl(px5(i))
        pysg(nstrg, npart) = sngl(py5(i))
        pzsg(nstrg, npart) = sngl(pz5(i))
        pesg(nstrg, npart) = sngl(e5(i))
        pmsg(nstrg, npart) = sngl(xmass5(i))
      End If
    End Do
    Call hjana2
  Else If (isoft==2) Then
    nsp = ihnt2(1)
    nst = ihnt2(3)
    nsi = nsg
    npar = 0
    istr = 0
    mstj(1) = 0
    ihpr2(1) = 0
    isflag = 0
    If (ihpr2(20)/=0) Then
      Do ntp = 1, 2
        Do jjtp = 1, ihnt2(2*ntp-1)
          istr = istr + 1
          Call hijfrg(jjtp, ntp, ierror)
          If (ntp==1) Then
            npj(jjtp) = max0(n-2, 0)
          Else
            ntj(jjtp) = max0(n-2, 0)
          End If
          Do ii = 1, n
            npar = npar + 1
            lstrg0(npar) = istr
            lpart0(npar) = ii
            ityp0(npar) = k(ii, 2)
            gz0(npar) = 0D0
            ft0(npar) = 0D0
            px0(npar) = dble(p(ii,1))
            py0(npar) = dble(p(ii,2))
            pz0(npar) = dble(p(ii,3))
            xmass0(npar) = dble(p(ii,5))
            e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
            If (ntp==1) Then
              gx0(npar) = dble(yp(1,jjtp)+0.5*bb*cos(phirp))
              gy0(npar) = dble(yp(2,jjtp)+0.5*bb*sin(phirp))
              iityp = ityp0(npar)
              nstrg = lstrg0(npar)
              If (iityp==2112 .Or. iityp==2212) Then
              Else If ((iityp==1 .Or. iityp==2) .And. (ii==1 .Or. ii==n)) Then
                pp(nstrg, 6) = sngl(px0(npar))
                pp(nstrg, 7) = sngl(py0(npar))
                pp(nstrg, 14) = sngl(xmass0(npar))
              Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (ii==1 .Or. ii==n)) Then
                pp(nstrg, 8) = sngl(px0(npar))
                pp(nstrg, 9) = sngl(py0(npar))
                pp(nstrg, 15) = sngl(xmass0(npar))
              Else
                npart = lpart0(npar) - 1
                kfpj(nstrg, npart) = ityp0(npar)
                pjpx(nstrg, npart) = sngl(px0(npar))
                pjpy(nstrg, npart) = sngl(py0(npar))
                pjpz(nstrg, npart) = sngl(pz0(npar))
                pjpe(nstrg, npart) = sngl(e0(npar))
                pjpm(nstrg, npart) = sngl(xmass0(npar))
              End If
            Else
              gx0(npar) = dble(yt(1,jjtp)-0.5*bb*cos(phirp))
              gy0(npar) = dble(yt(2,jjtp)-0.5*bb*sin(phirp))
              iityp = ityp0(npar)
              nstrg = lstrg0(npar) - nsp
              If (iityp==2112 .Or. iityp==2212) Then
              Else If ((iityp==1 .Or. iityp==2) .And. (ii==1 .Or. ii==n)) Then
                pt(nstrg, 6) = sngl(px0(npar))
                pt(nstrg, 7) = sngl(py0(npar))
                pt(nstrg, 14) = sngl(xmass0(npar))
              Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (ii==1 .Or. ii==n)) Then
                pt(nstrg, 8) = sngl(px0(npar))
                pt(nstrg, 9) = sngl(py0(npar))
                pt(nstrg, 15) = sngl(xmass0(npar))
              Else
                npart = lpart0(npar) - 1
                kftj(nstrg, npart) = ityp0(npar)
                pjtx(nstrg, npart) = sngl(px0(npar))
                pjty(nstrg, npart) = sngl(py0(npar))
                pjtz(nstrg, npart) = sngl(pz0(npar))
                pjte(nstrg, npart) = sngl(e0(npar))
                pjtm(nstrg, npart) = sngl(xmass0(npar))
              End If
            End If
          End Do
        End Do
      End Do
      Do isg = 1, nsg
        istr = istr + 1
        Call hijfrg(isg, 3, ierror)
        njsg(isg) = n
        Do ii = 1, n
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = ii
          ityp0(npar) = k(ii, 2)
          gx0(npar) = 0.5D0*dble(yp(1,iasg(isg,1))+yt(1,iasg(isg,2)))
          gy0(npar) = 0.5D0*dble(yp(2,iasg(isg,1))+yt(2,iasg(isg,2)))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(p(ii,1))
          py0(npar) = dble(p(ii,2))
          pz0(npar) = dble(p(ii,3))
          xmass0(npar) = dble(p(ii,5))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
        End Do
      End Do
    End If
    mul = npar
    Call hjana1
    If (ioscar==3) Write (95, *) iaevt, mul
    Call zpcmn
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    itest = itest + 1
    Do i = 1, mul
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If
    End Do
    Do nmom = 1, 5
      Do nstrg = 1, nsp
        pp(nstrg, nmom) = 0.
      End Do
      Do nstrg = 1, nst
        pt(nstrg, nmom) = 0.
      End Do
    End Do
    Do i = 1, mul
      iityp = ityp5(i)
      If (lstrg1(i)<=nsp) Then
        nstrg = lstrg1(i)
        If (iityp==2112 .Or. iityp==2212) Then
          pp(nstrg, 1) = sngl(px5(i))
          pp(nstrg, 2) = sngl(py5(i))
          pp(nstrg, 3) = sngl(pz5(i))
          pp(nstrg, 4) = sngl(e5(i))
          pp(nstrg, 5) = sngl(xmass5(i))
        Else If ((iityp==1 .Or. iityp==2) .And. (lpart1(i)==1 .Or. lpart1(i)==(npj(nstrg)+2))) Then
          pp(nstrg, 6) = sngl(px5(i))
          pp(nstrg, 7) = sngl(py5(i))
          pp(nstrg, 14) = sngl(xmass5(i))
          pp(nstrg, 1) = pp(nstrg, 1) + sngl(px5(i))
          pp(nstrg, 2) = pp(nstrg, 2) + sngl(py5(i))
          pp(nstrg, 3) = pp(nstrg, 3) + sngl(pz5(i))
          pp(nstrg, 4) = pp(nstrg, 4) + sngl(e5(i))
          pp(nstrg, 5) = sqrt(pp(nstrg,4)**2-pp(nstrg,1)**2-pp(nstrg,2)**2-pp(nstrg,3)**2)
        Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (lpart1(i)==1 .Or. lpart1(i)==(npj(nstrg)+2))) Then
          pp(nstrg, 8) = sngl(px5(i))
          pp(nstrg, 9) = sngl(py5(i))
          pp(nstrg, 15) = sngl(xmass5(i))
          pp(nstrg, 1) = pp(nstrg, 1) + sngl(px5(i))
          pp(nstrg, 2) = pp(nstrg, 2) + sngl(py5(i))
          pp(nstrg, 3) = pp(nstrg, 3) + sngl(pz5(i))
          pp(nstrg, 4) = pp(nstrg, 4) + sngl(e5(i))
          pp(nstrg, 5) = sqrt(pp(nstrg,4)**2-pp(nstrg,1)**2-pp(nstrg,2)**2-pp(nstrg,3)**2)
        Else
          npart = lpart1(i) - 1
          kfpj(nstrg, npart) = ityp5(i)
          pjpx(nstrg, npart) = sngl(px5(i))
          pjpy(nstrg, npart) = sngl(py5(i))
          pjpz(nstrg, npart) = sngl(pz5(i))
          pjpe(nstrg, npart) = sngl(e5(i))
          pjpm(nstrg, npart) = sngl(xmass5(i))
        End If
      Else If (lstrg1(i)<=nsp+nst) Then
        nstrg = lstrg1(i) - nsp
        If (iityp==2112 .Or. iityp==2212) Then
          pt(nstrg, 1) = sngl(px5(i))
          pt(nstrg, 2) = sngl(py5(i))
          pt(nstrg, 3) = sngl(pz5(i))
          pt(nstrg, 4) = sngl(e5(i))
          pt(nstrg, 5) = sngl(xmass5(i))
        Else If ((iityp==1 .Or. iityp==2) .And. (lpart1(i)==1 .Or. lpart1(i)==(ntj(nstrg)+2))) Then
          pt(nstrg, 6) = sngl(px5(i))
          pt(nstrg, 7) = sngl(py5(i))
          pt(nstrg, 14) = sngl(xmass5(i))
          pt(nstrg, 1) = pt(nstrg, 1) + sngl(px5(i))
          pt(nstrg, 2) = pt(nstrg, 2) + sngl(py5(i))
          pt(nstrg, 3) = pt(nstrg, 3) + sngl(pz5(i))
          pt(nstrg, 4) = pt(nstrg, 4) + sngl(e5(i))
          pt(nstrg, 5) = sqrt(pt(nstrg,4)**2-pt(nstrg,1)**2-pt(nstrg,2)**2-pt(nstrg,3)**2)
        Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (lpart1(i)==1 .Or. lpart1(i)==(ntj(nstrg)+2))) Then
          pt(nstrg, 8) = sngl(px5(i))
          pt(nstrg, 9) = sngl(py5(i))
          pt(nstrg, 15) = sngl(xmass5(i))
          pt(nstrg, 1) = pt(nstrg, 1) + sngl(px5(i))
          pt(nstrg, 2) = pt(nstrg, 2) + sngl(py5(i))
          pt(nstrg, 3) = pt(nstrg, 3) + sngl(pz5(i))
          pt(nstrg, 4) = pt(nstrg, 4) + sngl(e5(i))
          pt(nstrg, 5) = sqrt(pt(nstrg,4)**2-pt(nstrg,1)**2-pt(nstrg,2)**2-pt(nstrg,3)**2)
        Else
          npart = lpart1(i) - 1
          kftj(nstrg, npart) = ityp5(i)
          pjtx(nstrg, npart) = sngl(px5(i))
          pjty(nstrg, npart) = sngl(py5(i))
          pjtz(nstrg, npart) = sngl(pz5(i))
          pjte(nstrg, npart) = sngl(e5(i))
          pjtm(nstrg, npart) = sngl(xmass5(i))
        End If
      Else
        nstrg = lstrg1(i) - nsp - nst
        npart = lpart1(i)
        k2sg(nstrg, npart) = ityp5(i)
        pxsg(nstrg, npart) = sngl(px5(i))
        pysg(nstrg, npart) = sngl(py5(i))
        pzsg(nstrg, npart) = sngl(pz5(i))
        pesg(nstrg, npart) = sngl(e5(i))
        pmsg(nstrg, npart) = sngl(xmass5(i))
      End If
    End Do
    mstj(1) = 1
    ihpr2(1) = 1
    isflag = 1
    hipr1(1) = 0.94
    Call hjana2
  Else If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    isflag = 0
    If (ihpr2(20)/=0) Then
      Do isg = 1, nsg
        Call hijfrg(isg, 3, ierror)
        nsbst = 1
        idstr = 92
        If (ihpr2(21)==0) Then
          Call luedit(2)
        Else
          551 nsbst = nsbst + 1
          If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 551
          idstr = k(nsbst, 2)
          nsbst = nsbst + 1
        End If
        If (frame=='LAB') Then
          Call hboost
        End If
        nsbstr = 0
        Do i = nsbst, n
          If (k(i,2)==idstr) Then
            nsbstr = nsbstr + 1
            Goto 560
          End If
          k(i, 4) = nsbstr
          natt = natt + 1
          katt(natt, 1) = k(i, 2)
          katt(natt, 2) = 20
          katt(natt, 4) = k(i, 1)
          If (k(i,3)==0) Then
            katt(natt, 3) = 0
          Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
            katt(natt, 3) = 0
          Else
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
          End If
          patt(natt, 1) = p(i, 1)
          patt(natt, 2) = p(i, 2)
          patt(natt, 3) = p(i, 3)
          patt(natt, 4) = p(i, 4)
          eatt = eatt + p(i, 4)
          gxar(natt) = 0.5*(yp(1,iasg(isg,1))+yt(1,iasg(isg,2)))
          gyar(natt) = 0.5*(yp(2,iasg(isg,1))+yt(2,iasg(isg,2)))
          gzar(natt) = 0.
          ftar(natt) = 0.
          itypar(natt) = k(i, 2)
          pxar(natt) = p(i, 1)
          pyar(natt) = p(i, 2)
          pzar(natt) = p(i, 3)
          pear(natt) = p(i, 4)
          xmar(natt) = p(i, 5)
          xstrg0(natt) = dble(gxar(natt))
          ystrg0(natt) = dble(gyar(natt))
          istrg0(natt) = isg
        560 End Do
      End Do
      jtp(1) = ihnt2(1)
      jtp(2) = ihnt2(3)
      Do ntp = 1, 2
        Do jjtp = 1, jtp(ntp)
          Call hijfrg(jjtp, ntp, ierror)
          nsbst = 1
          idstr = 92
          If (ihpr2(21)==0) Then
            Call luedit(2)
          Else
            581 nsbst = nsbst + 1
            If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 581
            idstr = k(nsbst, 2)
            nsbst = nsbst + 1
          End If
          If (frame=='LAB') Then
            Call hboost
          End If
          nftp = nfp(jjtp, 5)
          If (ntp==2) nftp = 10 + nft(jjtp, 5)
          nsbstr = 0
          Do i = nsbst, n
            If (k(i,2)==idstr) Then
              nsbstr = nsbstr + 1
              Goto 590
            End If
            k(i, 4) = nsbstr
            natt = natt + 1
            katt(natt, 1) = k(i, 2)
            katt(natt, 2) = nftp
            katt(natt, 4) = k(i, 1)
            If (k(i,3)==0) Then
              katt(natt, 3) = 0
            Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
              katt(natt, 3) = 0
            Else
              katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
            End If
            patt(natt, 1) = p(i, 1)
            patt(natt, 2) = p(i, 2)
            patt(natt, 3) = p(i, 3)
            patt(natt, 4) = p(i, 4)
            eatt = eatt + p(i, 4)
            If (ntp==1) Then
              gxar(natt) = yp(1, jjtp) + 0.5*bb*cos(phirp)
              gyar(natt) = yp(2, jjtp) + 0.5*bb*sin(phirp)
            Else
              gxar(natt) = yt(1, jjtp) - 0.5*bb*cos(phirp)
              gyar(natt) = yt(2, jjtp) - 0.5*bb*sin(phirp)
            End If
            gzar(natt) = 0.
            ftar(natt) = 0.
            itypar(natt) = k(i, 2)
            pxar(natt) = p(i, 1)
            pyar(natt) = p(i, 2)
            pzar(natt) = p(i, 3)
            pear(natt) = p(i, 4)
            xmar(natt) = p(i, 5)
            xstrg0(natt) = dble(gxar(natt))
            ystrg0(natt) = dble(gyar(natt))
            istrg0(natt) = ntp*10000 + jjtp
          590 End Do
        End Do
      End Do
    End If
    If (ndr>=1) Then
      Do i = 1, ndr
        natt = natt + 1
        katt(natt, 1) = kfdr(i)
        katt(natt, 2) = 40
        katt(natt, 3) = 0
        patt(natt, 1) = pdr(i, 1)
        patt(natt, 2) = pdr(i, 2)
        patt(natt, 3) = pdr(i, 3)
        patt(natt, 4) = pdr(i, 4)
        eatt = eatt + pdr(i, 4)
        gxar(natt) = rtdr(i, 1)
        gyar(natt) = rtdr(i, 2)
        gzar(natt) = 0.
        ftar(natt) = 0.
        itypar(natt) = katt(natt, 1)
        pxar(natt) = patt(natt, 1)
        pyar(natt) = patt(natt, 2)
        pzar(natt) = patt(natt, 3)
        pear(natt) = patt(natt, 4)
        xmar(natt) = pdr(i, 5)
      End Do
    End If
    Call embedhighpt
    Call hjana1
    Call htop
    nsp = 0
    nst = 0
    nsg = natt
    nsi = nsg
    If (ioscar==3) Write (95, *) iaevt, mul
    Call zpcmn
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    itest = itest + 1
    Do i = 1, mul
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If
    End Do
    Do i = 1, maxstr
      Do j = 1, 3
        k1sgs(i, j) = 0
        k2sgs(i, j) = 0
        pxsgs(i, j) = 0D0
        pysgs(i, j) = 0D0
        pzsgs(i, j) = 0D0
        pesgs(i, j) = 0D0
        pmsgs(i, j) = 0D0
        gxsgs(i, j) = 0D0
        gysgs(i, j) = 0D0
        gzsgs(i, j) = 0D0
        ftsgs(i, j) = 0D0
      End Do
    End Do
    Do i = 1, mul
      iityp = ityp5(i)
      nstrg = lstrg1(i)
      npart = lpart1(i)
      k2sgs(nstrg, npart) = ityp5(i)
      pxsgs(nstrg, npart) = px5(i)
      pysgs(nstrg, npart) = py5(i)
      pzsgs(nstrg, npart) = pz5(i)
      pmsgs(nstrg, npart) = xmass5(i)
      e5(i) = dsqrt(px5(i)**2+py5(i)**2+pz5(i)**2+xmass5(i)**2)
      pesgs(nstrg, npart) = e5(i)
      gxsgs(nstrg, npart) = gx5(i)
      gysgs(nstrg, npart) = gy5(i)
      gzsgs(nstrg, npart) = gz5(i)
      ftsgs(nstrg, npart) = ft5(i)
    End Do
    Call hjana2
  End If
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    natt = 0
    eatt = 0.
    Call ptoh
    Do i = 1, nnozpc
      natt = natt + 1
      katt(natt, 1) = itypn(i)
      patt(natt, 1) = pxn(i)
      patt(natt, 2) = pyn(i)
      patt(natt, 3) = pzn(i)
      patt(natt, 4) = een(i)
      eatt = eatt + een(i)
      gxar(natt) = gxn(i)
      gyar(natt) = gyn(i)
      gzar(natt) = gzn(i)
      ftar(natt) = ftn(i)
      itypar(natt) = itypn(i)
      pxar(natt) = pxn(i)
      pyar(natt) = pyn(i)
      pzar(natt) = pzn(i)
      pear(natt) = een(i)
      xmar(natt) = xmn(i)
    End Do
    Goto 565
  End If
  If (ihpr2(20)/=0) Then
    Do isg = 1, nsg
      Call hijfrg(isg, 3, ierror)
      If (mstu(24)/=0 .Or. ierror>0) Then
        mstu(24) = 0
        mstu(28) = 0
        If (ihpr2(10)/=0) Then
          Write (6, *) 'error occured ISG, repeat the event'
          Write (6, *) isg
        End If
        Goto 50
      End If
      nsbst = 1
      idstr = 92
      If (ihpr2(21)==0) Then
        Call luedit(2)
      Else
        351 nsbst = nsbst + 1
        If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 351
        idstr = k(nsbst, 2)
        nsbst = nsbst + 1
      End If
      If (frame=='LAB') Then
        Call hboost
      End If
      nsbstr = 0
      Do i = nsbst, n
        If (k(i,2)==idstr) Then
          nsbstr = nsbstr + 1
          Goto 360
        End If
        k(i, 4) = nsbstr
        natt = natt + 1
        katt(natt, 1) = k(i, 2)
        katt(natt, 2) = 20
        katt(natt, 4) = k(i, 1)
        If (k(i,3)==0) Then
          katt(natt, 3) = 0
        Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
          katt(natt, 3) = 0
        Else
          katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
        End If
        patt(natt, 1) = p(i, 1)
        patt(natt, 2) = p(i, 2)
        patt(natt, 3) = p(i, 3)
        patt(natt, 4) = p(i, 4)
        eatt = eatt + p(i, 4)
        lsg = nsp + nst + isg
        gxar(natt) = sngl(zt1(lsg))
        gyar(natt) = sngl(zt2(lsg))
        gzar(natt) = sngl(zt3(lsg))
        ftar(natt) = sngl(ataui(lsg))
        itypar(natt) = k(i, 2)
        pxar(natt) = p(i, 1)
        pyar(natt) = p(i, 2)
        pzar(natt) = p(i, 3)
        pear(natt) = p(i, 4)
        xmar(natt) = p(i, 5)
      360 End Do
    End Do
    jtp(1) = ihnt2(1)
    jtp(2) = ihnt2(3)
    Do ntp = 1, 2
      Do jjtp = 1, jtp(ntp)
        Call hijfrg(jjtp, ntp, ierror)
        If (mstu(24)/=0 .Or. ierror>0) Then
          mstu(24) = 0
          mstu(28) = 0
          If (ihpr2(10)/=0) Then
            Write (6, *) 'error occured P&T, repeat the event'
            Write (6, *) ntp, jjtp
          End If
          Goto 50
        End If
        nsbst = 1
        idstr = 92
        If (ihpr2(21)==0) Then
          Call luedit(2)
        Else
          381 nsbst = nsbst + 1
          If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 381
          idstr = k(nsbst, 2)
          nsbst = nsbst + 1
        End If
        If (frame=='LAB') Then
          Call hboost
        End If
        nftp = nfp(jjtp, 5)
        If (ntp==2) nftp = 10 + nft(jjtp, 5)
        nsbstr = 0
        Do i = nsbst, n
          If (k(i,2)==idstr) Then
            nsbstr = nsbstr + 1
            Goto 390
          End If
          k(i, 4) = nsbstr
          natt = natt + 1
          katt(natt, 1) = k(i, 2)
          katt(natt, 2) = nftp
          katt(natt, 4) = k(i, 1)
          If (k(i,3)==0) Then
            katt(natt, 3) = 0
          Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
            katt(natt, 3) = 0
          Else
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
          End If
          patt(natt, 1) = p(i, 1)
          patt(natt, 2) = p(i, 2)
          patt(natt, 3) = p(i, 3)
          patt(natt, 4) = p(i, 4)
          eatt = eatt + p(i, 4)
          If (ntp==1) Then
            lsg = jjtp
          Else
            lsg = jjtp + nsp
          End If
          gxar(natt) = sngl(zt1(lsg))
          gyar(natt) = sngl(zt2(lsg))
          gzar(natt) = sngl(zt3(lsg))
          ftar(natt) = sngl(ataui(lsg))
          itypar(natt) = k(i, 2)
          pxar(natt) = p(i, 1)
          pyar(natt) = p(i, 2)
          pzar(natt) = p(i, 3)
          pear(natt) = p(i, 4)
          xmar(natt) = p(i, 5)
        390 End Do
      End Do
    End Do
  End If
  Do i = 1, ndr
    natt = natt + 1
    katt(natt, 1) = kfdr(i)
    katt(natt, 2) = 40
    katt(natt, 3) = 0
    patt(natt, 1) = pdr(i, 1)
    patt(natt, 2) = pdr(i, 2)
    patt(natt, 3) = pdr(i, 3)
    patt(natt, 4) = pdr(i, 4)
    eatt = eatt + pdr(i, 4)
    gxar(natt) = rtdr(i, 1)
    gyar(natt) = rtdr(i, 2)
    gzar(natt) = 0.
    ftar(natt) = 0.
    itypar(natt) = katt(natt, 1)
    pxar(natt) = patt(natt, 1)
    pyar(natt) = patt(natt, 2)
    pzar(natt) = patt(natt, 3)
    pear(natt) = patt(natt, 4)
    xmar(natt) = pdr(i, 5)
  End Do
  565 Continue
  dengy = eatt/(ihnt2(1)*hint1(6)+ihnt2(3)*hint1(7)) - 1.0
  If (abs(dengy)>hipr1(43) .And. ihpr2(20)/=0 .And. ihpr2(21)==0) Then
    If (ihpr2(10)/=0) Write (6, *) 'Energy not conserved, repeat the event'
    Write (6, *) 'violated:EATT(GeV),NATT,B(fm)=', eatt, natt, bimp
    Goto 50
  End If
  Write (6, *) 'satisfied:EATT(GeV),NATT,B(fm)=', eatt, natt, bimp
  Write (6, *) ' '
  Write (94, *) iaevt, miss, ihnt2(1), ihnt2(3), bimp
  Do jp = 1, ihnt2(1)
    Write (94, 243) yp(1, jp) + 0.5*bb*cos(phirp), yp(2, jp) + 0.5*bb*sin(phirp), jp, nfp(jp, 5), yp(3, jp), nfp(jp, 3), nfp(jp, 4)
  End Do
  Do jt = 1, ihnt2(3)
    Write (94, 243) yt(1, jt) - 0.5*bb*cos(phirp), yt(2, jt) - 0.5*bb*sin(phirp), -jt, nft(jt, 5), yt(3, jt), nft(jt, 3), nft(jt, 4)
  End Do
  Return
  210 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
  211 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
  395 Format (3I8, F10.4, 4I5)
  243 Format (F10.3, 1X, F10.3, 2(1X,I5), 1X, F10.3, 2(1X,I5))
End Subroutine hijing
