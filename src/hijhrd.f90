Subroutine hijhrd(jp, jt, jout, jflg, iopjt0)
  Parameter (maxstr=150001)
  Dimension ip(100, 2), ipq(50), ipb(50), it(100, 2), itq(50), itb(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
  Common /xydr/rtdr(maxstr, 2)
  Common /rndf77/nseed
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Common /pyint1/mint(400), vint(400)
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
  mxjt = 500
  mxsg = 900
  mxsj = 100
  jflg = 0
  ihnt2(11) = jp
  ihnt2(12) = jt
  iopjet = iopjt0
  If (iopjet==1 .And. (nfp(jp,6)/=0 .Or. nft(jt,6)/=0)) iopjet = 0
  If (jp>ihnt2(1) .Or. jt>ihnt2(3)) Return
  If (nfp(jp,6)<0 .Or. nft(jt,6)<0) Return
  If (jout==0) Then
    epp = pp(jp, 4) + pp(jp, 3)
    epm = pp(jp, 4) - pp(jp, 3)
    etp = pt(jt, 4) + pt(jt, 3)
    etm = pt(jt, 4) - pt(jt, 3)
    If (epp<0.0) Goto 1000
    If (epm<0.0) Goto 1000
    If (etp<0.0) Goto 1000
    If (etm<0.0) Goto 1000
    If (epp/(epm+0.01)<=etp/(etm+0.01)) Return
  End If
  ecut1 = hipr1(1) + hipr1(8) + pp(jp, 14) + pp(jp, 15)
  ecut2 = hipr1(1) + hipr1(8) + pt(jt, 14) + pt(jt, 15)
  If (pp(jp,4)<=ecut1) Then
    nfp(jp, 6) = -abs(nfp(jp,6))
    Return
  End If
  If (pt(jt,4)<=ecut2) Then
    nft(jt, 6) = -abs(nft(jt,6))
    Return
  End If
  miss = 0
  misp = 0
  mist = 0
  If (nfp(jp,10)==0 .And. nft(jt,10)==0) Then
    mint(44) = mint4
    mint(45) = mint5
    xsec(0, 1) = atxs(0)
    xsec(11, 1) = atxs(11)
    xsec(12, 1) = atxs(12)
    xsec(28, 1) = atxs(28)
    Do i = 1, 20
      coef(11, i) = atco(11, i)
      coef(12, i) = atco(12, i)
      coef(28, i) = atco(28, i)
    End Do
  Else
    isub11 = 0
    isub12 = 0
    isub28 = 0
    If (xsec(11,1)/=0) isub11 = 1
    If (xsec(12,1)/=0) isub12 = 1
    If (xsec(28,1)/=0) isub28 = 1
    mint(44) = mint4 - isub11 - isub12 - isub28
    mint(45) = mint5 - isub11 - isub12 - isub28
    xsec(0, 1) = atxs(0) - atxs(11) - atxs(12) - atxs(28)
    xsec(11, 1) = 0.0
    xsec(12, 1) = 0.0
    xsec(28, 1) = 0.0
    Do i = 1, 20
      coef(11, i) = 0.0
      coef(12, i) = 0.0
      coef(28, i) = 0.0
    End Do
  End If
  155 Call pythia
  jj = mint(31)
  If (jj/=1) Goto 155
  If (k(7,2)==-k(8,2)) Then
    qmass2 = (p(7,4)+p(8,4))**2 - (p(7,1)+p(8,1))**2 - (p(7,2)+p(8,2))**2 - (p(7,3)+p(8,3))**2
    qm = ulmass(k(7,2))
    If (qmass2<(2.0*qm+hipr1(1))**2) Goto 155
  End If
  pxp = pp(jp, 1) - p(3, 1)
  pyp = pp(jp, 2) - p(3, 2)
  pzp = pp(jp, 3) - p(3, 3)
  pep = pp(jp, 4) - p(3, 4)
  pxt = pt(jt, 1) - p(4, 1)
  pyt = pt(jt, 2) - p(4, 2)
  pzt = pt(jt, 3) - p(4, 3)
  pet = pt(jt, 4) - p(4, 4)
  If (pep<=ecut1) Then
    misp = misp + 1
    If (misp<50) Goto 155
    nfp(jp, 6) = -abs(nfp(jp,6))
    Return
  End If
  If (pet<=ecut2) Then
    mist = mist + 1
    If (mist<50) Goto 155
    nft(jt, 6) = -abs(nft(jt,6))
    Return
  End If
  wp = pep + pzp + pet + pzt
  wm = pep - pzp + pet - pzt
  If (wp<0.0 .Or. wm<0.0) Then
    miss = miss + 1
    If (pttrig>0) Then
      If (miss<maxmiss) Then
        Write (6, *) 'Failed to generate minijet Pt>', pttrig, 'GeV'
        Goto 155
      End If
    Else
      If (miss<50) Goto 155
    End If
    Return
  End If
  sw = wp*wm
  ampx = sqrt((ecut1-hipr1(8))**2+pxp**2+pyp**2+0.01)
  amtx = sqrt((ecut2-hipr1(8))**2+pxt**2+pyt**2+0.01)
  sxx = (ampx+amtx)**2
  If (sw<sxx .Or. vint(43)<hipr1(1)) Then
    miss = miss + 1
    If (miss>maxmiss) Goto 155
    Return
  End If
  hint1(41) = p(7, 1)
  hint1(42) = p(7, 2)
  hint1(43) = p(7, 3)
  hint1(44) = p(7, 4)
  hint1(45) = p(7, 5)
  hint1(46) = sqrt(p(7,1)**2+p(7,2)**2)
  hint1(51) = p(8, 1)
  hint1(52) = p(8, 2)
  hint1(53) = p(8, 3)
  hint1(54) = p(8, 4)
  hint1(55) = p(8, 5)
  hint1(56) = sqrt(p(8,1)**2+p(8,2)**2)
  ihnt2(14) = k(7, 2)
  ihnt2(15) = k(8, 2)
  pinird = (1.0-exp(-2.0*(vint(47)-hidat(1))))/(1.0+exp(-2.0*(vint(47)-hidat(1))))
  iinird = 0
  If (ranart(nseed)<=pinird) iinird = 1
  If (k(7,2)==-k(8,2)) Goto 190
  If (k(7,2)==21 .And. k(8,2)==21 .And. iopjet==1) Goto 190
  jflg = 2
  jpp = 0
  lpq = 0
  lpb = 0
  jtt = 0
  ltq = 0
  ltb = 0
  is7 = 0
  is8 = 0
  hint1(47) = 0.0
  hint1(48) = 0.0
  hint1(49) = 0.0
  hint1(50) = 0.0
  hint1(67) = 0.0
  hint1(68) = 0.0
  hint1(69) = 0.0
  hint1(70) = 0.0
  Do i = 9, n
    If (k(i,3)==1 .Or. k(i,3)==2 .Or. abs(k(i,2))>30) Goto 180
    If (k(i,3)==7) Then
      hint1(47) = hint1(47) + p(i, 1)
      hint1(48) = hint1(48) + p(i, 2)
      hint1(49) = hint1(49) + p(i, 3)
      hint1(50) = hint1(50) + p(i, 4)
    End If
    If (k(i,3)==8) Then
      hint1(67) = hint1(67) + p(i, 1)
      hint1(68) = hint1(68) + p(i, 2)
      hint1(69) = hint1(69) + p(i, 3)
      hint1(70) = hint1(70) + p(i, 4)
    End If
    If (k(i,2)>21 .And. k(i,2)<=30) Then
      ndr = ndr + 1
      iadr(ndr, 1) = jp
      iadr(ndr, 2) = jt
      kfdr(ndr) = k(i, 2)
      pdr(ndr, 1) = p(i, 1)
      pdr(ndr, 2) = p(i, 2)
      pdr(ndr, 3) = p(i, 3)
      pdr(ndr, 4) = p(i, 4)
      pdr(ndr, 5) = p(i, 5)
      rtdr(ndr, 1) = 0.5*(yp(1,jp)+yt(1,jt))
      rtdr(ndr, 2) = 0.5*(yp(2,jp)+yt(2,jt))
      Goto 180
    End If
    If (k(i,3)==7 .Or. k(i,3)==3) Then
      If (k(i,3)==7 .And. k(i,2)/=21 .And. k(i,2)==k(7,2) .And. is7==0) Then
        pp(jp, 10) = p(i, 1)
        pp(jp, 11) = p(i, 2)
        pp(jp, 12) = p(i, 3)
        pzp = pzp + p(i, 3)
        pep = pep + p(i, 4)
        nfp(jp, 10) = 1
        is7 = 1
        Goto 180
      End If
      If (k(i,3)==3 .And. (k(i,2)/=21 .Or. iinird==0)) Then
        pxp = pxp + p(i, 1)
        pyp = pyp + p(i, 2)
        pzp = pzp + p(i, 3)
        pep = pep + p(i, 4)
        Goto 180
      End If
      jpp = jpp + 1
      ip(jpp, 1) = i
      ip(jpp, 2) = 0
      If (k(i,2)/=21) Then
        If (k(i,2)>0) Then
          lpq = lpq + 1
          ipq(lpq) = jpp
          ip(jpp, 2) = lpq
        Else If (k(i,2)<0) Then
          lpb = lpb + 1
          ipb(lpb) = jpp
          ip(jpp, 2) = -lpb
        End If
      End If
    Else If (k(i,3)==8 .Or. k(i,3)==4) Then
      If (k(i,3)==8 .And. k(i,2)/=21 .And. k(i,2)==k(8,2) .And. is8==0) Then
        pt(jt, 10) = p(i, 1)
        pt(jt, 11) = p(i, 2)
        pt(jt, 12) = p(i, 3)
        pzt = pzt + p(i, 3)
        pet = pet + p(i, 4)
        nft(jt, 10) = 1
        is8 = 1
        Goto 180
      End If
      If (k(i,3)==4 .And. (k(i,2)/=21 .Or. iinird==0)) Then
        pxt = pxt + p(i, 1)
        pyt = pyt + p(i, 2)
        pzt = pzt + p(i, 3)
        pet = pet + p(i, 4)
        Goto 180
      End If
      jtt = jtt + 1
      it(jtt, 1) = i
      it(jtt, 2) = 0
      If (k(i,2)/=21) Then
        If (k(i,2)>0) Then
          ltq = ltq + 1
          itq(ltq) = jtt
          it(jtt, 2) = ltq
        Else If (k(i,2)<0) Then
          ltb = ltb + 1
          itb(ltb) = jtt
          it(jtt, 2) = -ltb
        End If
      End If
    End If
  180 End Do
  If (lpq/=lpb .Or. ltq/=ltb) Then
    miss = miss + 1
    If (miss<=maxmiss) Goto 155
    Write (6, *) ' Q -QBAR NOT MATCHED IN HIJHRD'
    jflg = 0
    Return
  End If
  j = 0
  181 j = j + 1
  If (j>jpp) Goto 182
  If (ip(j,2)==0) Then
    Goto 181
  Else If (ip(j,2)/=0) Then
    lp = abs(ip(j,2))
    ip1 = ip(j, 1)
    ip2 = ip(j, 2)
    ip(j, 1) = ip(ipq(lp), 1)
    ip(j, 2) = ip(ipq(lp), 2)
    ip(ipq(lp), 1) = ip1
    ip(ipq(lp), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipq(lp)
    Else If (ip2<0) Then
      ipb(-ip2) = ipq(lp)
    End If
    ip1 = ip(j+1, 1)
    ip2 = ip(j+1, 2)
    ip(j+1, 1) = ip(ipb(lp), 1)
    ip(j+1, 2) = ip(ipb(lp), 2)
    ip(ipb(lp), 1) = ip1
    ip(ipb(lp), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipb(lp)
    Else If (ip2<0) Then
      ipb(-ip2) = ipb(lp)
    End If
    j = j + 1
    Goto 181
  End If
  182 j = 0
  183 j = j + 1
  If (j>jtt) Goto 184
  If (it(j,2)==0) Then
    Goto 183
  Else If (it(j,2)/=0) Then
    lt = abs(it(j,2))
    it1 = it(j, 1)
    it2 = it(j, 2)
    it(j, 1) = it(itq(lt), 1)
    it(j, 2) = it(itq(lt), 2)
    it(itq(lt), 1) = it1
    it(itq(lt), 2) = it2
    If (it2>0) Then
      itq(it2) = itq(lt)
    Else If (it2<0) Then
      itb(-it2) = itq(lt)
    End If
    it1 = it(j+1, 1)
    it2 = it(j+1, 2)
    it(j+1, 1) = it(itb(lt), 1)
    it(j+1, 2) = it(itb(lt), 2)
    it(itb(lt), 1) = it1
    it(itb(lt), 2) = it2
    If (it2>0) Then
      itq(it2) = itb(lt)
    Else If (it2<0) Then
      itb(-it2) = itb(lt)
    End If
    j = j + 1
    Goto 183
  End If
  184 Continue
  If (npj(jp)+jpp>mxjt .Or. ntj(jt)+jtt>mxjt) Then
    jflg = 0
    Write (6, *) 'number of partons per string exceeds'
    Write (6, *) 'the common block size'
    Return
  End If
  Do j = 1, jpp
    kfpj(jp, npj(jp)+j) = k(ip(j,1), 2)
    pjpx(jp, npj(jp)+j) = p(ip(j,1), 1)
    pjpy(jp, npj(jp)+j) = p(ip(j,1), 2)
    pjpz(jp, npj(jp)+j) = p(ip(j,1), 3)
    pjpe(jp, npj(jp)+j) = p(ip(j,1), 4)
    pjpm(jp, npj(jp)+j) = p(ip(j,1), 5)
  End Do
  npj(jp) = npj(jp) + jpp
  Do j = 1, jtt
    kftj(jt, ntj(jt)+j) = k(it(j,1), 2)
    pjtx(jt, ntj(jt)+j) = p(it(j,1), 1)
    pjty(jt, ntj(jt)+j) = p(it(j,1), 2)
    pjtz(jt, ntj(jt)+j) = p(it(j,1), 3)
    pjte(jt, ntj(jt)+j) = p(it(j,1), 4)
    pjtm(jt, ntj(jt)+j) = p(it(j,1), 5)
  End Do
  ntj(jt) = ntj(jt) + jtt
  Goto 900
  190 jflg = 3
  If (k(7,2)/=21 .And. k(8,2)/=21 .And. k(7,2)*k(8,2)>0) Goto 155
  jpp = 0
  lpq = 0
  lpb = 0
  Do i = 9, n
    If (k(i,3)==1 .Or. k(i,3)==2 .Or. abs(k(i,2))>30) Goto 200
    If (k(i,2)>21 .And. k(i,2)<=30) Then
      ndr = ndr + 1
      iadr(ndr, 1) = jp
      iadr(ndr, 2) = jt
      kfdr(ndr) = k(i, 2)
      pdr(ndr, 1) = p(i, 1)
      pdr(ndr, 2) = p(i, 2)
      pdr(ndr, 3) = p(i, 3)
      pdr(ndr, 4) = p(i, 4)
      pdr(ndr, 5) = p(i, 5)
      rtdr(ndr, 1) = 0.5*(yp(1,jp)+yt(1,jt))
      rtdr(ndr, 2) = 0.5*(yp(2,jp)+yt(2,jt))
      Goto 200
    End If
    If (k(i,3)==3 .And. (k(i,2)/=21 .Or. iinird==0)) Then
      pxp = pxp + p(i, 1)
      pyp = pyp + p(i, 2)
      pzp = pzp + p(i, 3)
      pep = pep + p(i, 4)
      Goto 200
    End If
    If (k(i,3)==4 .And. (k(i,2)/=21 .Or. iinird==0)) Then
      pxt = pxt + p(i, 1)
      pyt = pyt + p(i, 2)
      pzt = pzt + p(i, 3)
      pet = pet + p(i, 4)
      Goto 200
    End If
    jpp = jpp + 1
    ip(jpp, 1) = i
    ip(jpp, 2) = 0
    If (k(i,2)/=21) Then
      If (k(i,2)>0) Then
        lpq = lpq + 1
        ipq(lpq) = jpp
        ip(jpp, 2) = lpq
      Else If (k(i,2)<0) Then
        lpb = lpb + 1
        ipb(lpb) = jpp
        ip(jpp, 2) = -lpb
      End If
    End If
  200 End Do
  If (lpq/=lpb) Then
    miss = miss + 1
    If (miss<=maxmiss) Goto 155
    Write (6, *) lpq, lpb, 'Q-QBAR NOT CONSERVED OR NOT MATCHED'
    jflg = 0
    Return
  End If
  j = 0
  220 j = j + 1
  If (j>jpp) Goto 222
  If (ip(j,2)==0) Goto 220
  lp = abs(ip(j,2))
  ip1 = ip(j, 1)
  ip2 = ip(j, 2)
  ip(j, 1) = ip(ipq(lp), 1)
  ip(j, 2) = ip(ipq(lp), 2)
  ip(ipq(lp), 1) = ip1
  ip(ipq(lp), 2) = ip2
  If (ip2>0) Then
    ipq(ip2) = ipq(lp)
  Else If (ip2<0) Then
    ipb(-ip2) = ipq(lp)
  End If
  ipq(lp) = j
  ip1 = ip(j+1, 1)
  ip2 = ip(j+1, 2)
  ip(j+1, 1) = ip(ipb(lp), 1)
  ip(j+1, 2) = ip(ipb(lp), 2)
  ip(ipb(lp), 1) = ip1
  ip(ipb(lp), 2) = ip2
  If (ip2>0) Then
    ipq(ip2) = ipb(lp)
  Else If (ip2<0) Then
    ipb(-ip2) = ipb(lp)
  End If
  ipb(lp) = j + 1
  j = j + 1
  Goto 220
  222 Continue
  If (lpq>=1) Then
    Do l0 = 2, lpq
      ip1 = ip(2*l0-3, 1)
      ip2 = ip(2*l0-3, 2)
      ip(2*l0-3, 1) = ip(ipq(l0), 1)
      ip(2*l0-3, 2) = ip(ipq(l0), 2)
      ip(ipq(l0), 1) = ip1
      ip(ipq(l0), 2) = ip2
      If (ip2>0) Then
        ipq(ip2) = ipq(l0)
      Else If (ip2<0) Then
        ipb(-ip2) = ipq(l0)
      End If
      ipq(l0) = 2*l0 - 3
      ip1 = ip(2*l0-2, 1)
      ip2 = ip(2*l0-2, 2)
      ip(2*l0-2, 1) = ip(ipb(l0), 1)
      ip(2*l0-2, 2) = ip(ipb(l0), 2)
      ip(ipb(l0), 1) = ip1
      ip(ipb(l0), 2) = ip2
      If (ip2>0) Then
        ipq(ip2) = ipb(l0)
      Else If (ip2<0) Then
        ipb(-ip2) = ipb(l0)
      End If
      ipb(l0) = 2*l0 - 2
    End Do
    ip1 = ip(2*lpq-1, 1)
    ip2 = ip(2*lpq-1, 2)
    ip(2*lpq-1, 1) = ip(ipq(1), 1)
    ip(2*lpq-1, 2) = ip(ipq(1), 2)
    ip(ipq(1), 1) = ip1
    ip(ipq(1), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipq(1)
    Else If (ip2<0) Then
      ipb(-ip2) = ipq(1)
    End If
    ipq(1) = 2*lpq - 1
    ip1 = ip(jpp, 1)
    ip2 = ip(jpp, 2)
    ip(jpp, 1) = ip(ipb(1), 1)
    ip(jpp, 2) = ip(ipb(1), 2)
    ip(ipb(1), 1) = ip1
    ip(ipb(1), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipb(1)
    Else If (ip2<0) Then
      ipb(-ip2) = ipb(1)
    End If
    ipb(1) = jpp
  End If
  If (nsg>=mxsg) Then
    jflg = 0
    Write (6, *) 'number of jets forming single strings exceeds'
    Write (6, *) 'the common block size'
    Return
  End If
  If (jpp>mxsj) Then
    jflg = 0
    Write (6, *) 'number of partons per single jet system'
    Write (6, *) 'exceeds the common block size'
    Return
  End If
  nsg = nsg + 1
  njsg(nsg) = jpp
  iasg(nsg, 1) = jp
  iasg(nsg, 2) = jt
  iasg(nsg, 3) = 0
  Do i = 1, jpp
    k1sg(nsg, i) = 2
    k2sg(nsg, i) = k(ip(i,1), 2)
    If (k2sg(nsg,i)<0) k1sg(nsg, i) = 1
    pxsg(nsg, i) = p(ip(i,1), 1)
    pysg(nsg, i) = p(ip(i,1), 2)
    pzsg(nsg, i) = p(ip(i,1), 3)
    pesg(nsg, i) = p(ip(i,1), 4)
    pmsg(nsg, i) = p(ip(i,1), 5)
  End Do
  k1sg(nsg, 1) = 2
  k1sg(nsg, jpp) = 1
  900 pp(jp, 1) = pxp
  pp(jp, 2) = pyp
  pp(jp, 3) = pzp
  pp(jp, 4) = pep
  pp(jp, 5) = 0.0
  pt(jt, 1) = pxt
  pt(jt, 2) = pyt
  pt(jt, 3) = pzt
  pt(jt, 4) = pet
  pt(jt, 5) = 0.0
  nfp(jp, 6) = nfp(jp, 6) + 1
  nft(jt, 6) = nft(jt, 6) + 1
  Return
  1000 jflg = -1
  If (ihpr2(10)==0) Return
  Write (6, *) 'Fatal HIJHRD error'
  Write (6, *) jp, ' proj E+,E-', epp, epm, ' status', nfp(jp, 5)
  Write (6, *) jt, ' targ E+,E_', etp, etm, ' status', nft(jt, 5)
  Return
End Subroutine hijhrd
