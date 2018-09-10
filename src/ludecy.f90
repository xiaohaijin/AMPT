Subroutine ludecy(ip)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension vdcy(4), kflo(4), kfl1(4), pv(10, 5), rord(10), ue(3), be(3), wtcor(10)
  Common /resdcy/nsav, iksdcy
  Save /resdcy/
  Data wtcor/2., 5., 15., 60., 250., 1500., 1.2E4, 1.2E5, 150., 16./
  pawt(a, b, c) = sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.*a)
  four(i, j) = p(i, 4)*p(j, 4) - p(i, 1)*p(j, 1) - p(i, 2)*p(j, 2) - p(i, 3)*p(j, 3)
  hmeps(ha) = ((1.-hrq-ha)**2+3.*ha*(1.+hrq-ha))*sqrt((1.-hrq-ha)**2-4.*hrq*ha)
  ntry = 0
  nsav = n
  kfa = iabs(k(ip,2))
  kfs = isign(1, k(ip,2))
  kc = lucomp(kfa)
  mstj(92) = 0
  If (k(ip,1)==5) Then
    v(ip, 5) = 0.
  Else If (k(ip,1)/=4) Then
    v(ip, 5) = -pmas(kc, 4)*log(rlu(0))
  End If
  Do j = 1, 4
    vdcy(j) = v(ip, j) + v(ip, 5)*p(ip, j)/p(ip, 5)
  End Do
  mout = 0
  If (mstj(22)==2) Then
    If (pmas(kc,4)>parj(71)) mout = 1
  Else If (mstj(22)==3) Then
    If (vdcy(1)**2+vdcy(2)**2+vdcy(3)**2>parj(72)**2) mout = 1
  Else If (mstj(22)==4) Then
    If (vdcy(1)**2+vdcy(2)**2>parj(73)**2) mout = 1
    If (abs(vdcy(3))>parj(74)) mout = 1
  End If
  If (mout==1 .And. k(ip,1)/=5) Then
    k(ip, 1) = 4
    Return
  End If
  kca = kc
  If (mdcy(kc,2)>0) Then
    mdmdcy = mdme(mdcy(kc,2), 2)
    If (mdmdcy>80 .And. mdmdcy<=90) kca = mdmdcy
  End If
  If (mdcy(kca,2)<=0 .Or. mdcy(kca,3)<=0) Then
    Call luerrm(9, '(LUDECY:) no decay channel defined')
    Return
  End If
  If (mod(kfa/1000,10)==0 .And. (kca==85 .Or. kca==87)) kfs = -kfs
  If (kchg(kc,3)==0) Then
    kfsp = 1
    kfsn = 0
    If (rlu(0)>0.5) kfs = -kfs
  Else If (kfs>0) Then
    kfsp = 1
    kfsn = 0
  Else
    kfsp = 0
    kfsn = 1
  End If
  nope = 0
  brsu = 0.
  Do idl = mdcy(kca, 2), mdcy(kca, 2) + mdcy(kca, 3) - 1
    If (mdme(idl,1)/=1 .And. kfsp*mdme(idl,1)/=2 .And. kfsn*mdme(idl,1)/=3) Goto 120
    If (mdme(idl,2)>100) Goto 120
    nope = nope + 1
    brsu = brsu + brat(idl)
  120 End Do
  If (nope==0) Then
    Call luerrm(2, '(LUDECY:) all decay channels closed by user')
    Return
  End If
  130 rbr = brsu*rlu(0)
  idl = mdcy(kca, 2) - 1
  140 idl = idl + 1
  If (mdme(idl,1)/=1 .And. kfsp*mdme(idl,1)/=2 .And. kfsn*mdme(idl,1)/=3) Then
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1) Goto 140
  Else If (mdme(idl,2)>100) Then
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1) Goto 140
  Else
    idc = idl
    rbr = rbr - brat(idl)
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1 .And. rbr>0.) Goto 140
  End If
  mmat = mdme(idc, 2)
  150 ntry = ntry + 1
  If (ntry>1000) Then
    Call luerrm(14, '(LUDECY:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = n
  np = 0
  nq = 0
  mbst = 0
  If (mmat>=11 .And. mmat/=46 .And. p(ip,4)>20.*p(ip,5)) mbst = 1
  Do j = 1, 4
    pv(1, j) = 0.
    If (mbst==0) pv(1, j) = p(ip, j)
  End Do
  If (mbst==1) pv(1, 4) = p(ip, 5)
  pv(1, 5) = p(ip, 5)
  ps = 0.
  psq = 0.
  mrem = 0
  jtmax = 5
  If (mdme(idc+1,2)==101) jtmax = 10
  Do jt = 1, jtmax
    If (jt<=5) kp = kfdp(idc, jt)
    If (jt>=6) kp = kfdp(idc+1, jt-5)
    If (kp==0) Goto 170
    kpa = iabs(kp)
    kcp = lucomp(kpa)
    If (kchg(kcp,3)==0 .And. kpa/=81 .And. kpa/=82) Then
      kfp = kp
    Else If (kpa/=81 .And. kpa/=82) Then
      kfp = kfs*kp
    Else If (kpa==81 .And. mod(kfa/1000,10)==0) Then
      kfp = -kfs*mod(kfa/10, 10)
    Else If (kpa==81 .And. mod(kfa/100,10)>=mod(kfa/10,10)) Then
      kfp = kfs*(100*mod(kfa/10,100)+3)
    Else If (kpa==81) Then
      kfp = kfs*(1000*mod(kfa/10,10)+100*mod(kfa/100,10)+1)
    Else If (kp==82) Then
      Call lukfdi(-kfs*int(1.+(2.+parj(2))*rlu(0)), 0, kfp, kdump)
      If (kfp==0) Goto 150
      mstj(93) = 1
      If (pv(1,5)<parj(32)+2.*ulmass(kfp)) Goto 150
    Else If (kp==-82) Then
      kfp = -kfp
      If (iabs(kfp)>10) kfp = kfp + isign(10000, kfp)
    End If
    If (kpa==81 .Or. kpa==82) kcp = lucomp(kfp)
    kfpa = iabs(kfp)
    kqp = kchg(kcp, 2)
    If (mmat>=11 .And. mmat<=30 .And. kqp/=0) Then
      nq = nq + 1
      kflo(nq) = kfp
      mstj(93) = 2
      psq = psq + ulmass(kflo(nq))
    Else If (mmat>=42 .And. mmat<=43 .And. np==3 .And. mod(nq,2)==1) Then
      nq = nq - 1
      ps = ps - p(i, 5)
      k(i, 1) = 1
      kfi = k(i, 2)
      Call lukfdi(kfp, kfi, kfldmp, k(i,2))
      If (k(i,2)==0) Goto 150
      mstj(93) = 1
      p(i, 5) = ulmass(k(i,2))
      ps = ps + p(i, 5)
    Else
      i = i + 1
      np = np + 1
      If (mmat/=33 .And. kqp/=0) nq = nq + 1
      If (mmat==33 .And. kqp/=0 .And. kqp/=2) nq = nq + 1
      k(i, 1) = 1 + mod(nq, 2)
      If (mmat==4 .And. jt<=2 .And. kfp==21) k(i, 1) = 2
      If (mmat==4 .And. jt==3) k(i, 1) = 1
      k(i, 2) = kfp
      k(i, 3) = ip
      k(i, 4) = 0
      k(i, 5) = 0
      p(i, 5) = ulmass(kfp)
      If (mmat==45 .And. kfpa==89) p(i, 5) = parj(32)
      ps = ps + p(i, 5)
    End If
  170 End Do
  180 If (mmat>=11 .And. mmat<=30) Then
    psp = ps
    cnde = parj(61)*log(max((pv(1,5)-ps-psq)/parj(62),1.1))
    If (mmat==12) cnde = cnde + parj(63)
    190 ntry = ntry + 1
    If (ntry>1000) Then
      Call luerrm(14, '(LUDECY:) caught in infinite loop')
      If (mstu(21)>=1) Return
    End If
    If (mmat<=20) Then
      gauss = sqrt(-2.*cnde*log(max(1E-10,rlu(0))))*sin(paru(2)*rlu(0))
      nd = int(0.5+0.5*np+0.25*nq+cnde+gauss)
      If (nd<np+nq/2 .Or. nd<2 .Or. nd>10) Goto 190
      If (mmat==13 .And. nd==2) Goto 190
      If (mmat==14 .And. nd<=3) Goto 190
      If (mmat==15 .And. nd<=4) Goto 190
    Else
      nd = mmat - 20
    End If
    Do jt = 1, 4
      kfl1(jt) = kflo(jt)
    End Do
    If (nd==np+nq/2) Goto 220
    Do i = n + np + 1, n + nd - nq/2
      jt = 1 + int((nq-1)*rlu(0))
      Call lukfdi(kfl1(jt), 0, kfl2, k(i,2))
      If (k(i,2)==0) Goto 190
      kfl1(jt) = -kfl2
    End Do
    220 jt = 2
    jt2 = 3
    jt3 = 4
    If (nq==4 .And. rlu(0)<parj(66)) jt = 4
    If (jt==4 .And. isign(1,kfl1(1)*(10-iabs(kfl1(1))))*isign(1,kfl1(jt)*(10-iabs(kfl1(jt))))>0) jt = 3
    If (jt==3) jt2 = 2
    If (jt==4) jt3 = 2
    Call lukfdi(kfl1(1), kfl1(jt), kfldmp, k(n+nd-nq/2+1,2))
    If (k(n+nd-nq/2+1,2)==0) Goto 190
    If (nq==4) Call lukfdi(kfl1(jt2), kfl1(jt3), kfldmp, k(n+nd,2))
    If (nq==4 .And. k(n+nd,2)==0) Goto 190
    ps = psp
    Do i = n + np + 1, n + nd
      k(i, 1) = 1
      k(i, 3) = ip
      k(i, 4) = 0
      k(i, 5) = 0
      p(i, 5) = ulmass(k(i,2))
      ps = ps + p(i, 5)
    End Do
    If (ps+parj(64)>pv(1,5)) Goto 190
  Else If ((mmat==31 .Or. mmat==33 .Or. mmat==44 .Or. mmat==45) .And. np>=3) Then
    ps = ps - p(n+np, 5)
    pqt = (p(n+np,5)+parj(65))/pv(1, 5)
    Do j = 1, 5
      p(n+np, j) = pqt*pv(1, j)
      pv(1, j) = (1.-pqt)*pv(1, j)
    End Do
    If (ps+parj(64)>pv(1,5)) Goto 150
    nd = np - 1
    mrem = 1
  Else If (mmat==46) Then
    mstj(93) = 1
    psmc = ulmass(k(n+1,2))
    mstj(93) = 1
    psmc = psmc + ulmass(k(n+2,2))
    If (max(ps,psmc)+parj(32)>pv(1,5)) Goto 130
    hr1 = (p(n+1,5)/pv(1,5))**2
    hr2 = (p(n+2,5)/pv(1,5))**2
    If ((1.-hr1-hr2)*(2.+hr1+hr2)*sqrt((1.-hr1-hr2)**2-4.*hr1*hr2)<2.*rlu(0)) Goto 130
    nd = np
  Else
    If (np>=2 .And. ps+parj(64)>pv(1,5)) Goto 150
    nd = np
  End If
  If (mmat==45 .And. mstj(25)<=0) Then
    hlq = (parj(32)/pv(1,5))**2
    huq = (1.-(p(n+2,5)+parj(64))/pv(1,5))**2
    hrq = (p(n+2,5)/pv(1,5))**2
    250 hw = hlq + rlu(0)*(huq-hlq)
    If (hmeps(hw)<rlu(0)) Goto 250
    p(n+1, 5) = pv(1, 5)*sqrt(hw)
  Else If (mmat==45) Then
    hqw = (pv(1,5)/pmas(24,1))**2
    hlw = (parj(32)/pmas(24,1))**2
    huw = ((pv(1,5)-p(n+2,5)-parj(64))/pmas(24,1))**2
    hrq = (p(n+2,5)/pv(1,5))**2
    hg = pmas(24, 2)/pmas(24, 1)
    hatl = atan((hlw-1.)/hg)
    hm = min(1., huw-0.001)
    hmv1 = hmeps(hm/hqw)/((hm-1.)**2+hg**2)
    260 hm = hm - hg
    hmv2 = hmeps(hm/hqw)/((hm-1.)**2+hg**2)
    hsav1 = hmeps(hm/hqw)
    hsav2 = 1./((hm-1.)**2+hg**2)
    If (hmv2>hmv1 .And. hm-hg>hlw) Then
      hmv1 = hmv2
      Goto 260
    End If
    hmv = min(2.*hmv1, hmeps(hm/hqw)/hg**2)
    hm1 = 1. - sqrt(1./hmv-hg**2)
    If (hm1>hlw .And. hm1<hm) Then
      hm = hm1
    Else If (hmv2<=hmv1) Then
      hm = max(hlw, hm-min(0.1,1.-hm))
    End If
    hatm = atan((hm-1.)/hg)
    hwt1 = (hatm-hatl)/hg
    hwt2 = hmv*(min(1.,huw)-hm)
    hwt3 = 0.
    If (huw>1.) Then
      hatu = atan((huw-1.)/hg)
      hmp1 = hmeps(1./hqw)
      hwt3 = hmp1*hatu/hg
    End If
    270 hreg = rlu(0)*(hwt1+hwt2+hwt3)
    If (hreg<=hwt1) Then
      hw = 1. + hg*tan(hatl+rlu(0)*(hatm-hatl))
      hacc = hmeps(hw/hqw)
    Else If (hreg<=hwt1+hwt2) Then
      hw = hm + rlu(0)*(min(1.,huw)-hm)
      hacc = hmeps(hw/hqw)/((hw-1.)**2+hg**2)/hmv
    Else
      hw = 1. + hg*tan(rlu(0)*hatu)
      hacc = hmeps(hw/hqw)/hmp1
    End If
    If (hacc<rlu(0)) Goto 270
    p(n+1, 5) = pmas(24, 1)*sqrt(hw)
  End If
  nm = 0
  msgn = 0
  If (mmat==3 .Or. mmat==46) Then
    im = k(ip, 3)
    If (im<0 .Or. im>=ip) im = 0
    If (im/=0) kfam = iabs(k(im,2))
    If (im/=0 .And. mmat==3) Then
      Do il = max(ip-2, im+1), min(ip+2, n)
        If (k(il,3)==im) nm = nm + 1
      End Do
      If (nm/=2 .Or. kfam<=100 .Or. mod(kfam,10)/=1 .Or. mod(kfam/1000,10)/=0) nm = 0
    Else If (im/=0 .And. mmat==46) Then
      msgn = isign(1, k(im,2)*k(ip,2))
      If (kfam>100 .And. mod(kfam/1000,10)==0) msgn = msgn*(-1)**mod(kfam/100, 10)
    End If
  End If
  If (nd==1) Then
    Do j = 1, 4
      p(n+1, j) = p(ip, j)
    End Do
    Goto 510
  End If
  pv(nd, 5) = p(n+nd, 5)
  If (nd>=3) Then
    wtmax = 1./wtcor(nd-2)
    pmax = pv(1, 5) - ps + p(n+nd, 5)
    pmin = 0.
    Do il = nd - 1, 1, -1
      pmax = pmax + p(n+il, 5)
      pmin = pmin + p(n+il+1, 5)
      wtmax = wtmax*pawt(pmax, pmin, p(n+il,5))
    End Do
  End If
  310 If (nd==2) Then
  Else If (mmat==2) Then
    pmes = 4.*pmas(11, 1)**2
    pmrho2 = pmas(131, 1)**2
    pgrho2 = pmas(131, 2)**2
    320 pmst = pmes*(p(ip,5)**2/pmes)**rlu(0)
    wt = (1+0.5*pmes/pmst)*sqrt(max(0.,1.-pmes/pmst))*(1.-pmst/p(ip,5)**2)**3*(1.+pgrho2/pmrho2)/((1.-pmst/pmrho2)**2+pgrho2/pmrho2)
    If (wt<rlu(0)) Goto 320
    pv(2, 5) = max(2.00001*pmas(11,1), sqrt(pmst))
  Else
    330 rord(1) = 1.
    Do il1 = 2, nd - 1
      rsav = rlu(0)
      Do il2 = il1 - 1, 1, -1
        If (rsav<=rord(il2)) Goto 350
        rord(il2+1) = rord(il2)
      End Do
      350 rord(il2+1) = rsav
    End Do
    rord(nd) = 0.
    wt = 1.
    Do il = nd - 1, 1, -1
      pv(il, 5) = pv(il+1, 5) + p(n+il, 5) + (rord(il)-rord(il+1))*(pv(1,5)-ps)
      wt = wt*pawt(pv(il,5), pv(il+1,5), p(n+il,5))
    End Do
    If (wt<rlu(0)*wtmax) Goto 330
  End If
  370 Do il = 1, nd - 1
    pa = pawt(pv(il,5), pv(il+1,5), p(n+il,5))
    ue(3) = 2.*rlu(0) - 1.
    phi = paru(2)*rlu(0)
    ue(1) = sqrt(1.-ue(3)**2)*cos(phi)
    ue(2) = sqrt(1.-ue(3)**2)*sin(phi)
    Do j = 1, 3
      p(n+il, j) = pa*ue(j)
      pv(il+1, j) = -pa*ue(j)
    End Do
    p(n+il, 4) = sqrt(pa**2+p(n+il,5)**2)
    pv(il+1, 4) = sqrt(pa**2+pv(il+1,5)**2)
  End Do
  Do j = 1, 4
    p(n+nd, j) = pv(nd, j)
  End Do
  Do il = nd - 1, 1, -1
    Do j = 1, 3
      be(j) = pv(il, j)/pv(il, 4)
    End Do
    ga = pv(il, 4)/pv(il, 5)
    Do i = n + il, n + nd
      bep = be(1)*p(i, 1) + be(2)*p(i, 2) + be(3)*p(i, 3)
      Do j = 1, 3
        p(i, j) = p(i, j) + ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
      End Do
      p(i, 4) = ga*(p(i,4)+bep)
    End Do
  End Do
  If (mmat==1) Then
    wt = (p(n+1,5)*p(n+2,5)*p(n+3,5))**2 - (p(n+1,5)*four(n+2,n+3))**2 - (p(n+2,5)*four(n+1,n+3))**2 - (p(n+3,5)*four(n+1,n+2))**2 + 2.*four(n+1, n+2)*four(n+1, n+3)*four(n+2, n+3)
    If (max(wt*wtcor(9)/p(ip,5)**6,0.001)<rlu(0)) Goto 310
  Else If (mmat==2) Then
    four12 = four(n+1, n+2)
    four13 = four(n+1, n+3)
    four23 = 0.5*pmst - 0.25*pmes
    wt = (pmst-0.5*pmes)*(four12**2+four13**2) + pmes*(four12*four13+four12**2+four13**2)
    If (wt<rlu(0)*0.25*pmst*(p(ip,5)**2-pmst)**2) Goto 370
  Else If (mmat==3 .And. nm==2) Then
    If ((p(ip,5)**2*four(im,n+1)-four(ip,im)*four(ip,n+1))**2<=rlu(0)*(four(ip,im)**2-(p(ip,5)*p(im,5))**2)*(four(ip,n+1)**2-(p(ip,5)*p(n+1,5))**2)) Goto 370
  Else If (mmat==4) Then
    hx1 = 2.*four(ip, n+1)/p(ip, 5)**2
    hx2 = 2.*four(ip, n+2)/p(ip, 5)**2
    hx3 = 2.*four(ip, n+3)/p(ip, 5)**2
    wt = ((1.-hx1)/(hx2*hx3))**2 + ((1.-hx2)/(hx1*hx3))**2 + ((1.-hx3)/(hx1*hx2))**2
    If (wt<2.*rlu(0)) Goto 310
    If (k(ip+1,2)==22 .And. (1.-hx1)*p(ip,5)**2<4.*parj(32)**2) Goto 310
  Else If (mmat==41) Then
    hx1 = 2.*four(ip, n+1)/p(ip, 5)**2
    If (8.*hx1*(3.-2.*hx1)/9.<rlu(0)) Goto 310
  Else If (mmat>=42 .And. mmat<=44 .And. nd==3) Then
    If (mbst==0) wt = four(ip, n+1)*four(n+2, n+3)
    If (mbst==1) wt = p(ip, 5)*p(n+1, 4)*four(n+2, n+3)
    If (wt<rlu(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) Goto 310
  Else If (mmat>=42 .And. mmat<=44) Then
    Do j = 1, 4
      p(n+np+1, j) = 0.
      Do is = n + 3, n + np
        p(n+np+1, j) = p(n+np+1, j) + p(is, j)
      End Do
    End Do
    If (mbst==0) wt = four(ip, n+1)*four(n+2, n+np+1)
    If (mbst==1) wt = p(ip, 5)*p(n+1, 4)*four(n+2, n+np+1)
    If (wt<rlu(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) Goto 310
  Else If (mmat==46 .And. msgn/=0) Then
    If (msgn>0) wt = four(im, n+1)*four(n+2, ip+1)
    If (msgn<0) wt = four(im, n+2)*four(n+1, ip+1)
    If (wt<rlu(0)*p(im,5)**4/wtcor(10)) Goto 370
  End If
  If (mrem==1) Then
    Do j = 1, 5
      pv(1, j) = pv(1, j)/(1.-pqt)
    End Do
    nd = nd + 1
    mrem = 0
  End If
  If ((mmat==31 .Or. mmat==45) .And. nd==3) Then
    mstj(93) = 1
    pm2 = ulmass(k(n+2,2))
    mstj(93) = 1
    pm3 = ulmass(k(n+3,2))
    If (p(n+2,5)**2+p(n+3,5)**2+2.*four(n+2,n+3)>=(parj(32)+pm2+pm3)**2) Goto 510
    k(n+2, 1) = 1
    kftemp = k(n+2, 2)
    Call lukfdi(kftemp, k(n+3,2), kfldmp, k(n+2,2))
    If (k(n+2,2)==0) Goto 150
    p(n+2, 5) = ulmass(k(n+2,2))
    ps = p(n+1, 5) + p(n+2, 5)
    pv(2, 5) = p(n+2, 5)
    mmat = 0
    nd = 2
    Goto 370
  Else If (mmat==44) Then
    mstj(93) = 1
    pm3 = ulmass(k(n+3,2))
    mstj(93) = 1
    pm4 = ulmass(k(n+4,2))
    If (p(n+3,5)**2+p(n+4,5)**2+2.*four(n+3,n+4)>=(parj(32)+pm3+pm4)**2) Goto 480
    k(n+3, 1) = 1
    kftemp = k(n+3, 2)
    Call lukfdi(kftemp, k(n+4,2), kfldmp, k(n+3,2))
    If (k(n+3,2)==0) Goto 150
    p(n+3, 5) = ulmass(k(n+3,2))
    Do j = 1, 3
      p(n+3, j) = p(n+3, j) + p(n+4, j)
    End Do
    p(n+3, 4) = sqrt(p(n+3,1)**2+p(n+3,2)**2+p(n+3,3)**2+p(n+3,5)**2)
    ha = p(n+1, 4)**2 - p(n+2, 4)**2
    hb = ha - (p(n+1,5)**2-p(n+2,5)**2)
    hc = (p(n+1,1)-p(n+2,1))**2 + (p(n+1,2)-p(n+2,2))**2 + (p(n+1,3)-p(n+2,3))**2
    hd = (pv(1,4)-p(n+3,4))**2
    he = ha**2 - 2.*hd*(p(n+1,4)**2+p(n+2,4)**2) + hd**2
    hf = hd*hc - hb**2
    hg = hd*hc - ha*hb
    hh = (sqrt(hg**2+he*hf)-hg)/(2.*hf)
    Do j = 1, 3
      pcor = hh*(p(n+1,j)-p(n+2,j))
      p(n+1, j) = p(n+1, j) + pcor
      p(n+2, j) = p(n+2, j) - pcor
    End Do
    p(n+1, 4) = sqrt(p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2+p(n+1,5)**2)
    p(n+2, 4) = sqrt(p(n+2,1)**2+p(n+2,2)**2+p(n+2,3)**2+p(n+2,5)**2)
    nd = nd - 1
  End If
  480 If (mmat>=42 .And. mmat<=44 .And. iabs(k(n+1,2))<10) Then
    pmr = sqrt(max(0.,p(n+1,5)**2+p(n+2,5)**2+2.*four(n+1,n+2)))
    mstj(93) = 1
    pm1 = ulmass(k(n+1,2))
    mstj(93) = 1
    pm2 = ulmass(k(n+2,2))
    If (pmr>parj(32)+pm1+pm2) Goto 490
    kfldum = int(1.5+rlu(0))
    Call lukfdi(k(n+1,2), -isign(kfldum,k(n+1,2)), kfldmp, kf1)
    Call lukfdi(k(n+2,2), -isign(kfldum,k(n+2,2)), kfldmp, kf2)
    If (kf1==0 .Or. kf2==0) Goto 150
    psm = ulmass(kf1) + ulmass(kf2)
    If (mmat==42 .And. pmr>parj(64)+psm) Goto 490
    If (mmat>=43 .And. pmr>0.2*parj(32)+psm) Goto 490
    If (nd==4 .Or. kfa==15) Goto 150
    k(n+1, 1) = 1
    kftemp = k(n+1, 2)
    Call lukfdi(kftemp, k(n+2,2), kfldmp, k(n+1,2))
    If (k(n+1,2)==0) Goto 150
    p(n+1, 5) = ulmass(k(n+1,2))
    k(n+2, 2) = k(n+3, 2)
    p(n+2, 5) = p(n+3, 5)
    ps = p(n+1, 5) + p(n+2, 5)
    pv(2, 5) = p(n+3, 5)
    mmat = 0
    nd = 2
    Goto 370
  End If
  490 If (mmat==42 .And. iabs(k(n+1,2))<10) Then
    kflo(1) = k(n+1, 2)
    kflo(2) = k(n+2, 2)
    k(n+1, 1) = k(n+3, 1)
    k(n+1, 2) = k(n+3, 2)
    Do j = 1, 5
      pv(1, j) = p(n+1, j) + p(n+2, j)
      p(n+1, j) = p(n+3, j)
    End Do
    pv(1, 5) = pmr
    n = n + 1
    np = 0
    nq = 2
    ps = 0.
    mstj(93) = 2
    psq = ulmass(kflo(1))
    mstj(93) = 2
    psq = psq + ulmass(kflo(2))
    mmat = 11
    Goto 180
  End If
  510 n = n + nd
  If (mbst==1) Then
    Do j = 1, 3
      be(j) = p(ip, j)/p(ip, 4)
    End Do
    ga = p(ip, 4)/p(ip, 5)
    Do i = nsav + 1, n
      bep = be(1)*p(i, 1) + be(2)*p(i, 2) + be(3)*p(i, 3)
      Do j = 1, 3
        p(i, j) = p(i, j) + ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
      End Do
      p(i, 4) = ga*(p(i,4)+bep)
    End Do
  End If
  Do i = nsav + 1, n
    Do j = 1, 4
      v(i, j) = vdcy(j)
    End Do
    v(i, 5) = 0.
  End Do
  If (mstj(23)>=1 .And. mmat==4 .And. k(nsav+1,2)==21) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+2)
    k(nsav+1, 5) = mstu(5)*(nsav+3)
    k(nsav+2, 4) = mstu(5)*(nsav+3)
    k(nsav+2, 5) = mstu(5)*(nsav+1)
    k(nsav+3, 4) = mstu(5)*(nsav+1)
    k(nsav+3, 5) = mstu(5)*(nsav+2)
    mstj(92) = -(nsav+1)
  Else If (mstj(23)>=1 .And. mmat==4) Then
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+2, 4) = mstu(5)*(nsav+3)
    k(nsav+2, 5) = mstu(5)*(nsav+3)
    k(nsav+3, 4) = mstu(5)*(nsav+2)
    k(nsav+3, 5) = mstu(5)*(nsav+2)
    mstj(92) = nsav + 2
  Else If (mstj(23)>=1 .And. (mmat==32 .Or. mmat==44 .Or. mmat==46) .And. iabs(k(nsav+1,2))<=10 .And. iabs(k(nsav+2,2))<=10) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+2)
    k(nsav+1, 5) = mstu(5)*(nsav+2)
    k(nsav+2, 4) = mstu(5)*(nsav+1)
    k(nsav+2, 5) = mstu(5)*(nsav+1)
    mstj(92) = nsav + 1
  Else If (mstj(23)>=1 .And. mmat==33 .And. iabs(k(nsav+2,2))==21) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    kcp = lucomp(k(nsav+1,2))
    kqp = kchg(kcp, 2)*isign(1, k(nsav+1,2))
    jcon = 4
    If (kqp<0) jcon = 5
    k(nsav+1, jcon) = mstu(5)*(nsav+2)
    k(nsav+2, 9-jcon) = mstu(5)*(nsav+1)
    k(nsav+2, jcon) = mstu(5)*(nsav+3)
    k(nsav+3, 9-jcon) = mstu(5)*(nsav+2)
    mstj(92) = nsav + 1
  Else If (mstj(23)>=1 .And. mmat==33) Then
    k(nsav+1, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+3)
    k(nsav+1, 5) = mstu(5)*(nsav+3)
    k(nsav+3, 4) = mstu(5)*(nsav+1)
    k(nsav+3, 5) = mstu(5)*(nsav+1)
    mstj(92) = nsav + 1
  End If
  If (k(ip,1)==5) k(ip, 1) = 15
  If (k(ip,1)<=10) k(ip, 1) = 11
  k(ip, 4) = nsav + 1
  k(ip, 5) = n
  Return
End Subroutine ludecy
