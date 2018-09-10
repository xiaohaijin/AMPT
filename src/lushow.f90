Subroutine lushow(ip1, ip2, qmax)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension pmth(5, 40), ps(5), pma(4), pmsd(4), iep(4), ipa(4), kfla(4), kfld(4), kfl(4), itry(4), isi(4), isl(4), dp(4), dpt(5, 4)
  If (mstj(41)<=0 .Or. (mstj(41)==1 .And. qmax<=parj(82)) .Or. qmax<=min(parj(82),parj(83)) .Or. mstj(41)>=3) Return
  pmth(1, 21) = ulmass(21)
  pmth(2, 21) = sqrt(pmth(1,21)**2+0.25*parj(82)**2)
  pmth(3, 21) = 2.*pmth(2, 21)
  pmth(4, 21) = pmth(3, 21)
  pmth(5, 21) = pmth(3, 21)
  pmth(1, 22) = ulmass(22)
  pmth(2, 22) = sqrt(pmth(1,22)**2+0.25*parj(83)**2)
  pmth(3, 22) = 2.*pmth(2, 22)
  pmth(4, 22) = pmth(3, 22)
  pmth(5, 22) = pmth(3, 22)
  pmqth1 = parj(82)
  If (mstj(41)==2) pmqth1 = min(parj(82), parj(83))
  pmqth2 = pmth(2, 21)
  If (mstj(41)==2) pmqth2 = min(pmth(2,21), pmth(2,22))
  Do if = 1, 8
    pmth(1, if) = ulmass(if)
    pmth(2, if) = sqrt(pmth(1,if)**2+0.25*pmqth1**2)
    pmth(3, if) = pmth(2, if) + pmqth2
    pmth(4, if) = sqrt(pmth(1,if)**2+0.25*parj(82)**2) + pmth(2, 21)
    pmth(5, if) = sqrt(pmth(1,if)**2+0.25*parj(83)**2) + pmth(2, 22)
  End Do
  pt2min = max(0.5*parj(82), 1.1*parj(81))**2
  alams = parj(81)**2
  alfm = log(pt2min/alams)
  m3jc = 0
  If (ip1>0 .And. ip1<=min(n,mstu(4)-mstu(32)) .And. ip2==0) Then
    npa = 1
    ipa(1) = ip1
  Else If (min(ip1,ip2)>0 .And. max(ip1,ip2)<=min(n,mstu(4)-mstu(32))) Then
    npa = 2
    ipa(1) = ip1
    ipa(2) = ip2
  Else If (ip1>0 .And. ip1<=min(n,mstu(4)-mstu(32)) .And. ip2<0 .And. ip2>=-3) Then
    npa = iabs(ip2)
    Do i = 1, npa
      ipa(i) = ip1 + i - 1
    End Do
  Else
    Call luerrm(12, '(LUSHOW:) failed to reconstruct showering system')
    If (mstu(21)>=1) Return
  End If
  irej = 0
  Do j = 1, 5
    ps(j) = 0.
  End Do
  pm = 0.
  Do i = 1, npa
    kfla(i) = iabs(k(ipa(i),2))
    pma(i) = p(ipa(i), 5)
    If (kfla(i)/=0 .And. (kfla(i)<=8 .Or. kfla(i)==21)) pma(i) = pmth(3, kfla(i))
    pm = pm + pma(i)
    If (kfla(i)==0 .Or. (kfla(i)>8 .And. kfla(i)/=21) .Or. pma(i)>qmax) irej = irej + 1
    Do j = 1, 4
      ps(j) = ps(j) + p(ipa(i), j)
    End Do
  End Do
  If (irej==npa) Return
  ps(5) = sqrt(max(0.,ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2))
  If (npa==1) ps(5) = ps(4)
  If (ps(5)<=pm+pmqth1) Return
  If (npa==2 .And. mstj(47)>=1) Then
    If (kfla(1)>=1 .And. kfla(1)<=8 .And. kfla(2)>=1 .And. kfla(2)<=8) m3jc = 1
    If (mstj(47)>=2) m3jc = 1
  End If
  ns = n
  If (n>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  If (npa>=2) Then
    k(n+1, 1) = 11
    k(n+1, 2) = 21
    k(n+1, 3) = 0
    k(n+1, 4) = 0
    k(n+1, 5) = 0
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = 0.
    p(n+1, 4) = ps(5)
    p(n+1, 5) = ps(5)
    v(n+1, 5) = ps(5)**2
    n = n + 1
  End If
  nep = npa
  im = ns
  If (npa==1) im = ns - 1
  140 im = im + 1
  If (n>ns) Then
    If (im>n) Goto 380
    kflm = iabs(k(im,2))
    If (kflm==0 .Or. (kflm>8 .And. kflm/=21)) Goto 140
    If (p(im,5)<pmth(2,kflm)) Goto 140
    igm = k(im, 3)
  Else
    igm = -1
  End If
  If (n+nep>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  iau = 0
  If (igm>0) Then
    If (k(im-1,3)==igm) iau = im - 1
    If (n>=im+1 .And. k(im+1,3)==igm) iau = im + 1
  End If
  If (igm>=0) Then
    k(im, 4) = n + 1
    Do i = 1, nep
      k(n+i, 3) = im
    End Do
  Else
    k(n+1, 3) = ipa(1)
  End If
  If (igm<=0) Then
    Do i = 1, nep
      k(n+i, 2) = k(ipa(i), 2)
    End Do
  Else If (kflm/=21) Then
    k(n+1, 2) = k(im, 2)
    k(n+2, 2) = k(im, 5)
  Else If (k(im,5)==21) Then
    k(n+1, 2) = 21
    k(n+2, 2) = 21
  Else
    k(n+1, 2) = k(im, 5)
    k(n+2, 2) = -k(im, 5)
  End If
  Do ip = 1, nep
    k(n+ip, 1) = 3
    k(n+ip, 4) = 0
    k(n+ip, 5) = 0
    kfld(ip) = iabs(k(n+ip,2))
    itry(ip) = 0
    isl(ip) = 0
    isi(ip) = 0
    If (kfld(ip)>0 .And. (kfld(ip)<=8 .Or. kfld(ip)==21)) isi(ip) = 1
  End Do
  islm = 0
  If (igm<=0) Then
    Do i = 1, npa
      If (npa>=3) p(n+i, 4) = (ps(4)*p(ipa(i),4)-ps(1)*p(ipa(i),1)-ps(2)*p(ipa(i),2)-ps(3)*p(ipa(i),3))/ps(5)
      p(n+i, 5) = min(qmax, ps(5))
      If (npa>=3) p(n+i, 5) = min(p(n+i,5), p(n+i,4))
      If (isi(i)==0) p(n+i, 5) = p(ipa(i), 5)
    End Do
  Else
    If (mstj(43)<=2) pem = v(im, 2)
    If (mstj(43)>=3) pem = p(im, 4)
    p(n+1, 5) = min(p(im,5), v(im,1)*pem)
    p(n+2, 5) = min(p(im,5), (1.-v(im,1))*pem)
    If (k(n+2,2)==22) p(n+2, 5) = pmth(1, 22)
  End If
  Do i = 1, nep
    pmsd(i) = p(n+i, 5)
    If (isi(i)==1) Then
      If (p(n+i,5)<=pmth(3,kfld(i))) p(n+i, 5) = pmth(1, kfld(i))
    End If
    v(n+i, 5) = p(n+i, 5)**2
  End Do
  200 inum = 0
  If (nep==1) inum = 1
  Do i = 1, nep
    If (inum==0 .And. isl(i)==1) inum = i
  End Do
  Do i = 1, nep
    If (inum==0 .And. itry(i)==0 .And. isi(i)==1) Then
      If (p(n+i,5)>=pmth(2,kfld(i))) inum = i
    End If
  End Do
  If (inum==0) Then
    rmax = 0.
    Do i = 1, nep
      If (isi(i)==1 .And. pmsd(i)>=pmqth2) Then
        rpm = p(n+i, 5)/pmsd(i)
        If (rpm>rmax .And. p(n+i,5)>=pmth(2,kfld(i))) Then
          rmax = rpm
          inum = i
        End If
      End If
    End Do
  End If
  inum = max(1, inum)
  iep(1) = n + inum
  Do i = 2, nep
    iep(i) = iep(i-1) + 1
    If (iep(i)>n+nep) iep(i) = n + 1
  End Do
  Do i = 1, nep
    kfl(i) = iabs(k(iep(i),2))
  End Do
  itry(inum) = itry(inum) + 1
  If (itry(inum)>200) Then
    Call luerrm(14, '(LUSHOW:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  z = 0.5
  If (kfl(1)==0 .Or. (kfl(1)>8 .And. kfl(1)/=21)) Goto 300
  If (p(iep(1),5)<pmth(2,kfl(1))) Goto 300
  If (nep==1) Then
    pmed = ps(4)
  Else If (igm==0 .Or. mstj(43)<=2) Then
    pmed = p(im, 5)
  Else
    If (inum==1) pmed = v(im, 1)*pem
    If (inum==2) pmed = (1.-v(im,1))*pem
  End If
  If (mod(mstj(43),2)==1) Then
    zc = pmth(2, 21)/pmed
    zce = pmth(2, 22)/pmed
  Else
    zc = 0.5*(1.-sqrt(max(0.,1.-(2.*pmth(2,21)/pmed)**2)))
    If (zc<1E-4) zc = (pmth(2,21)/pmed)**2
    zce = 0.5*(1.-sqrt(max(0.,1.-(2.*pmth(2,22)/pmed)**2)))
    If (zce<1E-4) zce = (pmth(2,22)/pmed)**2
  End If
  zc = min(zc, 0.491)
  zce = min(zce, 0.491)
  If ((mstj(41)==1 .And. zc>0.49) .Or. (mstj(41)==2 .And. min(zc,zce)>0.49)) Then
    p(iep(1), 5) = pmth(1, kfl(1))
    v(iep(1), 5) = p(iep(1), 5)**2
    Goto 300
  End If
  If (mstj(49)==0 .And. kfl(1)==21) Then
    fbr = 6.*log((1.-zc)/zc) + mstj(45)*(0.5-zc)
  Else If (mstj(49)==0) Then
    fbr = (8./3.)*log((1.-zc)/zc)
  Else If (mstj(49)==1 .And. kfl(1)==21) Then
    fbr = (parj(87)+mstj(45)*parj(88))*(1.-2.*zc)
  Else If (mstj(49)==1) Then
    fbr = (1.-2.*zc)/3.
    If (igm==0 .And. m3jc==1) fbr = 4.*fbr
  Else If (kfl(1)==21) Then
    fbr = 6.*mstj(45)*(0.5-zc)
  Else
    fbr = 2.*log((1.-zc)/zc)
  End If
  If (mstj(41)==2 .And. kfl(1)>=1 .And. kfl(1)<=8) fbre = (kchg(kfl(1),1)/3.)**2*2.*log((1.-zce)/zce)
  260 pms = v(iep(1), 5)
  If (igm>=0) Then
    pm2 = 0.
    Do i = 2, nep
      pm = p(iep(i), 5)
      If (kfl(i)>0 .And. (kfl(i)<=8 .Or. kfl(i)==21)) pm = pmth(2, kfl(i))
      pm2 = pm2 + pm
    End Do
    pms = min(pms, (p(im,5)-pm2)**2)
  End If
  b0 = 27./6.
  Do if = 4, mstj(45)
    If (pms>4.*pmth(2,if)**2) b0 = (33.-2.*if)/6.
  End Do
  If (mstj(44)<=0) Then
    pmsqcd = pms*exp(max(-100.,log(rlu(0))*paru(2)/(paru(111)*fbr)))
  Else If (mstj(44)==1) Then
    pmsqcd = 4.*alams*(0.25*pms/alams)**(rlu(0)**(b0/fbr))
  Else
    pmsqcd = pms*rlu(0)**(alfm*b0/fbr)
  End If
  If (zc>0.49 .Or. pmsqcd<=pmth(4,kfl(1))**2) pmsqcd = pmth(2, kfl(1))**2
  v(iep(1), 5) = pmsqcd
  mce = 1
  If (mstj(41)==2 .And. kfl(1)>=1 .And. kfl(1)<=8) Then
    pmsqed = pms*exp(max(-100.,log(rlu(0))*paru(2)/(paru(101)*fbre)))
    If (zce>0.49 .Or. pmsqed<=pmth(5,kfl(1))**2) pmsqed = pmth(2, kfl(1))**2
    If (pmsqed>pmsqcd) Then
      v(iep(1), 5) = pmsqed
      mce = 2
    End If
  End If
  p(iep(1), 5) = sqrt(v(iep(1),5))
  If (p(iep(1),5)<=pmth(3,kfl(1))) Then
    p(iep(1), 5) = pmth(1, kfl(1))
    v(iep(1), 5) = p(iep(1), 5)**2
    Goto 300
  End If
  If (mce==2) Then
    z = 1. - (1.-zce)*(zce/(1.-zce))**rlu(0)
    If (1.+z**2<2.*rlu(0)) Goto 260
    k(iep(1), 5) = 22
  Else If (mstj(49)/=1 .And. kfl(1)/=21) Then
    z = 1. - (1.-zc)*(zc/(1.-zc))**rlu(0)
    If (1.+z**2<2.*rlu(0)) Goto 260
    k(iep(1), 5) = 21
  Else If (mstj(49)==0 .And. mstj(45)*(0.5-zc)<rlu(0)*fbr) Then
    z = (1.-zc)*(zc/(1.-zc))**rlu(0)
    If (rlu(0)>0.5) z = 1. - z
    If ((1.-z*(1.-z))**2<rlu(0)) Goto 260
    k(iep(1), 5) = 21
  Else If (mstj(49)/=1) Then
    z = zc + (1.-2.*zc)*rlu(0)
    If (z**2+(1.-z)**2<rlu(0)) Goto 260
    kflb = 1 + int(mstj(45)*rlu(0))
    pmq = 4.*pmth(2, kflb)**2/v(iep(1), 5)
    If (pmq>=1.) Goto 260
    pmq0 = 4.*pmth(2, 21)**2/v(iep(1), 5)
    If (mod(mstj(43),2)==0 .And. (1.+0.5*pmq)*sqrt(1.-pmq)<rlu(0)*(1.+0.5*pmq0)*sqrt(1.-pmq0)) Goto 260
    k(iep(1), 5) = kflb
  Else If (kfl(1)/=21) Then
    z = 1. - sqrt(zc**2+rlu(0)*(1.-2.*zc))
    k(iep(1), 5) = 21
  Else If (rlu(0)*(parj(87)+mstj(45)*parj(88))<=parj(87)) Then
    z = zc + (1.-2.*zc)*rlu(0)
    k(iep(1), 5) = 21
  Else
    z = zc + (1.-2.*zc)*rlu(0)
    kflb = 1 + int(mstj(45)*rlu(0))
    pmq = 4.*pmth(2, kflb)**2/v(iep(1), 5)
    If (pmq>=1.) Goto 260
    k(iep(1), 5) = kflb
  End If
  If (mce==1 .And. mstj(44)>=2) Then
    If (z*(1.-z)*v(iep(1),5)<pt2min) Goto 260
    If (alfm/log(v(iep(1),5)*z*(1.-z)/alams)<rlu(0)) Goto 260
  End If
  If (kfl(1)==21) Then
    kflgd1 = iabs(k(iep(1),5))
    kflgd2 = kflgd1
  Else
    kflgd1 = kfl(1)
    kflgd2 = iabs(k(iep(1),5))
  End If
  If (nep==1) Then
    ped = ps(4)
  Else If (nep>=3) Then
    ped = p(iep(1), 4)
  Else If (igm==0 .Or. mstj(43)<=2) Then
    ped = 0.5*(v(im,5)+v(iep(1),5)-pm2**2)/p(im, 5)
  Else
    If (iep(1)==n+1) ped = v(im, 1)*pem
    If (iep(1)==n+2) ped = (1.-v(im,1))*pem
  End If
  If (mod(mstj(43),2)==1) Then
    pmqth3 = 0.5*parj(82)
    If (kflgd2==22) pmqth3 = 0.5*parj(83)
    pmq1 = (pmth(1,kflgd1)**2+pmqth3**2)/v(iep(1), 5)
    pmq2 = (pmth(1,kflgd2)**2+pmqth3**2)/v(iep(1), 5)
    zd = sqrt(max(0.,(1.-v(iep(1),5)/ped**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
    zh = 1. + pmq1 - pmq2
  Else
    zd = sqrt(max(0.,1.-v(iep(1),5)/ped**2))
    zh = 1.
  End If
  zl = 0.5*(zh-zd)
  zu = 0.5*(zh+zd)
  If (z<zl .Or. z>zu) Goto 260
  If (kfl(1)==21) v(iep(1), 3) = log(zu*(1.-zl)/max(1E-20,zl*(1.-zu)))
  If (kfl(1)/=21) v(iep(1), 3) = log((1.-zl)/max(1E-10,1.-zu))
  If (igm==0 .And. m3jc==1) Then
    x1 = z*(1.+v(iep(1),5)/v(ns+1,5))
    x2 = 1. - v(iep(1), 5)/v(ns+1, 5)
    x3 = (1.-x1) + (1.-x2)
    If (mce==2) Then
      ki1 = k(ipa(inum), 2)
      ki2 = k(ipa(3-inum), 2)
      qf1 = kchg(iabs(ki1), 1)*isign(1, ki1)/3.
      qf2 = kchg(iabs(ki2), 1)*isign(1, ki2)/3.
      wshow = qf1**2*(1.-x1)/x3*(1.+(x1/(2.-x2))**2) + qf2**2*(1.-x2)/x3*(1.+(x2/(2.-x1))**2)
      wme = (qf1*(1.-x1)/x3-qf2*(1.-x2)/x3)**2*(x1**2+x2**2)
    Else If (mstj(49)/=1) Then
      wshow = 1. + (1.-x1)/x3*(x1/(2.-x2))**2 + (1.-x2)/x3*(x2/(2.-x1))**2
      wme = x1**2 + x2**2
    Else
      wshow = 4.*x3*((1.-x1)/(2.-x2)**2+(1.-x2)/(2.-x1)**2)
      wme = x3**2
    End If
    If (wme<rlu(0)*wshow) Goto 260
  Else If (mce==1 .And. igm>0 .And. mstj(42)>=2) Then
    maom = 1
    zm = v(im, 1)
    If (iep(1)==n+2) zm = 1. - v(im, 1)
    the2id = z*(1.-z)*(zm*p(im,4))**2/v(iep(1), 5)
    iaom = im
    290 If (k(iaom,5)==22) Then
      iaom = k(iaom, 3)
      If (k(iaom,3)<=ns) maom = 0
      If (maom==1) Goto 290
    End If
    If (maom==1) Then
      the2im = v(iaom, 1)*(1.-v(iaom,1))*p(iaom, 4)**2/v(iaom, 5)
      If (the2id<the2im) Goto 260
    End If
  End If
  If (mstj(48)==1) Then
    If (nep==1 .And. im==ns) Then
      the2id = z*(1.-z)*ps(4)**2/v(iep(1), 5)
      If (the2id<1./parj(85)**2) Goto 260
    Else If (nep==2 .And. iep(1)==ns+2) Then
      the2id = z*(1.-z)*(0.5*p(im,4))**2/v(iep(1), 5)
      If (the2id<1./parj(85)**2) Goto 260
    Else If (nep==2 .And. iep(1)==ns+3) Then
      the2id = z*(1.-z)*(0.5*p(im,4))**2/v(iep(1), 5)
      If (the2id<1./parj(86)**2) Goto 260
    End If
  End If
  300 v(iep(1), 1) = z
  isl(1) = 0
  isl(2) = 0
  If (nep==1) Goto 330
  If (nep==2 .And. p(iep(1),5)+p(iep(2),5)>=p(im,5)) Goto 200
  Do i = 1, nep
    If (itry(i)==0 .And. kfld(i)>0 .And. (kfld(i)<=8 .Or. kfld(i)==21)) Then
      If (p(n+i,5)>=pmth(2,kfld(i))) Goto 200
    End If
  End Do
  If (nep==3) Then
    pa1s = (p(n+1,4)+p(n+1,5))*(p(n+1,4)-p(n+1,5))
    pa2s = (p(n+2,4)+p(n+2,5))*(p(n+2,4)-p(n+2,5))
    pa3s = (p(n+3,4)+p(n+3,5))*(p(n+3,4)-p(n+3,5))
    pts = 0.25*(2.*pa1s*pa2s+2.*pa1s*pa3s+2.*pa2s*pa3s-pa1s**2-pa2s**2-pa3s**2)/pa1s
    If (pts<=0.) Goto 200
  Else If (igm==0 .Or. mstj(43)<=2 .Or. mod(mstj(43),2)==0) Then
    Do i1 = n + 1, n + 2
      kflda = iabs(k(i1,2))
      If (kflda==0 .Or. (kflda>8 .And. kflda/=21)) Goto 320
      If (p(i1,5)<pmth(2,kflda)) Goto 320
      If (kflda==21) Then
        kflgd1 = iabs(k(i1,5))
        kflgd2 = kflgd1
      Else
        kflgd1 = kflda
        kflgd2 = iabs(k(i1,5))
      End If
      i2 = 2*n + 3 - i1
      If (igm==0 .Or. mstj(43)<=2) Then
        ped = 0.5*(v(im,5)+v(i1,5)-v(i2,5))/p(im, 5)
      Else
        If (i1==n+1) zm = v(im, 1)
        If (i1==n+2) zm = 1. - v(im, 1)
        pml = sqrt((v(im,5)-v(n+1,5)-v(n+2,5))**2-4.*v(n+1,5)*v(n+2,5))
        ped = pem*(0.5*(v(im,5)-pml+v(i1,5)-v(i2,5))+pml*zm)/v(im, 5)
      End If
      If (mod(mstj(43),2)==1) Then
        pmqth3 = 0.5*parj(82)
        If (kflgd2==22) pmqth3 = 0.5*parj(83)
        pmq1 = (pmth(1,kflgd1)**2+pmqth3**2)/v(i1, 5)
        pmq2 = (pmth(1,kflgd2)**2+pmqth3**2)/v(i1, 5)
        zd = sqrt(max(0.,(1.-v(i1,5)/ped**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
        zh = 1. + pmq1 - pmq2
      Else
        zd = sqrt(max(0.,1.-v(i1,5)/ped**2))
        zh = 1.
      End If
      zl = 0.5*(zh-zd)
      zu = 0.5*(zh+zd)
      If (i1==n+1 .And. (v(i1,1)<zl .Or. v(i1,1)>zu)) isl(1) = 1
      If (i1==n+2 .And. (v(i1,1)<zl .Or. v(i1,1)>zu)) isl(2) = 1
      If (kflda==21) v(i1, 4) = log(zu*(1.-zl)/max(1E-20,zl*(1.-zu)))
      If (kflda/=21) v(i1, 4) = log((1.-zl)/max(1E-10,1.-zu))
    320 End Do
    If (isl(1)==1 .And. isl(2)==1 .And. islm/=0) Then
      isl(3-islm) = 0
      islm = 3 - islm
    Else If (isl(1)==1 .And. isl(2)==1) Then
      zdr1 = max(0., v(n+1,3)/v(n+1,4)-1.)
      zdr2 = max(0., v(n+2,3)/v(n+2,4)-1.)
      If (zdr2>rlu(0)*(zdr1+zdr2)) isl(1) = 0
      If (isl(1)==1) isl(2) = 0
      If (isl(1)==0) islm = 1
      If (isl(2)==0) islm = 2
    End If
    If (isl(1)==1 .Or. isl(2)==1) Goto 200
  End If
  If (igm>0 .And. mod(mstj(43),2)==1 .And. (p(n+1,5)>=pmth(2,kfld(1)) .Or. p(n+2,5)>=pmth(2,kfld(2)))) Then
    pmq1 = v(n+1, 5)/v(im, 5)
    pmq2 = v(n+2, 5)/v(im, 5)
    zd = sqrt(max(0.,(1.-v(im,5)/pem**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
    zh = 1. + pmq1 - pmq2
    zl = 0.5*(zh-zd)
    zu = 0.5*(zh+zd)
    If (v(im,1)<zl .Or. v(im,1)>zu) Goto 200
  End If
  330 mazip = 0
  mazic = 0
  If (nep==1) Then
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,(p(ipa(1),4)+p(n+1,5))*(p(ipa(1),4)-p(n+1,5))))
    p(n+1, 4) = p(ipa(1), 4)
    v(n+1, 2) = p(n+1, 4)
  Else If (igm==0 .And. nep==2) Then
    ped1 = 0.5*(v(im,5)+v(n+1,5)-v(n+2,5))/p(im, 5)
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,(ped1+p(n+1,5))*(ped1-p(n+1,5))))
    p(n+1, 4) = ped1
    p(n+2, 1) = 0.
    p(n+2, 2) = 0.
    p(n+2, 3) = -p(n+1, 3)
    p(n+2, 4) = p(im, 5) - ped1
    v(n+1, 2) = p(n+1, 4)
    v(n+2, 2) = p(n+2, 4)
  Else If (nep==3) Then
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,pa1s))
    p(n+2, 1) = sqrt(pts)
    p(n+2, 2) = 0.
    p(n+2, 3) = 0.5*(pa3s-pa2s-pa1s)/p(n+1, 3)
    p(n+3, 1) = -p(n+2, 1)
    p(n+3, 2) = 0.
    p(n+3, 3) = -(p(n+1,3)+p(n+2,3))
    v(n+1, 2) = p(n+1, 4)
    v(n+2, 2) = p(n+2, 4)
    v(n+3, 2) = p(n+3, 4)
  Else
    zm = v(im, 1)
    pzm = sqrt(max(0.,(pem+p(im,5))*(pem-p(im,5))))
    pmls = (v(im,5)-v(n+1,5)-v(n+2,5))**2 - 4.*v(n+1, 5)*v(n+2, 5)
    If (pzm<=0.) Then
      pts = 0.
    Else If (mod(mstj(43),2)==1) Then
      pts = (pem**2*(zm*(1.-zm)*v(im,5)-(1.-zm)*v(n+1,5)-zm*v(n+2,5))-0.25*pmls)/pzm**2
    Else
      pts = pmls*(zm*(1.-zm)*pem**2/v(im,5)-0.25)/pzm**2
    End If
    pt = sqrt(max(0.,pts))
    hazip = 0.
    If (mstj(49)/=1 .And. mod(mstj(46),2)==1 .And. k(im,2)==21 .And. iau/=0) Then
      If (k(igm,3)/=0) mazip = 1
      zau = v(igm, 1)
      If (iau==im+1) zau = 1. - v(igm, 1)
      If (mazip==0) zau = 0.
      If (k(igm,2)/=21) Then
        hazip = 2.*zau/(1.+zau**2)
      Else
        hazip = (zau/(1.-zau*(1.-zau)))**2
      End If
      If (k(n+1,2)/=21) Then
        hazip = hazip*(-2.*zm*(1.-zm))/(1.-2.*zm*(1.-zm))
      Else
        hazip = hazip*(zm*(1.-zm)/(1.-zm*(1.-zm)))**2
      End If
    End If
    hazic = 0.
    If (mstj(46)>=2 .And. (k(n+1,2)==21 .Or. k(n+2,2)==21) .And. iau/=0) Then
      If (k(igm,3)/=0) mazic = n + 1
      If (k(igm,3)/=0 .And. k(n+1,2)/=21) mazic = n + 2
      If (k(igm,3)/=0 .And. k(n+1,2)==21 .And. k(n+2,2)==21 .And. zm>0.5) mazic = n + 2
      If (k(iau,2)==22) mazic = 0
      zs = zm
      If (mazic==n+2) zs = 1. - zm
      zgm = v(igm, 1)
      If (iau==im-1) zgm = 1. - v(igm, 1)
      If (mazic==0) zgm = 1.
      hazic = (p(im,5)/p(igm,5))*sqrt((1.-zs)*(1.-zgm)/(zs*zgm))
      hazic = min(0.95, hazic)
    End If
  End If
  340 If (nep==2 .And. igm>0) Then
    If (mod(mstj(43),2)==1) Then
      p(n+1, 4) = pem*v(im, 1)
    Else
      p(n+1, 4) = pem*(0.5*(v(im,5)-sqrt(pmls)+v(n+1,5)-v(n+2,5))+sqrt(pmls)*zm)/v(im, 5)
    End If
    phi = paru(2)*rlu(0)
    p(n+1, 1) = pt*cos(phi)
    p(n+1, 2) = pt*sin(phi)
    If (pzm>0.) Then
      p(n+1, 3) = 0.5*(v(n+2,5)-v(n+1,5)-v(im,5)+2.*pem*p(n+1,4))/pzm
    Else
      p(n+1, 3) = 0.
    End If
    p(n+2, 1) = -p(n+1, 1)
    p(n+2, 2) = -p(n+1, 2)
    p(n+2, 3) = pzm - p(n+1, 3)
    p(n+2, 4) = pem - p(n+1, 4)
    If (mstj(43)<=2) Then
      v(n+1, 2) = (pem*p(n+1,4)-pzm*p(n+1,3))/p(im, 5)
      v(n+2, 2) = (pem*p(n+2,4)-pzm*p(n+2,3))/p(im, 5)
    End If
  End If
  If (igm>0) Then
    If (mstj(43)<=2) Then
      bex = p(igm, 1)/p(igm, 4)
      bey = p(igm, 2)/p(igm, 4)
      bez = p(igm, 3)/p(igm, 4)
      ga = p(igm, 4)/p(igm, 5)
      gabep = ga*(ga*(bex*p(im,1)+bey*p(im,2)+bez*p(im,3))/(1.+ga)-p(im,4))
    Else
      bex = 0.
      bey = 0.
      bez = 0.
      ga = 1.
      gabep = 0.
    End If
    the = ulangl(p(im,3)+gabep*bez, sqrt((p(im,1)+gabep*bex)**2+(p(im,2)+gabep*bey)**2))
    phi = ulangl(p(im,1)+gabep*bex, p(im,2)+gabep*bey)
    Do i = n + 1, n + 2
      dp(1) = dble(cos(the)*cos(phi)*p(i,1)-sin(phi)*p(i,2)+sin(the)*cos(phi)*p(i,3))
      dp(2) = dble(cos(the)*sin(phi)*p(i,1)+cos(phi)*p(i,2)+sin(the)*sin(phi)*p(i,3))
      dp(3) = dble(-sin(the)*p(i,1)+cos(the)*p(i,3))
      dp(4) = dble(p(i,4))
      dbp = dble(bex)*dp(1) + dble(bey)*dp(2) + dble(bez)*dp(3)
      dgabp = dble(ga)*(dble(ga)*dbp/(1D0+dble(ga))+dp(4))
      p(i, 1) = sngl(dp(1)+dgabp*dble(bex))
      p(i, 2) = sngl(dp(2)+dgabp*dble(bey))
      p(i, 3) = sngl(dp(3)+dgabp*dble(bez))
      p(i, 4) = ga*sngl(dp(4)+dbp)
    End Do
  End If
  If (mazip/=0 .Or. mazic/=0) Then
    Do j = 1, 3
      dpt(1, j) = dble(p(im,j))
      dpt(2, j) = dble(p(iau,j))
      dpt(3, j) = dble(p(n+1,j))
    End Do
    dpma = dpt(1, 1)*dpt(2, 1) + dpt(1, 2)*dpt(2, 2) + dpt(1, 3)*dpt(2, 3)
    dpmd = dpt(1, 1)*dpt(3, 1) + dpt(1, 2)*dpt(3, 2) + dpt(1, 3)*dpt(3, 3)
    dpmm = dpt(1, 1)**2 + dpt(1, 2)**2 + dpt(1, 3)**2
    Do j = 1, 3
      dpt(4, j) = dpt(2, j) - dpma*dpt(1, j)/dpmm
      dpt(5, j) = dpt(3, j) - dpmd*dpt(1, j)/dpmm
    End Do
    dpt(4, 4) = dsqrt(dpt(4,1)**2+dpt(4,2)**2+dpt(4,3)**2)
    dpt(5, 4) = dsqrt(dpt(5,1)**2+dpt(5,2)**2+dpt(5,3)**2)
    If (sngl(min(dpt(4,4),dpt(5,4)))>(0.1*parj(82))) Then
      cad = sngl((dpt(4,1)*dpt(5,1)+dpt(4,2)*dpt(5,2)+dpt(4,3)*dpt(5,3))/(dpt(4,4)*dpt(5,4)))
      If (mazip/=0) Then
        If (1.+hazip*(2.*cad**2-1.)<rlu(0)*(1.+abs(hazip))) Goto 340
      End If
      If (mazic/=0) Then
        If (mazic==n+2) cad = -cad
        If ((1.-hazic)*(1.-hazic*cad)/(1.+hazic**2-2.*hazic*cad)<rlu(0)) Goto 340
      End If
    End If
  End If
  If (igm>=0) k(im, 1) = 14
  n = n + nep
  nep = 2
  If (n>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) n = ns
    If (mstu(21)>=1) Return
  End If
  Goto 140
  380 If (npa>=2) Then
    k(ns+1, 1) = 11
    k(ns+1, 2) = 94
    k(ns+1, 3) = ip1
    If (ip2>0 .And. ip2<ip1) k(ns+1, 3) = ip2
    k(ns+1, 4) = ns + 2
    k(ns+1, 5) = ns + 1 + npa
    iim = 1
  Else
    iim = 0
  End If
  Do i = ns + 1 + iim, n
    If (k(i,1)<=10 .And. k(i,2)==22) Then
      k(i, 1) = 1
    Else If (k(i,1)<=10) Then
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5))
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5))
    Else If (k(mod(k(i,4),mstu(5))+1,2)/=22) Then
      id1 = mod(k(i,4), mstu(5))
      If (k(i,2)>=1 .And. k(i,2)<=8) id1 = mod(k(i,4), mstu(5)) + 1
      id2 = 2*mod(k(i,4), mstu(5)) + 1 - id1
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5)) + id1
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5)) + id2
      k(id1, 4) = k(id1, 4) + mstu(5)*i
      k(id1, 5) = k(id1, 5) + mstu(5)*id2
      k(id2, 4) = k(id2, 4) + mstu(5)*id1
      k(id2, 5) = k(id2, 5) + mstu(5)*i
    Else
      id1 = mod(k(i,4), mstu(5))
      id2 = id1 + 1
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5)) + id1
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5)) + id1
      k(id1, 4) = k(id1, 4) + mstu(5)*i
      k(id1, 5) = k(id1, 5) + mstu(5)*i
      k(id2, 4) = 0
      k(id2, 5) = 0
    End If
  End Do
  If (npa>=2) Then
    bex = ps(1)/ps(4)
    bey = ps(2)/ps(4)
    bez = ps(3)/ps(4)
    ga = ps(4)/ps(5)
    gabep = ga*(ga*(bex*p(ipa(1),1)+bey*p(ipa(1),2)+bez*p(ipa(1),3))/(1.+ga)-p(ipa(1),4))
  Else
    bex = 0.
    bey = 0.
    bez = 0.
    gabep = 0.
  End If
  the = ulangl(p(ipa(1),3)+gabep*bez, sqrt((p(ipa(1),1)+gabep*bex)**2+(p(ipa(1),2)+gabep*bey)**2))
  phi = ulangl(p(ipa(1),1)+gabep*bex, p(ipa(1),2)+gabep*bey)
  If (npa==3) Then
    chi = ulangl(cos(the)*cos(phi)*(p(ipa(2),1)+gabep*bex)+cos(the)*sin(phi)*(p(ipa(2),2)+gabep*bey)-sin(the)*(p(ipa(2),3)+gabep*bez), -sin(phi)*(p(ipa(2),1)+gabep*bex)+cos(phi)*(p(ipa(2),2)+gabep*bey))
    Call ludbrb(ns+1, n, 0., chi, 0D0, 0D0, 0D0)
  End If
  dbex = dble(bex)
  dbey = dble(bey)
  dbez = dble(bez)
  Call ludbrb(ns+1, n, the, phi, dbex, dbey, dbez)
  Do i = ns + 1, n
    Do j = 1, 5
      v(i, j) = v(ip1, j)
    End Do
  End Do
  If (n==ns+npa+iim) Then
    n = ns
  Else
    Do ip = 1, npa
      k(ipa(ip), 1) = 14
      k(ipa(ip), 4) = k(ipa(ip), 4) + ns + iim + ip
      k(ipa(ip), 5) = k(ipa(ip), 5) + ns + iim + ip
      k(ns+iim+ip, 3) = ipa(ip)
      If (iim==1 .And. mstu(16)/=2) k(ns+iim+ip, 3) = ns + 1
      k(ns+iim+ip, 4) = mstu(5)*ipa(ip) + k(ns+iim+ip, 4)
      k(ns+iim+ip, 5) = mstu(5)*ipa(ip) + k(ns+iim+ip, 5)
    End Do
  End If
  Return
End Subroutine lushow
