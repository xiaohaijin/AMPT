Subroutine pysspa(ipu1, ipu2)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Dimension kfls(4), is(2), xs(2), zs(2), q2s(2), tevs(2), robo(5), xfs(2, -6:6), xfa(-6:6), xfb(-6:6), xfn(-6:6), wtap(-6:6), wtsf(-6:6), the2(2), alam(2), dq2(3), dpc(3), dpd(4), dpb(4)
  ipus1 = ipu1
  ipus2 = ipu2
  isub = mint(1)
  q2e = vint(52)
  If (iset(isub)==1) Then
    q2e = q2e/parp(67)
  Else If (iset(isub)==3 .Or. iset(isub)==4) Then
    q2e = pmas(23, 1)**2
    If (isub==8 .Or. isub==76 .Or. isub==77) q2e = pmas(24, 1)**2
  End If
  tmax = log(parp(67)*parp(63)*q2e/parp(61)**2)
  If (parp(67)*q2e<max(parp(62)**2,2.*parp(61)**2) .Or. tmax<0.2) Return
  xe0 = 2.*parp(65)/vint(1)
  alams = paru(111)
  paru(111) = parp(61)
  ns = n
  100 n = ns
  Do jt = 1, 2
    kfls(jt) = mint(14+jt)
    kfls(jt+2) = kfls(jt)
    xs(jt) = vint(40+jt)
    zs(jt) = 1.
    q2s(jt) = parp(67)*q2e
    tevs(jt) = tmax
    alam(jt) = parp(61)
    the2(jt) = 100.
    Do kfl = -6, 6
      xfs(jt, kfl) = xsfx(jt, kfl)
    End Do
  End Do
  dsh = dble(vint(44))
  If (iset(isub)==3 .Or. iset(isub)==4) dsh = dble(vint(26)*vint(2))
  120 n = n + 1
  jt = 1
  If (n>ns+1 .And. q2s(2)>q2s(1)) jt = 2
  kflb = kfls(jt)
  xb = xs(jt)
  Do kfl = -6, 6
    xfb(kfl) = xfs(jt, kfl)
  End Do
  dshr = 2D0*sqrt(dsh)
  dshz = dsh/dble(zs(jt))
  xe = max(xe0, xb*(1./(1.-parp(66))-1.))
  If (xb+xe>=0.999) Then
    q2b = 0.
    Goto 220
  End If
  If (mstp(62)<=1) Then
    q2b = 0.5*(1./zs(jt)+1.)*q2s(jt) + 0.5*(1./zs(jt)-1.)*(q2s(3-jt)-sngl(dsh)+sqrt((sngl(dsh)+q2s(1)+q2s(2))**2+8.*q2s(1)*q2s(2)*zs(jt)/(1.-zs(jt))))
    tevb = log(parp(63)*q2b/alam(jt)**2)
  Else
    q2b = q2s(jt)
    tevb = tevs(jt)
  End If
  alsdum = ulalps(parp(63)*q2b)
  tevb = tevb + 2.*log(alam(jt)/paru(117))
  tevbsv = tevb
  alam(jt) = paru(117)
  b0 = (33.-2.*mstu(118))/6.
  Do kfl = -6, 6
    wtap(kfl) = 0.
    wtsf(kfl) = 0.
  End Do
  If (kflb==21) Then
    wtapq = 16.*(1.-sqrt(xb+xe))/(3.*sqrt(xb))
    Do kfl = -mstp(54), mstp(54)
      If (kfl==0) wtap(kfl) = 6.*log((1.-xb)/xe)
      If (kfl/=0) wtap(kfl) = wtapq
    End Do
  Else
    wtap(0) = 0.5*xb*(1./(xb+xe)-1.)
    wtap(kflb) = 8.*log((1.-xb)*(xb+xe)/xe)/3.
  End If
  160 wtsum = 0.
  If (kflb/=21) xfbo = xfb(kflb)
  If (kflb==21) xfbo = xfb(0)
  If (xfbo==0.0) Then
    Write (mstu(11), 1000)
    Write (mstu(11), 1001) kflb, xfb(kflb)
    xfbo = 0.00001
  End If
  Do kfl = -mstp(54), mstp(54)
    wtsf(kfl) = xfb(kfl)/xfbo
    wtsum = wtsum + wtap(kfl)*wtsf(kfl)
  End Do
  wtsum = max(0.0001, wtsum)
  180 If (mstp(64)<=0) Then
    tevb = tevb + log(rlu(0))*paru(2)/(paru(111)*wtsum)
  Else If (mstp(64)==1) Then
    tevb = tevb*exp(max(-100.,log(rlu(0))*b0/wtsum))
  Else
    tevb = tevb*exp(max(-100.,log(rlu(0))*b0/(5.*wtsum)))
  End If
  190 q2ref = alam(jt)**2*exp(tevb)
  q2b = q2ref/parp(63)
  If (q2b<parp(62)**2) Then
    q2b = 0.
  Else
    wtran = rlu(0)*wtsum
    kfla = -mstp(54) - 1
    200 kfla = kfla + 1
    wtran = wtran - wtap(kfla)*wtsf(kfla)
    If (kfla<mstp(54) .And. wtran>0.) Goto 200
    If (kfla==0) kfla = 21
    If (kflb==21 .And. kfla==21) Then
      z = 1./(1.+((1.-xb)/xb)*(xe/(1.-xb))**rlu(0))
      wtz = (1.-z*(1.-z))**2
    Else If (kflb==21) Then
      z = xb/(1.-rlu(0)*(1.-sqrt(xb+xe)))**2
      wtz = 0.5*(1.+(1.-z)**2)*sqrt(z)
    Else If (kfla==21) Then
      z = xb*(1.+rlu(0)*(1./(xb+xe)-1.))
      wtz = 1. - 2.*z*(1.-z)
    Else
      z = 1. - (1.-xb)*(xe/((xb+xe)*(1.-xb)))**rlu(0)
      wtz = 0.5*(1.+z**2)
    End If
    If (mstp(65)>=1) Then
      rsoft = 6.
      If (kflb/=21) rsoft = 8./3.
      z = z*(tevb/tevs(jt))**(rsoft*xe/((xb+xe)*b0))
      If (z<=xb) Goto 180
    End If
    If (mstp(64)>=2) Then
      If ((1.-z)*q2b<parp(62)**2) Goto 180
      alprat = tevb/(tevb+log(1.-z))
      If (alprat<5.*rlu(0)) Goto 180
      If (alprat>5.) wtz = wtz*alprat/5.
    End If
    If (mstp(62)>=3) Then
      the2t = (4.*z**2*q2b)/(vint(2)*(1.-z)*xb**2)
      If (the2t>the2(jt)) Goto 180
    End If
    Call pystfu(mint(10+jt), xb, q2ref, xfn, jt)
    If (kflb/=21) xfbn = xfn(kflb)
    If (kflb==21) xfbn = xfn(0)
    If (xfbn<1E-20) Then
      If (kfla==kflb) Then
        tevb = tevbsv
        wtap(kflb) = 0.
        Goto 160
      Else If (tevbsv-tevb>0.2) Then
        tevb = 0.5*(tevbsv+tevb)
        Goto 190
      Else
        xfbn = 1E-10
      End If
    End If
    Do kfl = -mstp(54), mstp(54)
      xfb(kfl) = xfn(kfl)
    End Do
    xa = xb/z
    Call pystfu(mint(10+jt), xa, q2ref, xfa, jt)
    If (kfla/=21) xfan = xfa(kfla)
    If (kfla==21) xfan = xfa(0)
    If (xfan<1E-20) Goto 160
    If (kfla/=21) wtsfa = wtsf(kfla)
    If (kfla==21) wtsfa = wtsf(0)
    If (wtz*xfan/xfbn<rlu(0)*wtsfa) Goto 160
  End If
  220 If (n==ns+2) Then
    dq2(jt) = dble(q2b)
    dplcm = dsqrt((dsh+dq2(1)+dq2(2))**2-4D0*dq2(1)*dq2(2))/dshr
    Do jr = 1, 2
      i = ns + jr
      If (jr==1) ipo = ipus1
      If (jr==2) ipo = ipus2
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
      k(i, 1) = 14
      k(i, 2) = kfls(jr+2)
      k(i, 4) = ipo
      k(i, 5) = ipo
      p(i, 3) = sngl(dplcm)*(-1)**(jr+1)
      p(i, 4) = sngl((dsh+dq2(3-jr)-dq2(jr))/dshr)
      p(i, 5) = -sqrt(sngl(dq2(jr)))
      k(ipo, 1) = 14
      k(ipo, 3) = i
      k(ipo, 4) = mod(k(ipo,4), mstu(5)) + mstu(5)*i
      k(ipo, 5) = mod(k(ipo,5), mstu(5)) + mstu(5)*i
    End Do
  Else If (n>ns+2) Then
    jr = 3 - jt
    dq2(3) = dble(q2b)
    dpc(1) = dble(p(is(1),4))
    dpc(2) = dble(p(is(2),4))
    dpc(3) = dble(0.5*(abs(p(is(1),3))+abs(p(is(2),3))))
    dpd(1) = dsh + dq2(jr) + dq2(jt)
    dpd(2) = dshz + dq2(jr) + dq2(3)
    dpd(3) = sqrt(dpd(1)**2-4D0*dq2(jr)*dq2(jt))
    dpd(4) = sqrt(dpd(2)**2-4D0*dq2(jr)*dq2(3))
    ikin = 0
    If (q2s(jr)>=(0.5*parp(62))**2 .And. dpd(1)-dpd(3)>=1D-10*dpd(1)) ikin = 1
    If (ikin==0) dmsma = (dq2(jt)/dble(zs(jt))-dq2(3))*(dsh/(dsh+dq2(jt))-dsh/(dshz+dq2(3)))
    If (ikin==1) dmsma = (dpd(1)*dpd(2)-dpd(3)*dpd(4))/(2.D0*dq2(jr)) - dq2(jt) - dq2(3)
    it = n
    Do j = 1, 5
      k(it, j) = 0
      p(it, j) = 0.
      v(it, j) = 0.
    End Do
    k(it, 1) = 3
    k(it, 2) = 21
    If (kflb==21 .And. kfls(jt+2)/=21) k(it, 2) = -kfls(jt+2)
    If (kflb/=21 .And. kfls(jt+2)==21) k(it, 2) = kflb
    p(it, 5) = ulmass(k(it,2))
    If (sngl(dmsma)<=p(it,5)**2) Goto 100
    If (mstp(63)>=1) Then
      p(it, 4) = sngl((dshz-dsh-dble(p(it,5))**2)/dshr)
      p(it, 3) = sqrt(p(it,4)**2-p(it,5)**2)
      If (mstp(63)==1) Then
        q2tim = sngl(dmsma)
      Else If (mstp(63)==2) Then
        q2tim = min(sngl(dmsma), parp(71)*q2s(jt))
      Else
        q2tim = sngl(dmsma)
      End If
      Call lushow(it, 0, sqrt(q2tim))
      If (n>=it+1) p(it, 5) = p(it+1, 5)
    End If
    dms = dble(p(it,5)**2)
    If (ikin==0) dpt2 = (dmsma-dms)*(dshz+dq2(3))/(dsh+dq2(jt))
    If (ikin==1) dpt2 = (dmsma-dms)*(0.5D0*dpd(1)*dpd(2)+0.5D0*dpd(3)*dpd(4)-dq2(jr)*(dq2(jt)+dq2(3)+dms))/(4.D0*dsh*dpc(3)**2)
    If (dpt2<0.D0) Goto 100
    dpb(1) = (0.5D0*dpd(2)-dpc(jr)*(dshz+dq2(jr)-dq2(jt)-dms)/dshr)/dpc(3) - dpc(3)
    p(it, 1) = sqrt(sngl(dpt2))
    p(it, 3) = sngl(dpb(1))*(-1)**(jt+1)
    p(it, 4) = sngl((dshz-dsh-dms)/dshr)
    If (n>=it+1) Then
      dpb(1) = sqrt(dpb(1)**2+dpt2)
      dpb(2) = sqrt(dpb(1)**2+dms)
      dpb(3) = dble(p(it+1,3))
      dpb(4) = sqrt(dpb(3)**2+dms)
      dbez = (dpb(4)*dpb(1)-dpb(3)*dpb(2))/(dpb(4)*dpb(2)-dpb(3)*dpb(1))
      Call ludbrb(it+1, n, 0., 0., 0D0, 0D0, dbez)
      the = ulangl(p(it,3), p(it,1))
      Call ludbrb(it+1, n, the, 0., 0D0, 0D0, 0D0)
    End If
    Do j = 1, 5
      k(n+1, j) = 0
      p(n+1, j) = 0.
      v(n+1, j) = 0.
    End Do
    k(n+1, 1) = 14
    k(n+1, 2) = kflb
    p(n+1, 1) = p(it, 1)
    p(n+1, 3) = p(it, 3) + p(is(jt), 3)
    p(n+1, 4) = p(it, 4) + p(is(jt), 4)
    p(n+1, 5) = -sqrt(sngl(dq2(3)))
    k(is(jt), 3) = n + 1
    k(it, 3) = n + 1
    id1 = it
    If ((k(n+1,2)>0 .And. k(n+1,2)/=21 .And. k(id1,2)>0 .And. k(id1,2)/=21) .Or. (k(n+1,2)<0 .And. k(id1,2)==21) .Or. (k(n+1,2)==21 .And. k(id1,2)==21 .And. rlu(0)>0.5) .Or. (k(n+1,2)==21 .And. k(id1,2)<0)) id1 = is(jt)
    id2 = it + is(jt) - id1
    k(n+1, 4) = k(n+1, 4) + id1
    k(n+1, 5) = k(n+1, 5) + id2
    k(id1, 4) = k(id1, 4) + mstu(5)*(n+1)
    k(id1, 5) = k(id1, 5) + mstu(5)*id2
    k(id2, 4) = k(id2, 4) + mstu(5)*id1
    k(id2, 5) = k(id2, 5) + mstu(5)*(n+1)
    n = n + 1
    Call ludbrb(ns+1, n, 0., 0., -dble((p(n,1)+p(is(jr),1))/(p(n,4)+p(is(jr),4))), 0D0, -dble((p(n,3)+p(is(jr),3))/(p(n,4)+p(is(jr),4))))
    ir = n + (jt-1)*(is(1)-n)
    Call ludbrb(ns+1, n, -ulangl(p(ir,3),p(ir,1)), paru(2)*rlu(0), 0D0, 0D0, 0D0)
  End If
  is(jt) = n
  q2s(jt) = q2b
  dq2(jt) = dble(q2b)
  If (mstp(62)>=3) the2(jt) = the2t
  dsh = dshz
  If (q2b>=(0.5*parp(62))**2) Then
    kfls(jt+2) = kfls(jt)
    kfls(jt) = kfla
    xs(jt) = xa
    zs(jt) = z
    Do kfl = -6, 6
      xfs(jt, kfl) = xfa(kfl)
    End Do
    tevs(jt) = tevb
  Else
    If (jt==1) ipu1 = n
    If (jt==2) ipu2 = n
  End If
  If (n>mstu(4)-mstu(32)-10) Then
    Call luerrm(11, '(PYSSPA:) no more memory left in LUJETS')
    If (mstu(21)>=1) n = ns
    If (mstu(21)>=1) Return
  End If
  If (max(q2s(1),q2s(2))>=(0.5*parp(62))**2 .Or. n<=ns+1) Goto 120
  Do j = 1, 3
    robo(j+2) = (p(ns+1,j)+p(ns+2,j))/(p(ns+1,4)+p(ns+2,4))
  End Do
  Do j = 1, 5
    p(n+2, j) = p(ns+1, j)
  End Do
  robot = robo(3)**2 + robo(4)**2 + robo(5)**2
  If (robot>=0.999999) Then
    robot = 1.00001*sqrt(robot)
    robo(3) = robo(3)/robot
    robo(4) = robo(4)/robot
    robo(5) = robo(5)/robot
  End If
  Call ludbrb(n+2, n+2, 0., 0., -dble(robo(3)), -dble(robo(4)), -dble(robo(5)))
  robo(2) = ulangl(p(n+2,1), p(n+2,2))
  robo(1) = ulangl(p(n+2,3), sqrt(p(n+2,1)**2+p(n+2,2)**2))
  Call ludbrb(mint(83)+5, ns, robo(1), robo(2), dble(robo(3)), dble(robo(4)), dble(robo(5)))
  k(ipu1, 3) = mint(83) + 3
  k(ipu2, 3) = mint(83) + 4
  Do jt = 1, 2
    mint(12+jt) = kfls(jt)
    vint(140+jt) = xs(jt)
  End Do
  paru(111) = alams
  Return
  1000 Format (5X, 'structure function has a zero point here')
  1001 Format (5X, 'xf(x,i=', I5, ')=', F10.5)
End Subroutine pysspa
