Subroutine pyresd
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Dimension iref(10, 6), kdcy(2), kfl1(2), kfl2(2), nsd(2), ilin(6), coup(6, 4), pk(6, 4), pkk(6, 6), cthe(2), phi(2), wdtp(0:40), wdte(0:40, 0:5)
  Complex fgk, ha(6, 6), hc(6, 6)
  fgk(i1, i2, i3, i4, i5, i6) = 4.*ha(i1, i3)*hc(i2, i6)*(ha(i1,i5)*hc(i1,i4)+ha(i3,i5)*hc(i3,i4))
  digk(dt, du) = -4.D0*d34*d56 + dt*(3.D0*dt+4.D0*du) + dt**2*(dt*du/(d34*d56)-2.D0*(1.D0/d34+1.D0/d56)*(dt+du)+2.D0*(d34/d56+d56/d34))
  djgk(dt, du) = 8.D0*(d34+d56)**2 - 8.D0*(d34+d56)*(dt+du) - 6.D0*dt*du - 2.D0*dt*du*(dt*du/(d34*d56)-2.D0*(1.D0/d34+1.D0/d56)*(dt+du)+2.D0*(d34/d56+d56/d34))
  isub = mint(1)
  sh = vint(44)
  iref(1, 5) = 0
  iref(1, 6) = 0
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    iref(1, 1) = mint(84) + 2 + iset(isub)
    iref(1, 2) = 0
    iref(1, 3) = mint(83) + 6 + iset(isub)
    iref(1, 4) = 0
  Else If (iset(isub)==2 .Or. iset(isub)==4) Then
    iref(1, 1) = mint(84) + 1 + iset(isub)
    iref(1, 2) = mint(84) + 2 + iset(isub)
    iref(1, 3) = mint(83) + 5 + iset(isub)
    iref(1, 4) = mint(83) + 6 + iset(isub)
  End If
  np = 1
  ip = 0
  100 ip = ip + 1
  ninh = 0
  jtmax = 2
  If (ip==1 .And. (iset(isub)==1 .Or. iset(isub)==3)) jtmax = 1
  Do jt = 1, jtmax
    kdcy(jt) = 0
    kfl1(jt) = 0
    kfl2(jt) = 0
    nsd(jt) = iref(ip, jt)
    id = iref(ip, jt)
    If (id==0) Goto 140
    kfa = iabs(k(id,2))
    If (kfa<23 .Or. kfa>40) Goto 140
    If (mdcy(kfa,1)/=0) Then
      If (isub==1 .Or. isub==141) mint(61) = 1
      Call pywidt(kfa, p(id,5), wdtp, wdte)
      If (kchg(kfa,3)==0) Then
        ipm = 2
      Else
        ipm = (5+isign(1,k(id,2)))/2
      End If
      If (jtmax==1 .Or. iabs(k(iref(ip,1),2))/=iabs(k(iref(ip,2),2))) Then
        i12 = 4
      Else
        If (jt==1) i12 = int(4.5+rlu(0))
        i12 = 9 - i12
      End If
      rkfl = (wdte(0,1)+wdte(0,ipm)+wdte(0,i12))*rlu(0)
      Do i = 1, mdcy(kfa, 3)
        idc = i + mdcy(kfa, 2) - 1
        kfl1(jt) = kfdp(idc, 1)*isign(1, k(id,2))
        kfl2(jt) = kfdp(idc, 2)*isign(1, k(id,2))
        rkfl = rkfl - (wdte(i,1)+wdte(i,ipm)+wdte(i,i12))
        If (rkfl<=0.) Goto 130
      End Do
      130 Continue
    End If
    If ((kfa==23 .Or. kfa==24) .And. kfl1(jt)==0) ninh = ninh + 1
    If (kfl1(jt)==0) Goto 140
    kdcy(jt) = 2
    If (iabs(kfl1(jt))<=10 .Or. kfl1(jt)==21) kdcy(jt) = 1
    If ((iabs(kfl1(jt))>=23 .And. iabs(kfl1(jt))<=25) .Or. (iabs(kfl1(jt))==37)) kdcy(jt) = 3
    nsd(jt) = n
    pid5 = p(id, 5)
    If (kdcy(jt)==1) Then
      Call lu2ent(-(n+1), kfl1(jt), kfl2(jt), pid5)
    Else
      Call lu2ent(n+1, kfl1(jt), kfl2(jt), pid5)
    End If
    If (jtmax==1) Then
      cthe(jt) = vint(13) + (vint(33)-vint(13)+vint(34)-vint(14))*rlu(0)
      If (cthe(jt)>vint(33)) cthe(jt) = cthe(jt) + vint(14) - vint(33)
      phi(jt) = vint(24)
    Else
      cthe(jt) = 2.*rlu(0) - 1.
      phi(jt) = paru(2)*rlu(0)
    End If
  140 End Do
  If (mint(3)==1 .And. ip==1) Then
    mint(25) = kfl1(1)
    mint(26) = kfl2(1)
  End If
  If (jtmax==1 .And. kdcy(1)==0) Goto 530
  If (jtmax==2 .And. kdcy(1)==0 .And. kdcy(2)==0) Goto 530
  If (mstp(45)<=0 .Or. iref(ip,2)==0 .Or. ninh>=1) Goto 500
  If (k(iref(1,1),2)==25 .And. ip==1) Goto 500
  If (k(iref(1,1),2)==25 .And. kdcy(1)*kdcy(2)==0) Goto 500
  ilin(1) = mint(84) + 1
  If (k(mint(84)+1,2)>0) ilin(1) = mint(84) + 2
  If (k(ilin(1),2)==21) ilin(1) = 2*mint(84) + 3 - ilin(1)
  ilin(2) = 2*mint(84) + 3 - ilin(1)
  imin = 1
  If (iref(ip,5)==25) imin = 3
  imax = 2
  iord = 1
  If (k(iref(ip,1),2)==23) iord = 2
  If (k(iref(ip,1),2)==24 .And. k(iref(ip,2),2)==-24) iord = 2
  If (iabs(k(iref(ip,iord),2))==25) iord = 3 - iord
  If (kdcy(iord)==0) iord = 3 - iord
  Do jt = iord, 3 - iord, 3 - 2*iord
    If (kdcy(jt)==0) Then
      ilin(imax+1) = nsd(jt)
      imax = imax + 1
    Else If (k(nsd(jt)+1,2)>0) Then
      ilin(imax+1) = n + 2*jt - 1
      ilin(imax+2) = n + 2*jt
      imax = imax + 2
      k(n+2*jt-1, 2) = k(nsd(jt)+1, 2)
      k(n+2*jt, 2) = k(nsd(jt)+2, 2)
    Else
      ilin(imax+1) = n + 2*jt
      ilin(imax+2) = n + 2*jt - 1
      imax = imax + 2
      k(n+2*jt-1, 2) = k(nsd(jt)+1, 2)
      k(n+2*jt, 2) = k(nsd(jt)+2, 2)
    End If
  End Do
  xw = paru(102)
  Do i = imin, imax
    Do j = 1, 4
      coup(i, j) = 0.
    End Do
    kfa = iabs(k(ilin(i),2))
    If (kfa>20) Goto 410
    coup(i, 1) = luchge(kfa)/3.
    coup(i, 2) = (-1)**mod(kfa, 2)
    coup(i, 4) = -2.*coup(i, 1)*xw
    coup(i, 3) = coup(i, 2) + coup(i, 4)
  410 End Do
  sqmz = pmas(23, 1)**2
  gzmz = pmas(23, 1)*pmas(23, 2)
  sqmw = pmas(24, 1)**2
  gzmw = pmas(24, 1)*pmas(24, 2)
  sqmzp = pmas(32, 1)**2
  gzmzp = pmas(32, 1)*pmas(32, 2)
  420 Do i = n + 1, n + 4
    k(i, 1) = 1
    Do j = 1, 5
      p(i, j) = 0.
    End Do
  End Do
  Do jt = 1, jtmax
    If (kdcy(jt)==0) Goto 440
    id = iref(ip, jt)
    p(n+2*jt-1, 3) = 0.5*p(id, 5)
    p(n+2*jt-1, 4) = 0.5*p(id, 5)
    p(n+2*jt, 3) = -0.5*p(id, 5)
    p(n+2*jt, 4) = 0.5*p(id, 5)
    cthe(jt) = 2.*rlu(0) - 1.
    phi(jt) = paru(2)*rlu(0)
    Call ludbrb(n+2*jt-1, n+2*jt, acos(cthe(jt)), phi(jt), dble(p(id,1)/p(id,4)), dble(p(id,2)/p(id,4)), dble(p(id,3)/p(id,4)))
  440 End Do
  Do i = 1, imax
    k(n+4+i, 1) = 1
    p(n+4+i, 4) = sqrt(p(ilin(i),1)**2+p(ilin(i),2)**2+p(ilin(i),3)**2+p(ilin(i),5)**2)
    p(n+4+i, 5) = p(ilin(i), 5)
    Do j = 1, 3
      p(n+4+i, j) = p(ilin(i), j)
    End Do
  End Do
  therr = acos(2.*rlu(0)-1.)
  phirr = paru(2)*rlu(0)
  Call ludbrb(n+5, n+4+imax, therr, phirr, 0D0, 0D0, 0D0)
  Do i = 1, imax
    Do j = 1, 4
      pk(i, j) = p(n+4+i, j)
    End Do
  End Do
  If (isub==22 .Or. isub==23 .Or. isub==25) Then
    Do i1 = imin, imax - 1
      Do i2 = i1 + 1, imax
        ha(i1, i2) = sqrt((pk(i1,4)-pk(i1,3))*(pk(i2,4)+pk(i2,3))/(1E-20+pk(i1,1)**2+pk(i1,2)**2))*cmplx(pk(i1,1), pk(i1,2)) - sqrt((pk(i1,4)+pk(i1,3))*(pk(i2,4)-pk(i2,3))/(1E-20+pk(i2,1)**2+pk(i2,2)**2))*cmplx(pk(i2,1), pk(i2,2))
        hc(i1, i2) = conjg(ha(i1,i2))
        If (i1<=2) ha(i1, i2) = cmplx(0., 1.)*ha(i1, i2)
        If (i1<=2) hc(i1, i2) = cmplx(0., 1.)*hc(i1, i2)
        ha(i2, i1) = -ha(i1, i2)
        hc(i2, i1) = -hc(i1, i2)
      End Do
    End Do
  End If
  Do i = 1, 2
    Do j = 1, 4
      pk(i, j) = -pk(i, j)
    End Do
  End Do
  Do i1 = imin, imax - 1
    Do i2 = i1 + 1, imax
      pkk(i1, i2) = 2.*(pk(i1,4)*pk(i2,4)-pk(i1,1)*pk(i2,1)-pk(i1,2)*pk(i2,2)-pk(i1,3)*pk(i2,3))
      pkk(i2, i1) = pkk(i1, i2)
    End Do
  End Do
  If (iref(ip,5)==25) Then
    wt = 16.*pkk(3, 5)*pkk(4, 6)
    If (ip==1) wtmax = sh**2
    If (ip>=2) wtmax = p(iref(ip,6), 5)**4
  Else If (isub==1) Then
    If (kfa/=37) Then
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei+0.1)
      vi = ai - 4.*ei*xw
      ef = kchg(kfa, 1)/3.
      af = sign(1., ef+0.1)
      vf = af - 4.*ef*xw
      gg = 1.
      gz = 1./(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gzmz**2)
      zz = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmz)**2+gzmz**2)
      If (mstp(43)==1) Then
        gz = 0.
        zz = 0.
      Else If (mstp(43)==2) Then
        gg = 0.
        gz = 0.
      End If
      asym = 2.*(ei*ai*gz*ef*af+4.*vi*ai*zz*vf*af)/(ei**2*gg*ef**2+ei*vi*gz*ef*vf+(vi**2+ai**2)*zz*(vf**2+af**2))
      wt = 1. + asym*cthe(jt) + cthe(jt)**2
      wtmax = 2. + abs(asym)
    Else
      wt = 1. - cthe(jt)**2
      wtmax = 1.
    End If
  Else If (isub==2) Then
    wt = (1.+cthe(jt))**2
    wtmax = 4.
  Else If (isub==15 .Or. isub==19) Then
    wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,3)**2+pkk(2,4)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,4)**2+pkk(2,3)**2)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*((pkk(1,3)+pkk(1,4))**2+(pkk(2,3)+pkk(2,4))**2)
  Else If (isub==16 .Or. isub==20) Then
    wt = pkk(1, 3)**2 + pkk(2, 4)**2
    wtmax = (pkk(1,3)+pkk(1,4))**2 + (pkk(2,3)+pkk(2,4))**2
  Else If (isub==22) Then
    s34 = p(iref(ip,iord), 5)**2
    s56 = p(iref(ip,3-iord), 5)**2
    ti = pkk(1, 3) + pkk(1, 4) + s34
    ui = pkk(1, 5) + pkk(1, 6) + s56
    wt = coup(1, 3)**4*((coup(3,3)*coup(5,3)*abs(fgk(1,2,3,4,5,6)/ti+fgk(1,2,5,6,3,4)/ui))**2+(coup(3,4)*coup(5,3)*abs(fgk(1,2,4,3,5,6)/ti+fgk(1,2,5,6,4,3)/ui))**2+(coup(3,3)*coup(5,4)*abs(fgk(1,2,3,4,6,5)/ti+fgk(1,2,6,5,3,4)/ui))**2+(coup(3,4)*coup(5,4)*abs(fgk(1,2,4,3,6,5)/ti+fgk(1,2,6,5,4,3)/ui))**2) + coup(1, 4)**4*((coup(3,3)*coup(5,3)*abs(fgk(2,1,5,6,3,4)/ti+fgk(2,1,3,4,5,6)/ui))**2+(coup(3,4)*coup(5,3)*abs(fgk(2,1,6,5,3,4)/ti+fgk(2,1,3,4,6,5)/ui))**2+(coup(3,3)*coup(5,4)*abs(fgk(2,1,5,6,4,3)/ti+fgk(2,1,4,3,5,6)/ui))**2+(coup(3,4)*coup(5,4)*abs(fgk(2,1,6,5,4,3)/ti+fgk(2,1,4,3,6,5)/ui))**2)
    wtmax = 4.*s34*s56*(coup(1,3)**4+coup(1,4)**4)*(coup(3,3)**2+coup(3,4)**2)*(coup(5,3)**2+coup(5,4)**2)*4.*(ti/ui+ui/ti+2.*sh*(s34+s56)/(ti*ui)-s34*s56*(1./ti**2+1./ui**2))
  Else If (isub==23) Then
    d34 = dble(p(iref(ip,iord),5)**2)
    d56 = dble(p(iref(ip,3-iord),5)**2)
    dt = dble(pkk(1,3)+pkk(1,4)) + d34
    du = dble(pkk(1,5)+pkk(1,6)) + d56
    cawz = coup(2, 3)/sngl(dt) - 2.*(1.-xw)*coup(1, 2)/(sh-sqmw)
    cbwz = coup(1, 3)/sngl(du) + 2.*(1.-xw)*coup(1, 2)/(sh-sqmw)
    wt = coup(5, 3)**2*abs(cawz*fgk(1,2,3,4,5,6)+cbwz*fgk(1,2,5,6,3,4))**2 + coup(5, 4)**2*abs(cawz*fgk(1,2,3,4,6,5)+cbwz*fgk(1,2,6,5,3,4))**2
    wtmax = 4.*sngl(d34*d56)*(coup(5,3)**2+coup(5,4)**2)*(cawz**2*sngl(digk(dt,du))+cbwz**2*sngl(digk(du,dt))+cawz*cbwz*sngl(djgk(dt,du)))
  Else If (isub==24) Then
    wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*pkk(1, 3)*pkk(2, 4) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*pkk(1, 4)*pkk(2, 3)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*(pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))
  Else If (isub==25) Then
    d34 = dble(p(iref(ip,iord),5)**2)
    d56 = dble(p(iref(ip,3-iord),5)**2)
    dt = dble(pkk(1,3)+pkk(1,4)) + d34
    du = dble(pkk(1,5)+pkk(1,6)) + d56
    cdww = (coup(1,3)*sqmz/(sh-sqmz)+coup(1,2))/sh
    caww = cdww + 0.5*(coup(1,2)+1.)/sngl(dt)
    cbww = cdww + 0.5*(coup(1,2)-1.)/sngl(du)
    ccww = coup(1, 4)*sqmz/(sh-sqmz)/sh
    wt = abs(caww*fgk(1,2,3,4,5,6)-cbww*fgk(1,2,5,6,3,4))**2 + ccww**2*abs(fgk(2,1,5,6,3,4)-fgk(2,1,3,4,5,6))**2
    wtmax = 4.*sngl(d34*d56)*(caww**2*sngl(digk(dt,du))+cbww**2*sngl(digk(du,dt))-caww*cbww*sngl(djgk(dt,du))+ccww**2*sngl(digk(dt,du)+digk(du,dt)-djgk(dt,du)))
  Else If (isub==26) Then
    wt = pkk(1, 3)*pkk(2, 4)
    wtmax = (pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))
  Else If (isub==30) Then
    If (k(ilin(1),2)>0) wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,4)**2+pkk(3,5)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,3)**2+pkk(4,5)**2)
    If (k(ilin(1),2)<0) wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,3)**2+pkk(4,5)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,4)**2+pkk(3,5)**2)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*((pkk(1,3)+pkk(1,4))**2+(pkk(3,5)+pkk(4,5))**2)
  Else If (isub==31) Then
    If (k(ilin(1),2)>0) wt = pkk(1, 4)**2 + pkk(3, 5)**2
    If (k(ilin(1),2)<0) wt = pkk(1, 3)**2 + pkk(4, 5)**2
    wtmax = (pkk(1,3)+pkk(1,4))**2 + (pkk(3,5)+pkk(4,5))**2
  Else If (isub==141) Then
    ei = kchg(iabs(mint(15)), 1)/3.
    ai = sign(1., ei+0.1)
    vi = ai - 4.*ei*xw
    api = sign(1., ei+0.1)
    vpi = api - 4.*ei*xw
    ef = kchg(kfa, 1)/3.
    af = sign(1., ef+0.1)
    vf = af - 4.*ef*xw
    apf = sign(1., ef+0.1)
    vpf = apf - 4.*ef*xw
    gg = 1.
    gz = 1./(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gzmz**2)
    gzp = 1./(8.*xw*(1.-xw))*sh*(sh-sqmzp)/((sh-sqmzp)**2+gzmzp**2)
    zz = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmz)**2+gzmz**2)
    zzp = 2./(16.*xw*(1.-xw))**2*sh**2*((sh-sqmz)*(sh-sqmzp)+gzmz*gzmzp)/(((sh-sqmz)**2+gzmz**2)*((sh-sqmzp)**2+gzmzp**2))
    zpzp = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmzp)**2+gzmzp**2)
    If (mstp(44)==1) Then
      gz = 0.
      gzp = 0.
      zz = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==2) Then
      gg = 0.
      gz = 0.
      gzp = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==3) Then
      gg = 0.
      gz = 0.
      gzp = 0.
      zz = 0.
      zzp = 0.
    Else If (mstp(44)==4) Then
      gzp = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==5) Then
      gz = 0.
      zz = 0.
      zzp = 0.
    Else If (mstp(44)==6) Then
      gg = 0.
      gz = 0.
      gzp = 0.
    End If
    asym = 2.*(ei*ai*gz*ef*af+ei*api*gzp*ef*apf+4.*vi*ai*zz*vf*af+(vi*api+vpi*ai)*zzp*(vf*apf+vpf*af)+4.*vpi*api*zpzp*vpf*apf)/(ei**2*gg*ef**2+ei*vi*gz*ef*vf+ei*vpi*gzp*ef*vpf+(vi**2+ai**2)*zz*(vf**2+af**2)+(vi*vpi+ai*api)*zzp*(vf*vpf+af*apf)+(vpi**2+api**2)*zpzp*(vpf**2+apf**2))
    wt = 1. + asym*cthe(jt) + cthe(jt)**2
    wtmax = 2. + abs(asym)
  Else
    wt = 1.
    wtmax = 1.
  End If
  If (wt<rlu(0)*wtmax) Goto 420
  500 Do jt = 1, jtmax
    If (kdcy(jt)==0) Goto 520
    id = iref(ip, jt)
    Call ludbrb(nsd(jt)+1, nsd(jt)+2, acos(cthe(jt)), phi(jt), dble(p(id,1)/p(id,4)), dble(p(id,2)/p(id,4)), dble(p(id,3)/p(id,4)))
    k(id, 1) = k(id, 1) + 10
    k(id, 4) = nsd(jt) + 1
    k(id, 5) = nsd(jt) + 2
    idoc = mint(83) + mint(4)
    Do i = nsd(jt) + 1, nsd(jt) + 2
      mint(4) = mint(4) + 1
      i1 = mint(83) + mint(4)
      k(i, 3) = i1
      k(i1, 1) = 21
      k(i1, 2) = k(i, 2)
      k(i1, 3) = iref(ip, jt+2)
      Do j = 1, 5
        p(i1, j) = p(i, j)
      End Do
    End Do
    If (jtmax==1) Then
      mint(7) = mint(83) + 6 + 2*iset(isub)
      mint(8) = mint(83) + 7 + 2*iset(isub)
    End If
    pid5 = p(id, 5)
    If (mstp(71)>=1 .And. kdcy(jt)==1) Call lushow(nsd(jt)+1, nsd(jt)+2, pid5)
    If (kdcy(jt)/=3) Goto 520
    np = np + 1
    iref(np, 1) = nsd(jt) + 1
    iref(np, 2) = nsd(jt) + 2
    iref(np, 3) = idoc + 1
    iref(np, 4) = idoc + 2
    iref(np, 5) = k(iref(ip,jt), 2)
    iref(np, 6) = iref(ip, jt)
  520 End Do
  530 If (ip<np) Goto 100
  Return
End Subroutine pyresd
