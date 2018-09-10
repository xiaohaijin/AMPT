Subroutine pysigh(nchn, sigs)
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
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension x(2), xpq(-6:6), kfac(2, -40:40), wdtp(0:40), wdte(0:40, 0:5)
  nchn = 0
  sigs = 0.
  isub = mint(1)
  taumin = vint(11)
  ystmin = vint(12)
  ctnmin = vint(13)
  ctpmin = vint(14)
  xt2min = vint(15)
  taupmn = vint(16)
  tau = vint(21)
  yst = vint(22)
  cth = vint(23)
  xt2 = vint(25)
  taup = vint(26)
  taumax = vint(31)
  ystmax = vint(32)
  ctnmax = vint(33)
  ctpmax = vint(34)
  xt2max = vint(35)
  taupmx = vint(36)
  If (iset(isub)<=2 .Or. iset(isub)==5) Then
    x(1) = sqrt(tau)*exp(yst)
    x(2) = sqrt(tau)*exp(-yst)
  Else
    x(1) = sqrt(taup)*exp(yst)
    x(2) = sqrt(taup)*exp(-yst)
  End If
  If (mint(43)==4 .And. iset(isub)>=1 .And. (x(1)>0.999 .Or. x(2)>0.999)) Return
  sh = tau*vint(2)
  sqm3 = vint(63)
  sqm4 = vint(64)
  rm3 = sqm3/sh
  rm4 = sqm4/sh
  be34 = sqrt((1.-rm3-rm4)**2-4.*rm3*rm4)
  rpts = 4.*vint(71)**2/sh
  be34l = sqrt(max(0.,(1.-rm3-rm4)**2-4.*rm3*rm4-rpts))
  rm34 = 2.*rm3*rm4
  rsqm = 1. + rm34
  rthm = (4.*rm3*rm4+rpts)/(1.-rm3-rm4+be34l)
  th = -0.5*sh*max(rthm, 1.-rm3-rm4-be34*cth)
  uh = -0.5*sh*max(rthm, 1.-rm3-rm4+be34*cth)
  sqpth = 0.25*sh*be34**2*(1.-cth**2)
  sh2 = sh**2
  th2 = th**2
  uh2 = uh**2
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    q2 = sh
  Else If (mod(iset(isub),2)==0 .Or. iset(isub)==5) Then
    If (mstp(32)==1) Then
      q2 = 2.*sh*th*uh/(sh**2+th**2+uh**2)
    Else If (mstp(32)==2) Then
      q2 = sqpth + 0.5*(sqm3+sqm4)
    Else If (mstp(32)==3) Then
      q2 = min(-th, -uh)
    Else If (mstp(32)==4) Then
      q2 = sh
    End If
    If (iset(isub)==5 .And. mstp(82)>=2) q2 = q2 + parp(82)**2
  End If
  vint(41) = x(1)
  vint(42) = x(2)
  vint(44) = sh
  vint(43) = sqrt(sh)
  vint(45) = th
  vint(46) = uh
  vint(48) = sqpth
  vint(47) = sqrt(sqpth)
  vint(50) = taup*vint(2)
  vint(49) = sqrt(max(0.,vint(50)))
  vint(52) = q2
  vint(51) = sqrt(q2)
  If (iset(isub)<=0) Goto 145
  If (mint(43)>=2) Then
    q2sf = q2
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      q2sf = pmas(23, 1)**2
      If (isub==8 .Or. isub==76 .Or. isub==77) q2sf = pmas(24, 1)**2
    End If
    Do i = 3 - mint(41), mint(42)
      xsf = x(i)
      If (iset(isub)==5) xsf = x(i)/vint(142+i)
      Call pystfu(mint(10+i), xsf, q2sf, xpq, i)
      Do kfl = -6, 6
        xsfx(i, kfl) = xpq(kfl)
      End Do
    End Do
  End If
  If (mstp(33)/=3) as = ulalps(q2)
  fack = 1.
  faca = 1.
  If (mstp(33)==1) Then
    fack = parp(31)
  Else If (mstp(33)==2) Then
    fack = parp(31)
    faca = parp(32)/parp(31)
  Else If (mstp(33)==3) Then
    q2as = parp(33)*q2
    If (iset(isub)==5 .And. mstp(82)>=2) q2as = q2as + paru(112)*parp(82)
    as = ulalps(q2as)
  End If
  radc = 1. + as/paru(1)
  Do i = 1, 2
    Do j = -40, 40
      kfac(i, j) = 0
    End Do
    If (mint(40+i)==1) Then
      kfac(i, mint(10+i)) = 1
    Else
      Do j = -40, 40
        kfac(i, j) = kfin(i, j)
        If (abs(j)>mstp(54) .And. j/=21) kfac(i, j) = 0
        If (abs(j)<=6) Then
          If (xsfx(i,j)<1.E-10) kfac(i, j) = 0
        Else If (j==21) Then
          If (xsfx(i,0)<1.E-10) kfac(i, 21) = 0
        End If
      End Do
    End If
  End Do
  min1 = 0
  max1 = 0
  min2 = 0
  max2 = 0
  Do j = -20, 20
    If (kfac(1,-j)==1) min1 = -j
    If (kfac(1,j)==1) max1 = j
    If (kfac(2,-j)==1) min2 = -j
    If (kfac(2,j)==1) max2 = j
  End Do
  mina = min(min1, min2)
  maxa = max(max1, max2)
  sqmz = pmas(23, 1)**2
  gmmz = pmas(23, 1)*pmas(23, 2)
  sqmw = pmas(24, 1)**2
  gmmw = pmas(24, 1)*pmas(24, 2)
  sqmh = pmas(25, 1)**2
  gmmh = pmas(25, 1)*pmas(25, 2)
  sqmzp = pmas(32, 1)**2
  gmmzp = pmas(32, 1)*pmas(32, 2)
  sqmhc = pmas(37, 1)**2
  gmmhc = pmas(37, 1)*pmas(37, 2)
  sqmr = pmas(40, 1)**2
  gmmr = pmas(40, 1)*pmas(40, 2)
  aem = paru(101)
  xw = paru(102)
  comfac = paru(1)*paru(5)/vint(2)
  If (mint(43)==4) comfac = comfac*fack
  If ((mint(43)>=2 .Or. iset(isub)==3 .Or. iset(isub)==4) .And. iset(isub)/=5) Then
    atau0 = log(taumax/taumin)
    atau1 = (taumax-taumin)/(taumax*taumin)
    h1 = coef(isub, 1) + (atau0/atau1)*coef(isub, 2)/tau
    If (mint(72)>=1) Then
      taur1 = vint(73)
      gamr1 = vint(74)
      atau2 = log(taumax/taumin*(taumin+taur1)/(taumax+taur1))/taur1
      atau3 = (atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1))/gamr1
      h1 = h1 + (atau0/atau2)*coef(isub, 3)/(tau+taur1) + (atau0/atau3)*coef(isub, 4)*tau/((tau-taur1)**2+gamr1**2)
    End If
    If (mint(72)==2) Then
      taur2 = vint(75)
      gamr2 = vint(76)
      atau4 = log(taumax/taumin*(taumin+taur2)/(taumax+taur2))/taur2
      atau5 = (atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2))/gamr2
      h1 = h1 + (atau0/atau4)*coef(isub, 5)/(tau+taur2) + (atau0/atau5)*coef(isub, 6)*tau/((tau-taur2)**2+gamr2**2)
    End If
    comfac = comfac*atau0/(tau*h1)
  End If
  If (mint(43)==4 .And. iset(isub)/=5) Then
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst2 = ayst1
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))
    h2 = (ayst0/ayst1)*coef(isub, 7)*(yst-ystmin) + (ayst0/ayst2)*coef(isub, 8)*(ystmax-yst) + (ayst0/ayst3)*coef(isub, 9)/cosh(yst)
    comfac = comfac*ayst0/h2
  End If
  acth0 = ctnmax - ctnmin + ctpmax - ctpmin
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    If (mdcy(lucomp(kfpr(isub,1)),1)==1) Then
      If (kfpr(isub,1)==25 .Or. kfpr(isub,1)==37) Then
        comfac = comfac*0.5*acth0
      Else
        comfac = comfac*0.125*(3.*acth0+ctnmax**3-ctnmin**3+ctpmax**3-ctpmin**3)
      End If
    End If
  Else If (iset(isub)==2 .Or. iset(isub)==4) Then
    acth1 = log((max(rm34,rsqm-ctnmin)*max(rm34,rsqm-ctpmin))/(max(rm34,rsqm-ctnmax)*max(rm34,rsqm-ctpmax)))
    acth2 = log((max(rm34,rsqm+ctnmax)*max(rm34,rsqm+ctpmax))/(max(rm34,rsqm+ctnmin)*max(rm34,rsqm+ctpmin)))
    acth3 = 1./max(rm34, rsqm-ctnmax) - 1./max(rm34, rsqm-ctnmin) + 1./max(rm34, rsqm-ctpmax) - 1./max(rm34, rsqm-ctpmin)
    acth4 = 1./max(rm34, rsqm+ctnmin) - 1./max(rm34, rsqm+ctnmax) + 1./max(rm34, rsqm+ctpmin) - 1./max(rm34, rsqm+ctpmax)
    h3 = coef(isub, 10) + (acth0/acth1)*coef(isub, 11)/max(rm34, rsqm-cth) + (acth0/acth2)*coef(isub, 12)/max(rm34, rsqm+cth) + (acth0/acth3)*coef(isub, 13)/max(rm34, rsqm-cth)**2 + (acth0/acth4)*coef(isub, 14)/max(rm34, rsqm+cth)**2
    comfac = comfac*acth0*0.5*be34/h3
  End If
  If (mint(43)>=2 .And. (iset(isub)==3 .Or. iset(isub)==4)) Then
    ataup0 = log(taupmx/taupmn)
    ataup1 = ((1.-tau/taupmx)**4-(1.-tau/taupmn)**4)/(4.*tau)
    h4 = coef(isub, 15) + ataup0/ataup1*coef(isub, 16)/taup*(1.-tau/taup)**3
    If (1.-tau/taup>1.E-4) Then
      fzw = (1.+tau/taup)*log(taup/tau) - 2.*(1.-tau/taup)
    Else
      fzw = 1./6.*(1.-tau/taup)**3*tau/taup
    End If
    comfac = comfac*ataup0*fzw/h4
  End If
  If (iset(isub)==5) Then
    comfac = paru(1)*paru(5)*fack*0.5*vint(2)/sh2
    atau0 = log(2.*(1.+sqrt(1.-xt2))/xt2-1.)
    atau1 = 2.*atan(1./xt2-1.)/sqrt(xt2)
    h1 = coef(isub, 1) + (atau0/atau1)*coef(isub, 2)/sqrt(tau)
    comfac = comfac*atau0/h1
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))
    h2 = (ayst0/ayst1)*coef(isub, 7)*(yst-ystmin) + (ayst0/ayst1)*coef(isub, 8)*(ystmax-yst) + (ayst0/ayst3)*coef(isub, 9)/cosh(yst)
    comfac = comfac*ayst0/h2
    If (mstp(82)<=1) comfac = comfac*xt2**2*(1./vint(149)-1.)
    If (mstp(82)>=2) comfac = comfac*xt2**2/(vint(149)*(1.+vint(149)))
  End If
  145 If (isub<=10) Then
    If (isub==1) Then
      mint(61) = 2
      Call pywidt(23, sqrt(sh), wdtp, wdte)
      facz = comfac*aem**2*4./3.
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 150
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        facf = 1.
        If (iabs(i)<=10) facf = faca/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facf*facz*(ei**2*vint(111)+ei*vi/(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)*vint(112)+(vi**2+ai**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmz)**2+gmmz**2)*vint(114))
      150 End Do
    Else If (isub==2) Then
      Call pywidt(24, sqrt(sh), wdtp, wdte)
      facw = comfac*(aem/xw)**2*1./24*sh2/((sh-sqmw)**2+gmmw**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 170
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 160
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 160
          If ((ia<=10 .And. ja>10) .Or. (ia>10 .And. ja<=10)) Goto 160
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          facf = 1.
          If (ia<=10) facf = vckm((ia+1)/2, (ja+1)/2)*faca/3.
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facf*facw*(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))
        160 End Do
      170 End Do
    Else If (isub==3) Then
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*(aem/xw)**2*1./48.*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 180
        rmq = pmas(iabs(i), 1)**2/sh
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fach*rmq*sqrt(max(0.,1.-4.*rmq))
      180 End Do
    Else If (isub==4) Then
    Else If (isub==5) Then
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*1./(128.*paru(1)**2*16.*(1.-xw)**3)*(aem/xw)**4*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 200
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 190
          ei = kchg(iabs(i), 1)/3.
          ai = sign(1., ei)
          vi = ai - 4.*ei*xw
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*(vi**2+ai**2)*(vj**2+aj**2)
        190 End Do
      200 End Do
    Else If (isub==6) Then
    Else If (isub==7) Then
    Else If (isub==8) Then
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*1./(128*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 220
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 210
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 210
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        210 End Do
      220 End Do
    End If
  Else If (isub<=20) Then
    If (isub==11) Then
      facqq1 = comfac*as**2*4./9.*(sh2+uh2)/th2
      facqqb = comfac*as**2*4./9.*((sh2+uh2)/th2*faca-mstp(34)*2./3.*uh2/(sh*th))
      facqq2 = comfac*as**2*4./9.*((sh2+th2)/uh2-mstp(34)*2./3.*sh2/(th*uh))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 240
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 230
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facqq1
          If (i==-j) sigh(nchn) = facqqb
          If (i==j) Then
            sigh(nchn) = 0.5*sigh(nchn)
            nchn = nchn + 1
            isig(nchn, 1) = i
            isig(nchn, 2) = j
            isig(nchn, 3) = 2
            sigh(nchn) = 0.5*facqq2
          End If
        230 End Do
      240 End Do
    Else If (isub==12) Then
      Call pywidt(21, sqrt(sh), wdtp, wdte)
      facqqb = comfac*as**2*4./9.*(th2+uh2)/sh2*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 250
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facqqb
      250 End Do
    Else If (isub==13) Then
      facgg1 = comfac*as**2*32./27.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)
      facgg2 = comfac*as**2*32./27.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 260
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = 0.5*facgg1
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 2
        sigh(nchn) = 0.5*facgg2
      260 End Do
    Else If (isub==14) Then
      facgg = comfac*as*aem*8./9.*(th2+uh2)/(th*uh)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 270
        ei = kchg(iabs(i), 1)/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgg*ei**2
      270 End Do
    Else If (isub==15) Then
      faczg = comfac*as*aem/(xw*(1.-xw))*1./18.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      faczg = faczg*wids(23, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 280
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczg*(vi**2+ai**2)
      280 End Do
    Else If (isub==16) Then
      facwg = comfac*as*aem/xw*2./9.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 300
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 290
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 290
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facwg*fckm*wids(24, (5-kchw)/2)
        290 End Do
      300 End Do
    Else If (isub==17) Then
    Else If (isub==18) Then
      facgg = comfac*faca*aem**2*1./3.*(th2+uh2)/(th*uh)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 310
        ei = kchg(iabs(i), 1)/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgg*ei**4
      310 End Do
    Else If (isub==19) Then
      facgz = comfac*faca*aem**2/(xw*(1.-xw))*1./24.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      facgz = facgz*wids(23, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 320
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgz*ei**2*(vi**2+ai**2)
      320 End Do
    Else If (isub==20) Then
      facgw = comfac*faca*aem**2/xw*1./6.*((2.*uh-th)/(3.*(sh-sqm4)))**2*(th2+uh2+2.*sqm4*sh)/(th*uh)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 340
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 330
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 330
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facgw*fckm*wids(24, (5-kchw)/2)
        330 End Do
      340 End Do
    End If
  Else If (isub<=30) Then
    If (isub==21) Then
    Else If (isub==22) Then
      faczz = comfac*faca*(aem/(xw*(1.-xw)))**2*1./768.*(uh/th+th/uh+2.*(sqm3+sqm4)*sh/(th*uh)-sqm3*sqm4*(1./th2+1./uh2))
      faczz = faczz*wids(23, 1)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 350
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczz*(vi**4+6.*vi**2*ai**2+ai**4)
      350 End Do
    Else If (isub==23) Then
      faczw = comfac*faca*(aem/xw)**2*1./6.
      faczw = faczw*wids(23, 2)
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 370
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 360
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 360
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          ei = kchg(ia, 1)/3.
          ai = sign(1., ei)
          vi = ai - 4.*ei*xw
          ej = kchg(ja, 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          If (vi+ai>0) Then
            visav = vi
            aisav = ai
            vi = vj
            ai = aj
            vj = visav
            aj = aisav
          End If
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = faczw*fckm*(1./(sh-sqmw)**2*((9.-8.*xw)/4.*thuh+(8.*xw-6.)/4.*sh*(sqm3+sqm4))+(thuh-sh*(sqm3+sqm4))/(2.*(sh-sqmw))*((vj+aj)/th-(vi+ai)/uh)+thuh/(16.*(1.-xw))*((vj+aj)**2/th2+(vi+ai)**2/uh2)+sh*(sqm3+sqm4)/(8.*(1.-xw))*(vi+ai)*(vj+aj)/(th*uh))*wids(24, (5-kchw)/2)
        360 End Do
      370 End Do
    Else If (isub==24) Then
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      fachz = comfac*faca*(aem/(xw*(1.-xw)))**2*1./96.*(thuh+2.*sh*sqmz)/(sh-sqmz)**2
      fachz = fachz*wids(23, 2)*wids(25, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 380
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fachz*(vi**2+ai**2)
      380 End Do
    Else If (isub==25) Then
      facww = comfac*faca*(aem/xw)**2*1./12.
      facww = facww*wids(24, 1)
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 390
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        dsigww = thuh/sh2*(3.-(sh-3.*(sqm3+sqm4))/(sh-sqmz)*(vi+ai)/(2.*ai*(1.-xw))+(sh/(sh-sqmz))**2*(1.-2.*(sqm3+sqm4)/sh+12.*sqm3*sqm4/sh2)*(vi**2+ai**2)/(8.*(1.-xw)**2)) - 2.*sqmz/(sh-sqmz)*(vi+ai)/ai + sqmz*sh/(sh-sqmz)**2*(1.-2.*(sqm3+sqm4)/sh)*(vi**2+ai**2)/(2.*(1.-xw))
        If (kchg(iabs(i),1)<0) Then
          dsigww = dsigww + 2.*(1.+sqmz/(sh-sqmz)*(vi+ai)/(2.*ai))*(thuh/(sh*th)-(sqm3+sqm4)/th) + thuh/th2
        Else
          dsigww = dsigww + 2.*(1.+sqmz/(sh-sqmz)*(vi+ai)/(2.*ai))*(thuh/(sh*uh)-(sqm3+sqm4)/uh) + thuh/uh2
        End If
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facww*dsigww
      390 End Do
    Else If (isub==26) Then
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      fachw = comfac*faca*(aem/xw)**2*1./24.*(thuh+2.*sh*sqmw)/(sh-sqmw)**2
      fachw = fachw*wids(25, 2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 410
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(1,j)==0) Goto 400
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 400
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fachw*fckm*wids(24, (5-kchw)/2)
        400 End Do
      410 End Do
    Else If (isub==27) Then
    Else If (isub==28) Then
      facqg1 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*uh2/th2-uh/sh)*faca
      facqg2 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*sh2/th2-sh/uh)
      Do i = mina, maxa
        If (i==0) Goto 430
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 420
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 420
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facqg1
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 2
          sigh(nchn) = facqg2
        420 End Do
      430 End Do
    Else If (isub==29) Then
      fgq = comfac*faca*as*aem*1./3.*(sh2+uh2)/(-sh*uh)
      Do i = mina, maxa
        If (i==0) Goto 450
        ei = kchg(iabs(i), 1)/3.
        facgq = fgq*ei**2
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 440
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 440
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facgq
        440 End Do
      450 End Do
    Else If (isub==30) Then
      fzq = comfac*faca*as*aem/(xw*(1.-xw))*1./48.*(sh2+uh2+2.*sqm4*th)/(-sh*uh)
      fzq = fzq*wids(23, 2)
      Do i = mina, maxa
        If (i==0) Goto 470
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        faczq = fzq*(vi**2+ai**2)
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 460
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 460
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = faczq
        460 End Do
      470 End Do
    End If
  Else If (isub<=40) Then
    If (isub==31) Then
      facwq = comfac*faca*as*aem/xw*1./12.*(sh2+uh2+2.*sqm4*th)/(-sh*uh)
      Do i = mina, maxa
        If (i==0) Goto 490
        ia = iabs(i)
        kchw = isign(1, kchg(ia,1)*isign(1,i))
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 480
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 480
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facwq*vint(180+i)*wids(24, (5-kchw)/2)
        480 End Do
      490 End Do
    Else If (isub==32) Then
    Else If (isub==33) Then
    Else If (isub==34) Then
    Else If (isub==35) Then
    Else If (isub==36) Then
    Else If (isub==37) Then
    Else If (isub==38) Then
    Else If (isub==39) Then
    Else If (isub==40) Then
    End If
  Else If (isub<=50) Then
    If (isub==41) Then
    Else If (isub==42) Then
    Else If (isub==43) Then
    Else If (isub==44) Then
    Else If (isub==45) Then
    Else If (isub==46) Then
    Else If (isub==47) Then
    Else If (isub==48) Then
    Else If (isub==49) Then
    Else If (isub==50) Then
    End If
  Else If (isub<=60) Then
    If (isub==51) Then
    Else If (isub==52) Then
    Else If (isub==53) Then
      Call pywidt(21, sqrt(sh), wdtp, wdte)
      facqq1 = comfac*as**2*1./6.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facqq2 = comfac*as**2*1./6.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      If (kfac(1,21)*kfac(2,21)==0) Goto 500
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = facqq2
      500 Continue
    Else If (isub==54) Then
    Else If (isub==55) Then
    Else If (isub==56) Then
    Else If (isub==57) Then
    Else If (isub==58) Then
    Else If (isub==59) Then
    Else If (isub==60) Then
    End If
  Else If (isub<=70) Then
    If (isub==61) Then
    Else If (isub==62) Then
    Else If (isub==63) Then
    Else If (isub==64) Then
    Else If (isub==65) Then
    Else If (isub==66) Then
    Else If (isub==67) Then
    Else If (isub==68) Then
      facgg1 = comfac*as**2*9./4.*(sh2/th2+2.*sh/th+3.+2.*th/sh+th2/sh2)*faca
      facgg2 = comfac*as**2*9./4.*(uh2/sh2+2.*uh/sh+3.+2.*sh/uh+sh2/uh2)*faca
      facgg3 = comfac*as**2*9./4.*(th2/uh2+2.*th/uh+3+2.*uh/th+uh2/th2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 510
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = 0.5*facgg1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = 0.5*facgg2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 3
      sigh(nchn) = 0.5*facgg3
      510 Continue
    Else If (isub==69) Then
    Else If (isub==70) Then
    End If
  Else If (isub<=80) Then
    If (isub==71) Then
      be2 = 1. - 4.*sqmz/sh
      th = -0.5*sh*be2*(1.-cth)
      uh = -0.5*sh*be2*(1.+cth)
      shang = 1./(1.-xw)*sqmw/sqmz*(1.+be2)**2
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      thang = 1./(1.-xw)*sqmw/sqmz*(be2-cth)**2
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      uhang = 1./(1.-xw)*sqmw/sqmz*(be2+cth)**2
      auhre = (uh-sqmh)/((uh-sqmh)**2+gmmh**2)*uhang
      auhim = -gmmh/((uh-sqmh)**2+gmmh**2)*uhang
      fach = 0.5*comfac*1./(4096.*paru(1)**2*16.*(1.-xw)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+athre+auhre)**2+(ashim+athim+auhim)**2)*sqmz/sqmw
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 530
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 520
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*avi*avj
        520 End Do
      530 End Do
    Else If (isub==72) Then
      be2 = sqrt((1.-4.*sqmw/sh)*(1.-4.*sqmz/sh))
      cth2 = cth**2
      th = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh-be2*cth)
      uh = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh+be2*cth)
      shang = 4.*sqrt(sqmw/(sqmz*(1.-xw)))*(1.-2.*sqmw/sh)*(1.-2.*sqmz/sh)
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      atwre = (1.-xw)/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3./2.+be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2+2.*(sqmw+sqmz)/sh*be2*cth))
      atwim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3./2.-be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2-2.*(sqmw+sqmz)/sh*be2*cth))
      auwim = 0.
      a4re = 2.*(1.-xw)/sqmz*(3.-cth2-4.*(sqmw+sqmz)/sh)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2*16.*(1.-xw)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+atwre+auwre+a4re)**2+(ashim+atwim+auwim+a4im)**2)*sqmz/sqmw
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 550
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 540
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*avi*avj
        540 End Do
      550 End Do
    Else If (isub==73) Then
      be2 = 1. - 2.*(sqmz+sqmw)/sh + ((sqmz-sqmw)/sh)**2
      ep1 = 1. + (sqmz-sqmw)/sh
      ep2 = 1. - (sqmz-sqmw)/sh
      th = -0.5*sh*be2*(1.-cth)
      uh = (sqmz-sqmw)**2/sh - 0.5*sh*be2*(1.+cth)
      thang = sqrt(sqmw/(sqmz*(1.-xw)))*(be2-ep1*cth)*(be2-ep2*cth)
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      aswre = (1.-xw)/sqmz*sh/(sh-sqmw)*(-be2*(ep1+ep2)**4*cth+1./4.*(be2+ep1*ep2)**2*((ep1-ep2)**2-4.*be2*cth)+2.*be2*(be2+ep1*ep2)*(ep1+ep2)**2*cth-1./16.*sh/sqmw*(ep1**2-ep2**2)**2*(be2+ep1*ep2)**2)
      aswim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*(-be2*(ep2+ep1*cth)*(ep1+ep2*cth)*(be2+ep1*ep2)+be2*(ep2+ep1*cth)*(be2+ep1*ep2*cth)*(2.*ep2-ep2*cth+ep1)-be2*(ep2+ep1*cth)**2*(be2-ep2**2*cth)-1./8.*(be2+ep1*ep2*cth)**2*((ep1+ep2)**2+2.*be2*(1.-cth))+1./32.*sh/sqmw*(be2+ep1*ep2*cth)**2*(ep1**2-ep2**2)**2-be2*(ep1+ep2*cth)*(ep2+ep1*cth)*(be2+ep1*ep2)+be2*(ep1+ep2*cth)*(be2+ep1*ep2*cth)*(2.*ep1-ep1*cth+ep2)-be2*(ep1+ep2*cth)**2*(be2-ep1**2*cth)-1./8.*(be2+ep1*ep2*cth)**2*((ep1+ep2)**2+2.*be2*(1.-cth))+1./32.*sh/sqmw*(be2+ep1*ep2*cth)**2*(ep1**2-ep2**2)**2)
      auwim = 0.
      a4re = (1.-xw)/sqmz*(ep1**2*ep2**2*(cth**2-1.)-2.*be2*(ep1**2+ep2**2+ep1*ep2)*cth-2.*be2*ep1*ep2)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2*4.*(1.-xw))*(aem/xw)**4*(sh/sqmw)**2*((athre+aswre+auwre+a4re)**2+(athim+aswim+auwim+a4im)**2)*sqrt(sqmz/sqmw)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 570
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 560
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = ai - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*(avi*vint(180+j)+vint(180+i)*avj)
        560 End Do
      570 End Do
    Else If (isub==75) Then
    Else If (isub==76) Then
      be2 = sqrt((1.-4.*sqmw/sh)*(1.-4.*sqmz/sh))
      cth2 = cth**2
      th = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh-be2*cth)
      uh = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh+be2*cth)
      shang = 4.*sqrt(sqmw/(sqmz*(1.-xw)))*(1.-2.*sqmw/sh)*(1.-2.*sqmz/sh)
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      atwre = (1.-xw)/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3./2.+be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2+2.*(sqmw+sqmz)/sh*be2*cth))
      atwim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3./2.-be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2-2.*(sqmw+sqmz)/sh*be2*cth))
      auwim = 0.
      a4re = 2.*(1.-xw)/sqmz*(3.-cth2-4.*(sqmw+sqmz)/sh)
      a4im = 0.
      fach = 0.5*comfac*1./(4096.*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+atwre+auwre+a4re)**2+(ashim+atwim+auwim+a4im)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 590
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 580
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 580
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        580 End Do
      590 End Do
    Else If (isub==77) Then
      be2 = 1. - 4.*sqmw/sh
      be4 = be2**2
      cth2 = cth**2
      cth3 = cth**3
      th = -0.5*sh*be2*(1.-cth)
      uh = -0.5*sh*be2*(1.+cth)
      shang = (1.+be2)**2
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      thang = (be2-cth)**2
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      sgzang = 1./sqmw*be2*(3.-be2)**2*cth
      asgre = xw*sgzang
      asgim = 0.
      aszre = (1.-xw)*sh/(sh-sqmz)*sgzang
      aszim = 0.
      tgzang = 1./sqmw*(be2*(4.-2.*be2+be4)+be2*(4.-10.*be2+be4)*cth+(2.-11.*be2+10.*be4)*cth2+be2*cth3)
      atgre = 0.5*xw*sh/th*tgzang
      atgim = 0.
      atzre = 0.5*(1.-xw)*sh/(th-sqmz)*tgzang
      atzim = 0.
      a4re = 1./sqmw*(1.+2.*be2-6.*be2*cth-cth2)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+athre+asgre+aszre+atgre+atzre+a4re)**2+(ashim+athim+asgim+aszim+atgim+atzim+a4im)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 610
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 600
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 600
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        600 End Do
      610 End Do
    Else If (isub==78) Then
    Else If (isub==79) Then
    End If
  Else If (isub<=90) Then
    If (isub==81) Then
      facqqb = comfac*as**2*4./9.*(((th-sqm3)**2+(uh-sqm3)**2)/sh2+2.*sqm3/sh)
      If (mstp(35)>=1) Then
        If (mstp(35)==1) Then
          alssg = parp(35)
        Else
          mst115 = mstu(115)
          mstu(115) = mstp(36)
          q2bn = sqrt(sqm3*((sqrt(sh)-2.*sqrt(sqm3))**2+parp(36)**2))
          alssg = ulalps(q2bn)
          mstu(115) = mst115
        End If
        xrepu = paru(1)*alssg/(6.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        frepu = xrepu/(exp(min(100.,xrepu))-1.)
        pari(81) = frepu
        facqqb = facqqb*frepu
      End If
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 620
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facqqb
      620 End Do
    Else If (isub==82) Then
      facqq1 = comfac*faca*as**2*1./6.*((uh-sqm3)/(th-sqm3)-2.*(uh-sqm3)**2/sh2+4.*sqm3/sh*(th*uh-sqm3**2)/(th-sqm3)**2)
      facqq2 = comfac*faca*as**2*1./6.*((th-sqm3)/(uh-sqm3)-2.*(th-sqm3)**2/sh2+4.*sqm3/sh*(th*uh-sqm3**2)/(uh-sqm3)**2)
      If (mstp(35)>=1) Then
        If (mstp(35)==1) Then
          alssg = parp(35)
        Else
          mst115 = mstu(115)
          mstu(115) = mstp(36)
          q2bn = sqrt(sqm3*((sqrt(sh)-2.*sqrt(sqm3))**2+parp(36)**2))
          alssg = ulalps(q2bn)
          mstu(115) = mst115
        End If
        xattr = 4.*paru(1)*alssg/(3.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        fattr = xattr/(1.-exp(-min(100.,xattr)))
        xrepu = paru(1)*alssg/(6.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        frepu = xrepu/(exp(min(100.,xrepu))-1.)
        fatre = (2.*fattr+5.*frepu)/7.
        pari(81) = fatre
        facqq1 = facqq1*fatre
        facqq2 = facqq2*fatre
      End If
      If (kfac(1,21)*kfac(2,21)==0) Goto 630
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = facqq2
      630 Continue
    End If
  Else If (isub<=100) Then
    If (isub==91) Then
      sigs = xsec(isub, 1)
    Else If (isub==92) Then
      sigs = xsec(isub, 1)
    Else If (isub==93) Then
      sigs = xsec(isub, 1)
    Else If (isub==94) Then
      sigs = xsec(isub, 1)
    Else If (isub==95) Then
      sigs = xsec(isub, 1)
    Else If (isub==96) Then
      Call pywidt(21, sqrt(sh), wdtp, wdte)
      facqq1 = comfac*as**2*4./9.*(sh2+uh2)/th2
      facqqb = comfac*as**2*4./9.*((sh2+uh2)/th2*faca-mstp(34)*2./3.*uh2/(sh*th))
      facqq2 = comfac*as**2*4./9.*((sh2+th2)/uh2-mstp(34)*2./3.*sh2/(th*uh))
      Do i = -3, 3
        If (i==0) Goto 650
        Do j = -3, 3
          If (j==0) Goto 640
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 111
          sigh(nchn) = facqq1
          If (i==-j) sigh(nchn) = facqqb
          If (i==j) Then
            sigh(nchn) = 0.5*sigh(nchn)
            nchn = nchn + 1
            isig(nchn, 1) = i
            isig(nchn, 2) = j
            isig(nchn, 3) = 112
            sigh(nchn) = 0.5*facqq2
          End If
        640 End Do
      650 End Do
      facqqb = comfac*as**2*4./9.*(th2+uh2)/sh2*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))
      facgg1 = comfac*as**2*32./27.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)
      facgg2 = comfac*as**2*32./27.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)
      Do i = -3, 3
        If (i==0) Goto 660
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 121
        sigh(nchn) = facqqb
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 131
        sigh(nchn) = 0.5*facgg1
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 132
        sigh(nchn) = 0.5*facgg2
      660 End Do
      facqg1 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*uh2/th2-uh/sh)*faca
      facqg2 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*sh2/th2-sh/uh)
      Do i = -3, 3
        If (i==0) Goto 680
        Do isde = 1, 2
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 281
          sigh(nchn) = facqg1
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 282
          sigh(nchn) = facqg2
        End Do
      680 End Do
      facqq1 = comfac*as**2*1./6.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facqq2 = comfac*as**2*1./6.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facgg1 = comfac*as**2*9./4.*(sh2/th2+2.*sh/th+3.+2.*th/sh+th2/sh2)*faca
      facgg2 = comfac*as**2*9./4.*(uh2/sh2+2.*uh/sh+3.+2.*sh/uh+sh2/uh2)*faca
      facgg3 = comfac*as**2*9./4.*(th2/uh2+2.*th/uh+3+2.*uh/th+uh2/th2)
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 531
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 532
      sigh(nchn) = facqq2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 681
      sigh(nchn) = 0.5*facgg1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 682
      sigh(nchn) = 0.5*facgg2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 683
      sigh(nchn) = 0.5*facgg3
    End If
  Else If (isub<=110) Then
    If (isub==101) Then
    Else If (isub==102) Then
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      etare = 0.
      etaim = 0.
      Do i = 1, 2*mstp(1)
        eps = 4.*pmas(i, 1)**2/sh
        If (eps<=1.) Then
          If (eps>1.E-4) Then
            root = sqrt(1.-eps)
            rln = log((1.+root)/(1.-root))
          Else
            rln = log(4./eps-2.)
          End If
          phire = 0.25*(rln**2-paru(1)**2)
          phiim = 0.5*paru(1)*rln
        Else
          phire = -(asin(1./sqrt(eps)))**2
          phiim = 0.
        End If
        etare = etare + 0.5*eps*(1.+(eps-1.)*phire)
        etaim = etaim + 0.5*eps*(eps-1.)*phiim
      End Do
      eta2 = etare**2 + etaim**2
      fach = comfac*faca*(as/paru(1)*aem/xw)**2*1./512.*(sh/sqmw)**2*eta2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      If (kfac(1,21)*kfac(2,21)==0) Goto 700
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = fach
      700 Continue
    End If
  Else If (isub<=120) Then
    If (isub==111) Then
      a5stur = 0.
      a5stui = 0.
      Do i = 1, 2*mstp(1)
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epsh = 4.*sqmq/sqmh
        a5stur = a5stur + sqmq/sqmh*(4.+4.*sh/(th+uh)*(pyw1au(epss,1)-pyw1au(epsh,1))+(1.-4.*sqmq/(th+uh))*(pyw2au(epss,1)-pyw2au(epsh,1)))
        a5stui = a5stui + sqmq/sqmh*(4.*sh/(th+uh)*(pyw1au(epss,2)-pyw1au(epsh,2))+(1.-4.*sqmq/(th+uh))*(pyw2au(epss,2)-pyw2au(epsh,2)))
      End Do
      facgh = comfac*faca/(144.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh/sh*(uh**2+th**2)/(uh+th)**2*(a5stur**2+a5stui**2)
      facgh = facgh*wids(25, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 720
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgh
      720 End Do
    Else If (isub==112) Then
      a5tsur = 0.
      a5tsui = 0.
      Do i = 1, 2*mstp(1)
        sqmq = pmas(i, 1)**2
        epst = 4.*sqmq/th
        epsh = 4.*sqmq/sqmh
        a5tsur = a5tsur + sqmq/sqmh*(4.+4.*th/(sh+uh)*(pyw1au(epst,1)-pyw1au(epsh,1))+(1.-4.*sqmq/(sh+uh))*(pyw2au(epst,1)-pyw2au(epsh,1)))
        a5tsui = a5tsui + sqmq/sqmh*(4.*th/(sh+uh)*(pyw1au(epst,2)-pyw1au(epsh,2))+(1.-4.*sqmq/(sh+uh))*(pyw2au(epst,2)-pyw2au(epsh,2)))
      End Do
      facqh = comfac*faca/(384.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh/(-th)*(uh**2+sh**2)/(uh+sh)**2*(a5tsur**2+a5tsui**2)
      facqh = facqh*wids(25, 2)
      Do i = mina, maxa
        If (i==0) Goto 750
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 740
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 740
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facqh
        740 End Do
      750 End Do
    Else If (isub==113) Then
      a2stur = 0.
      a2stui = 0.
      a2ustr = 0.
      a2usti = 0.
      a2tusr = 0.
      a2tusi = 0.
      a4stur = 0.
      a4stui = 0.
      Do i = 6, 2*mstp(1)
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epst = 4.*sqmq/th
        epsu = 4.*sqmq/uh
        epsh = 4.*sqmq/sqmh
        If (epsh<1.E-6) Goto 760
        bestu = 0.5*(1.+sqrt(1.+epss*th/uh))
        beust = 0.5*(1.+sqrt(1.+epsu*sh/th))
        betus = 0.5*(1.+sqrt(1.+epst*uh/sh))
        beuts = bestu
        betsu = beust
        besut = betus
        w3stur = pyi3au(bestu, epsh, 1) - pyi3au(bestu, epss, 1) - pyi3au(bestu, epsu, 1)
        w3stui = pyi3au(bestu, epsh, 2) - pyi3au(bestu, epss, 2) - pyi3au(bestu, epsu, 2)
        w3sutr = pyi3au(besut, epsh, 1) - pyi3au(besut, epss, 1) - pyi3au(besut, epst, 1)
        w3suti = pyi3au(besut, epsh, 2) - pyi3au(besut, epss, 2) - pyi3au(besut, epst, 2)
        w3tsur = pyi3au(betsu, epsh, 1) - pyi3au(betsu, epst, 1) - pyi3au(betsu, epsu, 1)
        w3tsui = pyi3au(betsu, epsh, 2) - pyi3au(betsu, epst, 2) - pyi3au(betsu, epsu, 2)
        w3tusr = pyi3au(betus, epsh, 1) - pyi3au(betus, epst, 1) - pyi3au(betus, epss, 1)
        w3tusi = pyi3au(betus, epsh, 2) - pyi3au(betus, epst, 2) - pyi3au(betus, epss, 2)
        w3ustr = pyi3au(beust, epsh, 1) - pyi3au(beust, epsu, 1) - pyi3au(beust, epst, 1)
        w3usti = pyi3au(beust, epsh, 2) - pyi3au(beust, epsu, 2) - pyi3au(beust, epst, 2)
        w3utsr = pyi3au(beuts, epsh, 1) - pyi3au(beuts, epsu, 1) - pyi3au(beuts, epss, 1)
        w3utsi = pyi3au(beuts, epsh, 2) - pyi3au(beuts, epsu, 2) - pyi3au(beuts, epss, 2)
        b2stur = sqmq/sqmh**2*(sh*(uh-sh)/(sh+uh)+2.*th*uh*(uh+2.*sh)/(sh+uh)**2*(pyw1au(epst,1)-pyw1au(epsh,1))+(sqmq-sh/4.)*(0.5*pyw2au(epss,1)+0.5*pyw2au(epsh,1)-pyw2au(epst,1)+w3stur)+sh**2*(2.*sqmq/(sh+uh)**2-0.5/(sh+uh))*(pyw2au(epst,1)-pyw2au(epsh,1))+0.5*th*uh/sh*(pyw2au(epsh,1)-2.*pyw2au(epst,1))+0.125*(sh-12.*sqmq-4.*th*uh/sh)*w3tsur)
        b2stui = sqmq/sqmh**2*(2.*th*uh*(uh+2.*sh)/(sh+uh)**2*(pyw1au(epst,2)-pyw1au(epsh,2))+(sqmq-sh/4.)*(0.5*pyw2au(epss,2)+0.5*pyw2au(epsh,2)-pyw2au(epst,2)+w3stui)+sh**2*(2.*sqmq/(sh+uh)**2-0.5/(sh+uh))*(pyw2au(epst,2)-pyw2au(epsh,2))+0.5*th*uh/sh*(pyw2au(epsh,2)-2.*pyw2au(epst,2))+0.125*(sh-12.*sqmq-4.*th*uh/sh)*w3tsui)
        b2sutr = sqmq/sqmh**2*(sh*(th-sh)/(sh+th)+2.*uh*th*(th+2.*sh)/(sh+th)**2*(pyw1au(epsu,1)-pyw1au(epsh,1))+(sqmq-sh/4.)*(0.5*pyw2au(epss,1)+0.5*pyw2au(epsh,1)-pyw2au(epsu,1)+w3sutr)+sh**2*(2.*sqmq/(sh+th)**2-0.5/(sh+th))*(pyw2au(epsu,1)-pyw2au(epsh,1))+0.5*uh*th/sh*(pyw2au(epsh,1)-2.*pyw2au(epsu,1))+0.125*(sh-12.*sqmq-4.*uh*th/sh)*w3ustr)
        b2suti = sqmq/sqmh**2*(2.*uh*th*(th+2.*sh)/(sh+th)**2*(pyw1au(epsu,2)-pyw1au(epsh,2))+(sqmq-sh/4.)*(0.5*pyw2au(epss,2)+0.5*pyw2au(epsh,2)-pyw2au(epsu,2)+w3suti)+sh**2*(2.*sqmq/(sh+th)**2-0.5/(sh+th))*(pyw2au(epsu,2)-pyw2au(epsh,2))+0.5*uh*th/sh*(pyw2au(epsh,2)-2.*pyw2au(epsu,2))+0.125*(sh-12.*sqmq-4.*uh*th/sh)*w3usti)
        b2tsur = sqmq/sqmh**2*(th*(uh-th)/(th+uh)+2.*sh*uh*(uh+2.*th)/(th+uh)**2*(pyw1au(epss,1)-pyw1au(epsh,1))+(sqmq-th/4.)*(0.5*pyw2au(epst,1)+0.5*pyw2au(epsh,1)-pyw2au(epss,1)+w3tsur)+th**2*(2.*sqmq/(th+uh)**2-0.5/(th+uh))*(pyw2au(epss,1)-pyw2au(epsh,1))+0.5*sh*uh/th*(pyw2au(epsh,1)-2.*pyw2au(epss,1))+0.125*(th-12.*sqmq-4.*sh*uh/th)*w3stur)
        b2tsui = sqmq/sqmh**2*(2.*sh*uh*(uh+2.*th)/(th+uh)**2*(pyw1au(epss,2)-pyw1au(epsh,2))+(sqmq-th/4.)*(0.5*pyw2au(epst,2)+0.5*pyw2au(epsh,2)-pyw2au(epss,2)+w3tsui)+th**2*(2.*sqmq/(th+uh)**2-0.5/(th+uh))*(pyw2au(epss,2)-pyw2au(epsh,2))+0.5*sh*uh/th*(pyw2au(epsh,2)-2.*pyw2au(epss,2))+0.125*(th-12.*sqmq-4.*sh*uh/th)*w3stui)
        b2tusr = sqmq/sqmh**2*(th*(sh-th)/(th+sh)+2.*uh*sh*(sh+2.*th)/(th+sh)**2*(pyw1au(epsu,1)-pyw1au(epsh,1))+(sqmq-th/4.)*(0.5*pyw2au(epst,1)+0.5*pyw2au(epsh,1)-pyw2au(epsu,1)+w3tusr)+th**2*(2.*sqmq/(th+sh)**2-0.5/(th+sh))*(pyw2au(epsu,1)-pyw2au(epsh,1))+0.5*uh*sh/th*(pyw2au(epsh,1)-2.*pyw2au(epsu,1))+0.125*(th-12.*sqmq-4.*uh*sh/th)*w3utsr)
        b2tusi = sqmq/sqmh**2*(2.*uh*sh*(sh+2.*th)/(th+sh)**2*(pyw1au(epsu,2)-pyw1au(epsh,2))+(sqmq-th/4.)*(0.5*pyw2au(epst,2)+0.5*pyw2au(epsh,2)-pyw2au(epsu,2)+w3tusi)+th**2*(2.*sqmq/(th+sh)**2-0.5/(th+sh))*(pyw2au(epsu,2)-pyw2au(epsh,2))+0.5*uh*sh/th*(pyw2au(epsh,2)-2.*pyw2au(epsu,2))+0.125*(th-12.*sqmq-4.*uh*sh/th)*w3utsi)
        b2ustr = sqmq/sqmh**2*(uh*(th-uh)/(uh+th)+2.*sh*th*(th+2.*uh)/(uh+th)**2*(pyw1au(epss,1)-pyw1au(epsh,1))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,1)+0.5*pyw2au(epsh,1)-pyw2au(epss,1)+w3ustr)+uh**2*(2.*sqmq/(uh+th)**2-0.5/(uh+th))*(pyw2au(epss,1)-pyw2au(epsh,1))+0.5*sh*th/uh*(pyw2au(epsh,1)-2.*pyw2au(epss,1))+0.125*(uh-12.*sqmq-4.*sh*th/uh)*w3sutr)
        b2usti = sqmq/sqmh**2*(2.*sh*th*(th+2.*uh)/(uh+th)**2*(pyw1au(epss,2)-pyw1au(epsh,2))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,2)+0.5*pyw2au(epsh,2)-pyw2au(epss,2)+w3usti)+uh**2*(2.*sqmq/(uh+th)**2-0.5/(uh+th))*(pyw2au(epss,2)-pyw2au(epsh,2))+0.5*sh*th/uh*(pyw2au(epsh,2)-2.*pyw2au(epss,2))+0.125*(uh-12.*sqmq-4.*sh*th/uh)*w3suti)
        b2utsr = sqmq/sqmh**2*(uh*(sh-uh)/(uh+sh)+2.*th*sh*(sh+2.*uh)/(uh+sh)**2*(pyw1au(epst,1)-pyw1au(epsh,1))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,1)+0.5*pyw2au(epsh,1)-pyw2au(epst,1)+w3utsr)+uh**2*(2.*sqmq/(uh+sh)**2-0.5/(uh+sh))*(pyw2au(epst,1)-pyw2au(epsh,1))+0.5*th*sh/uh*(pyw2au(epsh,1)-2.*pyw2au(epst,1))+0.125*(uh-12.*sqmq-4.*th*sh/uh)*w3tusr)
        b2utsi = sqmq/sqmh**2*(2.*th*sh*(sh+2.*uh)/(uh+sh)**2*(pyw1au(epst,2)-pyw1au(epsh,2))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,2)+0.5*pyw2au(epsh,2)-pyw2au(epst,2)+w3utsi)+uh**2*(2.*sqmq/(uh+sh)**2-0.5/(uh+sh))*(pyw2au(epst,2)-pyw2au(epsh,2))+0.5*th*sh/uh*(pyw2au(epsh,2)-2.*pyw2au(epst,2))+0.125*(uh-12.*sqmq-4.*th*sh/uh)*w3tusi)
        b4stur = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epss,1)-pyw2au(epsh,1)+w3stur))
        b4stui = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epss,2)-pyw2au(epsh,2)+w3stui)
        b4tusr = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epst,1)-pyw2au(epsh,1)+w3tusr))
        b4tusi = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epst,2)-pyw2au(epsh,2)+w3tusi)
        b4ustr = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epsu,1)-pyw2au(epsh,1)+w3ustr))
        b4usti = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epsu,2)-pyw2au(epsh,2)+w3usti)
        a2stur = a2stur + b2stur + b2sutr
        a2stui = a2stui + b2stui + b2suti
        a2ustr = a2ustr + b2ustr + b2utsr
        a2usti = a2usti + b2usti + b2utsi
        a2tusr = a2tusr + b2tusr + b2tsur
        a2tusi = a2tusi + b2tusi + b2tsui
        a4stur = a4stur + b4stur + b4ustr + b4tusr
        a4stui = a4stui + b4stui + b4usti + b4tusi
      760 End Do
      facgh = comfac*faca*3./(128.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh**3/(sh*th*uh)*(a2stur**2+a2stui**2+a2ustr**2+a2usti**2+a2tusr**2+a2tusi**2+a4stur**2+a4stui**2)
      facgh = facgh*wids(25, 2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 770
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facgh
      770 Continue
    Else If (isub==114) Then
      asre = 0.
      asim = 0.
      Do i = 1, 2*mstp(1)
        ei = kchg(iabs(i), 1)/3.
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epst = 4.*sqmq/th
        epsu = 4.*sqmq/uh
        If (epss+abs(epst)+abs(epsu)<3.E-6) Then
          a0stur = 1. + (th-uh)/sh*log(th/uh) + 0.5*(th2+uh2)/sh2*(log(th/uh)**2+paru(1)**2)
          a0stui = 0.
          a0tsur = 1. + (sh-uh)/th*log(-sh/uh) + 0.5*(sh2+uh2)/th2*log(-sh/uh)**2
          a0tsui = -paru(1)*((sh-uh)/th+(sh2+uh2)/th2*log(-sh/uh))
          a0utsr = 1. + (th-sh)/uh*log(-th/sh) + 0.5*(th2+sh2)/uh2*log(-th/sh)**2
          a0utsi = paru(1)*((th-sh)/uh+(th2+sh2)/uh2*log(-th/sh))
          a1stur = -1.
          a1stui = 0.
          a2stur = -1.
          a2stui = 0.
        Else
          bestu = 0.5*(1.+sqrt(1.+epss*th/uh))
          beust = 0.5*(1.+sqrt(1.+epsu*sh/th))
          betus = 0.5*(1.+sqrt(1.+epst*uh/sh))
          beuts = bestu
          betsu = beust
          besut = betus
          a0stur = 1. + (1.+2.*th/sh)*pyw1au(epst, 1) + (1.+2.*uh/sh)*pyw1au(epsu, 1) + 0.5*((th2+uh2)/sh2-epss)*(pyw2au(epst,1)+pyw2au(epsu,1)) - 0.25*epst*(1.-0.5*epss)*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) - 0.25*epsu*(1.-0.5*epss)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.25*(-2.*(th2+uh2)/sh2+4.*epss+epst+epsu+0.5*epst*epsu)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a0stui = (1.+2.*th/sh)*pyw1au(epst, 2) + (1.+2.*uh/sh)*pyw1au(epsu, 2) + 0.5*((th2+uh2)/sh2-epss)*(pyw2au(epst,2)+pyw2au(epsu,2)) - 0.25*epst*(1.-0.5*epss)*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) - 0.25*epsu*(1.-0.5*epss)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.25*(-2.*(th2+uh2)/sh2+4.*epss+epst+epsu+0.5*epst*epsu)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
          a0tsur = 1. + (1.+2.*sh/th)*pyw1au(epss, 1) + (1.+2.*uh/th)*pyw1au(epsu, 1) + 0.5*((sh2+uh2)/th2-epst)*(pyw2au(epss,1)+pyw2au(epsu,1)) - 0.25*epss*(1.-0.5*epst)*(pyi3au(betus,epst,1)+pyi3au(betus,epss,1)) - 0.25*epsu*(1.-0.5*epst)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1)) + 0.25*(-2.*(sh2+uh2)/th2+4.*epst+epss+epsu+0.5*epss*epsu)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1))
          a0tsui = (1.+2.*sh/th)*pyw1au(epss, 2) + (1.+2.*uh/th)*pyw1au(epsu, 2) + 0.5*((sh2+uh2)/th2-epst)*(pyw2au(epss,2)+pyw2au(epsu,2)) - 0.25*epss*(1.-0.5*epst)*(pyi3au(betus,epst,2)+pyi3au(betus,epss,2)) - 0.25*epsu*(1.-0.5*epst)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2)) + 0.25*(-2.*(sh2+uh2)/th2+4.*epst+epss+epsu+0.5*epss*epsu)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2))
          a0utsr = 1. + (1.+2.*th/uh)*pyw1au(epst, 1) + (1.+2.*sh/uh)*pyw1au(epss, 1) + 0.5*((th2+sh2)/uh2-epsu)*(pyw2au(epst,1)+pyw2au(epss,1)) - 0.25*epst*(1.-0.5*epsu)*(pyi3au(beust,epsu,1)+pyi3au(beust,epst,1)) - 0.25*epss*(1.-0.5*epsu)*(pyi3au(beuts,epsu,1)+pyi3au(beuts,epss,1)) + 0.25*(-2.*(th2+sh2)/uh2+4.*epsu+epst+epss+0.5*epst*epss)*(pyi3au(betus,epst,1)+pyi3au(betus,epss,1))
          a0utsi = (1.+2.*th/uh)*pyw1au(epst, 2) + (1.+2.*sh/uh)*pyw1au(epss, 2) + 0.5*((th2+sh2)/uh2-epsu)*(pyw2au(epst,2)+pyw2au(epss,2)) - 0.25*epst*(1.-0.5*epsu)*(pyi3au(beust,epsu,2)+pyi3au(beust,epst,2)) - 0.25*epss*(1.-0.5*epsu)*(pyi3au(beuts,epsu,2)+pyi3au(beuts,epss,2)) + 0.25*(-2.*(th2+sh2)/uh2+4.*epsu+epst+epss+0.5*epst*epss)*(pyi3au(betus,epst,2)+pyi3au(betus,epss,2))
          a1stur = -1. - 0.25*(epss+epst+epsu)*(pyw2au(epss,1)+pyw2au(epst,1)+pyw2au(epsu,1)) + 0.25*(epsu+0.5*epss*epst)*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) + 0.25*(epst+0.5*epss*epsu)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.25*(epss+0.5*epst*epsu)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a1stui = -0.25*(epss+epst+epsu)*(pyw2au(epss,2)+pyw2au(epst,2)+pyw2au(epsu,2)) + 0.25*(epsu+0.5*epss*epst)*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) + 0.25*(epst+0.5*epss*epsu)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.25*(epss+0.5*epst*epsu)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
          a2stur = -1. + 0.125*epss*epst*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) + 0.125*epss*epsu*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.125*epst*epsu*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a2stui = 0.125*epss*epst*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) + 0.125*epss*epsu*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.125*epst*epsu*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
        End If
        asre = asre + ei**2*(a0stur+a0tsur+a0utsr+4.*a1stur+a2stur)
        asim = asim + ei**2*(a0stui+a0tsui+a0utsi+4.*a1stui+a2stui)
      End Do
      facgg = comfac*faca/(8.*paru(1)**2)*as**2*aem**2*(asre**2+asim**2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 790
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facgg
      790 Continue
    Else If (isub==115) Then
    Else If (isub==116) Then
    Else If (isub==117) Then
    End If
  Else If (isub<=140) Then
    If (isub==121) Then
    End If
  Else If (isub<=160) Then
    If (isub==141) Then
      mint(61) = 2
      Call pywidt(32, sqrt(sh), wdtp, wdte)
      faczp = comfac*aem**2*4./9.
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 800
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        api = sign(1., ei)
        vpi = api - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczp*(ei**2*vint(111)+ei*vi/(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)*vint(112)+ei*vpi/(8.*xw*(1.-xw))*sh*(sh-sqmzp)/((sh-sqmzp)**2+gmmzp**2)*vint(113)+(vi**2+ai**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmz)**2+gmmz**2)*vint(114)+2.*(vi*vpi+ai*api)/(16.*xw*(1.-xw))**2*sh2*((sh-sqmz)*(sh-sqmzp)+gmmz*gmmzp)/(((sh-sqmz)**2+gmmz**2)*((sh-sqmzp)**2+gmmzp**2))*vint(115)+(vpi**2+api**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmzp)**2+gmmzp**2)*vint(116))
      800 End Do
    Else If (isub==142) Then
      Call pywidt(37, sqrt(sh), wdtp, wdte)
      fhc = comfac*(aem/xw)**2*1./48.*(sh/sqmw)**2*sh2/((sh-sqmhc)**2+gmmhc**2)
      Do i = 1, mstp(54)/2
        il = 2*i - 1
        iu = 2*i
        rmql = pmas(il, 1)**2/sh
        rmqu = pmas(iu, 1)**2/sh
        fachc = fhc*((rmql*paru(121)+rmqu/paru(121))*(1.-rmql-rmqu)-4.*rmql*rmqu)/sqrt(max(0.,(1.-rmql-rmqu)**2-4.*rmql*rmqu))
        If (kfac(1,il)*kfac(2,-iu)==0) Goto 810
        kchhc = (kchg(il,1)-kchg(iu,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = il
        isig(nchn, 2) = -iu
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        810 If (kfac(1,-il)*kfac(2,iu)==0) Goto 820
        kchhc = (-kchg(il,1)+kchg(iu,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = -il
        isig(nchn, 2) = iu
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        820 If (kfac(1,iu)*kfac(2,-il)==0) Goto 830
        kchhc = (kchg(iu,1)-kchg(il,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = iu
        isig(nchn, 2) = -il
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        830 If (kfac(1,-iu)*kfac(2,il)==0) Goto 840
        kchhc = (-kchg(iu,1)+kchg(il,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = -iu
        isig(nchn, 2) = il
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
      840 End Do
    Else If (isub==143) Then
      Call pywidt(40, sqrt(sh), wdtp, wdte)
      facr = comfac*(aem/xw)**2*1./9.*sh2/((sh-sqmr)**2+gmmr**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 860
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 850
          ja = iabs(j)
          If (i*j>0 .Or. iabs(ia-ja)/=2) Goto 850
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facr*(wdte(0,1)+wdte(0,(10-(i+j))/4)+wdte(0,4))
        850 End Do
      860 End Do
    End If
  Else
    If (isub==161) Then
      fhcq = comfac*faca*as*aem/xw*1./24
      Do i = 1, mstp(54)
        iu = i + mod(i, 2)
        sqmq = pmas(iu, 1)**2
        fachcq = fhcq/paru(121)*sqmq/sqmw*(sh/(sqmq-uh)+2.*sqmq*(sqmhc-uh)/(sqmq-uh)**2+(sqmq-uh)/sh+2.*sqmq/(sqmq-uh)+2.*(sqmhc-uh)/(sqmq-uh)*(sqmhc-sqmq-sh)/sh)
        If (kfac(1,-i)*kfac(2,21)==0) Goto 870
        kchhc = isign(1, -kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = -i
        isig(nchn, 2) = 21
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        870 If (kfac(1,i)*kfac(2,21)==0) Goto 880
        kchhc = isign(1, kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = 21
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        880 If (kfac(1,21)*kfac(2,-i)==0) Goto 890
        kchhc = isign(1, -kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = 21
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        890 If (kfac(1,21)*kfac(2,i)==0) Goto 900
        kchhc = isign(1, kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = 21
        isig(nchn, 2) = i
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
      900 End Do
    End If
  End If
  If (isub<=90 .Or. isub>=96) Then
    Do ichn = 1, nchn
      If (mint(41)==2) Then
        kfl1 = isig(ichn, 1)
        If (kfl1==21) kfl1 = 0
        sigh(ichn) = sigh(ichn)*xsfx(1, kfl1)
      End If
      If (mint(42)==2) Then
        kfl2 = isig(ichn, 2)
        If (kfl2==21) kfl2 = 0
        sigh(ichn) = sigh(ichn)*xsfx(2, kfl2)
      End If
      sigs = sigs + sigh(ichn)
    End Do
  End If
  Return
End Subroutine pysigh
