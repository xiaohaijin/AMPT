Subroutine pyscat
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
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension wdtp(0:40), wdte(0:40, 0:5), pmq(2), z(2), cthe(2), phi(2)
  isub = mint(1)
  idoc = 6 + iset(isub)
  If (isub==95) idoc = 8
  mint(3) = idoc - 6
  If (idoc>=9) idoc = idoc + 2
  mint(4) = idoc
  ipu1 = mint(84) + 1
  ipu2 = mint(84) + 2
  ipu3 = mint(84) + 3
  ipu4 = mint(84) + 4
  ipu5 = mint(84) + 5
  ipu6 = mint(84) + 6
  Do jt = 1, mstp(126) + 10
    i = mint(83) + jt
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do
  Do jt = 1, 2
    i = mint(83) + jt
    k(i, 1) = 21
    k(i, 2) = mint(10+jt)
    p(i, 1) = 0.
    p(i, 2) = 0.
    p(i, 5) = vint(2+jt)
    p(i, 3) = vint(5)*(-1)**(jt+1)
    p(i, 4) = sqrt(p(i,3)**2+p(i,5)**2)
  End Do
  mint(6) = 2
  kfres = 0
  sh = vint(44)
  shr = sqrt(sh)
  shp = vint(26)*vint(2)
  shpr = sqrt(shp)
  shuser = shr
  If (iset(isub)>=3) shuser = shpr
  Do jt = 1, 2
    i = mint(84) + jt
    k(i, 1) = 14
    k(i, 2) = mint(14+jt)
    k(i, 3) = mint(83) + 2 + jt
    p(i, 5) = ulmass(k(i,2))
  End Do
  If (p(ipu1,5)+p(ipu2,5)>=shuser) Then
    p(ipu1, 5) = 0.
    p(ipu2, 5) = 0.
  End If
  p(ipu1, 4) = 0.5*(shuser+(p(ipu1,5)**2-p(ipu2,5)**2)/shuser)
  p(ipu1, 3) = sqrt(max(0.,p(ipu1,4)**2-p(ipu1,5)**2))
  p(ipu2, 4) = shuser - p(ipu1, 4)
  p(ipu2, 3) = -p(ipu1, 3)
  Do jt = 1, 2
    i1 = mint(83) + 4 + jt
    i2 = mint(84) + jt
    k(i1, 1) = 21
    k(i1, 2) = k(i2, 2)
    k(i1, 3) = i1 - 2
    Do j = 1, 5
      p(i1, j) = p(i2, j)
    End Do
  End Do
  If (isub==12 .Or. isub==53) Then
    Call pywidt(21, shr, wdtp, wdte)
    rkfl = (wdte(0,1)+wdte(0,2)+wdte(0,4))*rlu(0)
    Do i = 1, 2*mstp(1)
      kflq = i
      rkfl = rkfl - (wdte(i,1)+wdte(i,2)+wdte(i,4))
      If (rkfl<=0.) Goto 150
    End Do
    150 Continue
  End If
  js = 1
  mint(21) = mint(15)
  mint(22) = mint(16)
  mint(23) = 0
  mint(24) = 0
  kcc = 20
  kcs = isign(1, mint(15))
  If (isub<=10) Then
    If (isub==1) Then
      kfres = 23
    Else If (isub==2) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      kfres = isign(24, kch1+kch2)
    Else If (isub==3) Then
      kfres = 25
    Else If (isub==4) Then
    Else If (isub==5) Then
      xh = sh/shp
      mint(21) = mint(15)
      mint(22) = mint(16)
      pmq(1) = ulmass(mint(21))
      pmq(2) = ulmass(mint(22))
      240 jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 240
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 240
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 240
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 240
      kcc = 22
      kfres = 25
    Else If (isub==6) Then
    Else If (isub==7) Then
    Else If (isub==8) Then
      xh = sh/shp
      250 Do jt = 1, 2
        i = mint(14+jt)
        ia = iabs(i)
        If (ia<=10) Then
          rvckm = vint(180+i)*rlu(0)
          Do j = 1, mstp(1)
            ib = 2*j - 1 + mod(ia, 2)
            ipm = (5-isign(1,i))/2
            idc = j + mdcy(ia, 2) + 2
            If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 270
            mint(20+jt) = isign(ib, i)
            rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
            If (rvckm<=0.) Goto 280
          270 End Do
        Else
          ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
          mint(20+jt) = isign(ib, i)
        End If
        280 pmq(jt) = ulmass(mint(20+jt))
      End Do
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 250
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 250
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 250
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 250
      kcc = 22
      kfres = 25
    End If
  Else If (isub<=20) Then
    If (isub==11) Then
      kcc = mint(2)
      If (mint(15)*mint(16)<0) kcc = kcc + 2
    Else If (isub==12) Then
      mint(21) = isign(kflq, mint(15))
      mint(22) = -mint(21)
      kcc = 4
    Else If (isub==13) Then
      mint(21) = 21
      mint(22) = 21
      kcc = mint(2) + 4
    Else If (isub==14) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 22
      kcc = 17 + js
    Else If (isub==15) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 23
      kcc = 17 + js
    Else If (isub==16) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 21
      mint(23-js) = isign(24, kch1+kch2)
      kcc = 17 + js
    Else If (isub==17) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 25
      kcc = 17 + js
    Else If (isub==18) Then
      mint(21) = 22
      mint(22) = 22
    Else If (isub==19) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 22
      mint(23-js) = 23
    Else If (isub==20) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 22
      mint(23-js) = isign(24, kch1+kch2)
    End If
  Else If (isub<=30) Then
    If (isub==21) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 22
      mint(23-js) = 25
    Else If (isub==22) Then
      mint(21) = 23
      mint(22) = 23
    Else If (isub==23) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 23
      mint(23-js) = isign(24, kch1+kch2)
    Else If (isub==24) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 23
      mint(23-js) = 25
    Else If (isub==25) Then
      mint(21) = -isign(24, mint(15))
      mint(22) = -mint(21)
    Else If (isub==26) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)>0) js = 2
      mint(20+js) = isign(24, kch1+kch2)
      mint(23-js) = 25
    Else If (isub==27) Then
    Else If (isub==28) Then
      kcc = mint(2) + 6
      If (mint(15)==21) kcc = kcc + 2
      If (mint(15)/=21) kcs = isign(1, mint(15))
      If (mint(16)/=21) kcs = isign(1, mint(16))
    Else If (isub==29) Then
      If (mint(15)==21) js = 2
      mint(23-js) = 22
      kcc = 15 + js
      kcs = isign(1, mint(14+js))
    Else If (isub==30) Then
      If (mint(15)==21) js = 2
      mint(23-js) = 23
      kcc = 15 + js
      kcs = isign(1, mint(14+js))
    End If
  Else If (isub<=40) Then
    If (isub==31) Then
      If (mint(15)==21) js = 2
      i = mint(14+js)
      ia = iabs(i)
      mint(23-js) = isign(24, kchg(ia,1)*i)
      rvckm = vint(180+i)*rlu(0)
      Do j = 1, mstp(1)
        ib = 2*j - 1 + mod(ia, 2)
        ipm = (5-isign(1,i))/2
        idc = j + mdcy(ia, 2) + 2
        If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 220
        mint(20+js) = isign(ib, i)
        rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
        If (rvckm<=0.) Goto 230
      220 End Do
      230 kcc = 15 + js
      kcs = isign(1, mint(14+js))
    Else If (isub==32) Then
      If (mint(15)==21) js = 2
      mint(23-js) = 25
      kcc = 15 + js
      kcs = isign(1, mint(14+js))
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
      kcs = (-1)**int(1.5+rlu(0))
      mint(21) = isign(kflq, kcs)
      mint(22) = -mint(21)
      kcc = mint(2) + 10
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
      kcc = mint(2) + 12
      kcs = (-1)**int(1.5+rlu(0))
    Else If (isub==69) Then
    Else If (isub==70) Then
    End If
  Else If (isub<=80) Then
    If (isub==71 .Or. isub==72) Then
      xh = sh/shp
      mint(21) = mint(15)
      mint(22) = mint(16)
      pmq(1) = ulmass(mint(21))
      pmq(2) = ulmass(mint(22))
      290 jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 290
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 290
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 290
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 290
      kcc = 22
    Else If (isub==73) Then
      xh = sh/shp
      300 jt = int(1.5+rlu(0))
      i = mint(14+jt)
      ia = iabs(i)
      If (ia<=10) Then
        rvckm = vint(180+i)*rlu(0)
        Do j = 1, mstp(1)
          ib = 2*j - 1 + mod(ia, 2)
          ipm = (5-isign(1,i))/2
          idc = j + mdcy(ia, 2) + 2
          If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 320
          mint(20+jt) = isign(ib, i)
          rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
          If (rvckm<=0.) Goto 330
        320 End Do
      Else
        ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
        mint(20+jt) = isign(ib, i)
      End If
      330 pmq(jt) = ulmass(mint(20+jt))
      mint(23-jt) = mint(17-jt)
      pmq(3-jt) = ulmass(mint(23-jt))
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 300
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 300
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 300
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 300
      kcc = 22
    Else If (isub==74) Then
    Else If (isub==75) Then
    Else If (isub==76 .Or. isub==77) Then
      xh = sh/shp
      340 Do jt = 1, 2
        i = mint(14+jt)
        ia = iabs(i)
        If (ia<=10) Then
          rvckm = vint(180+i)*rlu(0)
          Do j = 1, mstp(1)
            ib = 2*j - 1 + mod(ia, 2)
            ipm = (5-isign(1,i))/2
            idc = j + mdcy(ia, 2) + 2
            If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 360
            mint(20+jt) = isign(ib, i)
            rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
            If (rvckm<=0.) Goto 370
          360 End Do
        Else
          ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
          mint(20+jt) = isign(ib, i)
        End If
        370 pmq(jt) = ulmass(mint(20+jt))
      End Do
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 340
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 340
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 340
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 340
      kcc = 22
    Else If (isub==78) Then
    Else If (isub==79) Then
    End If
  Else If (isub<=90) Then
    If (isub==81) Then
      mint(21) = isign(mint(46), mint(15))
      mint(22) = -mint(21)
      kcc = 4
    Else If (isub==82) Then
      kcs = (-1)**int(1.5+rlu(0))
      mint(21) = isign(mint(46), kcs)
      mint(22) = -mint(21)
      kcc = mint(2) + 10
    End If
  Else If (isub<=100) Then
    If (isub==95) Then
      kcc = mint(2) + 12
      kcs = (-1)**int(1.5+rlu(0))
    Else If (isub==96) Then
    End If
  Else If (isub<=110) Then
    If (isub==101) Then
      kcc = 21
      kfres = 22
    Else If (isub==102) Then
      kcc = 21
      kfres = 25
    End If
  Else If (isub<=120) Then
    If (isub==111) Then
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 25
      kcc = 17 + js
    Else If (isub==112) Then
      If (mint(15)==21) js = 2
      mint(23-js) = 25
      kcc = 15 + js
      kcs = isign(1, mint(14+js))
    Else If (isub==113) Then
      If (rlu(0)>0.5) js = 2
      mint(23-js) = 25
      kcc = 22 + js
      kcs = (-1)**int(1.5+rlu(0))
    Else If (isub==114) Then
      If (rlu(0)>0.5) js = 2
      mint(21) = 22
      mint(22) = 22
      kcc = 21
    Else If (isub==115) Then
    Else If (isub==116) Then
    Else If (isub==117) Then
    End If
  Else If (isub<=140) Then
    If (isub==121) Then
    End If
  Else If (isub<=160) Then
    If (isub==141) Then
      kfres = 32
    Else If (isub==142) Then
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      kfres = isign(37, kch1+kch2)
    Else If (isub==143) Then
      kfres = isign(40, mint(15)+mint(16))
    End If
  Else
    If (isub==161) Then
      If (mint(16)==21) js = 2
      ia = iabs(mint(17-js))
      mint(20+js) = isign(37, kchg(ia,1)*mint(17-js))
      ja = ia + mod(ia, 2) - mod(ia+1, 2)
      mint(23-js) = isign(ja, mint(17-js))
      kcc = 18 - js
      If (mint(15)/=21) kcs = isign(1, mint(15))
      If (mint(16)/=21) kcs = isign(1, mint(16))
    End If
  End If
  If (idoc==7) Then
    i = mint(83) + 7
    k(ipu3, 1) = 1
    k(ipu3, 2) = kfres
    k(ipu3, 3) = i
    p(ipu3, 4) = shuser
    p(ipu3, 5) = shuser
    k(ipu1, 4) = ipu2
    k(ipu1, 5) = ipu2
    k(ipu2, 4) = ipu1
    k(ipu2, 5) = ipu1
    k(i, 1) = 21
    k(i, 2) = kfres
    p(i, 4) = shuser
    p(i, 5) = shuser
    n = ipu3
    mint(21) = kfres
    mint(22) = 0
  Else If (idoc==8) Then
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      If (iabs(k(i,2))<=10 .Or. k(i,2)==21) Then
        p(i, 5) = ulmass(k(i,2))
      Else
        p(i, 5) = sqrt(vint(63+mod(js+jt,2)))
      End If
    End Do
    If (p(ipu3,5)+p(ipu4,5)>=shr) Then
      kfa1 = iabs(mint(21))
      kfa2 = iabs(mint(22))
      If ((kfa1>3 .And. kfa1/=21) .Or. (kfa2>3 .And. kfa2/=21)) Then
        mint(51) = 1
        Return
      End If
      p(ipu3, 5) = 0.
      p(ipu4, 5) = 0.
    End If
    p(ipu3, 4) = 0.5*(shr+(p(ipu3,5)**2-p(ipu4,5)**2)/shr)
    p(ipu3, 3) = sqrt(max(0.,p(ipu3,4)**2-p(ipu3,5)**2))
    p(ipu4, 4) = shr - p(ipu3, 4)
    p(ipu4, 3) = -p(ipu3, 3)
    n = ipu4
    mint(7) = mint(83) + 7
    mint(8) = mint(83) + 8
    Call ludbrb(ipu3, ipu4, acos(vint(23)), vint(24), 0D0, 0D0, 0D0)
  Else If (idoc==9) Then
  Else If (idoc==11) Then
    phi(1) = paru(2)*rlu(0)
    phi(2) = phi(1) - phir
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      p(i, 5) = ulmass(k(i,2))
      If (0.5*shpr*z(jt)<=p(i,5)) p(i, 5) = 0.
      pabs = sqrt(max(0.,(0.5*shpr*z(jt))**2-p(i,5)**2))
      ptabs = pabs*sqrt(max(0.,1.-cthe(jt)**2))
      p(i, 1) = ptabs*cos(phi(jt))
      p(i, 2) = ptabs*sin(phi(jt))
      p(i, 3) = pabs*cthe(jt)*(-1)**(jt+1)
      p(i, 4) = 0.5*shpr*z(jt)
      izw = mint(83) + 6 + jt
      k(izw, 1) = 21
      k(izw, 2) = 23
      If (isub==8) k(izw, 2) = isign(24, luchge(mint(14+jt)))
      k(izw, 3) = izw - 2
      p(izw, 1) = -p(i, 1)
      p(izw, 2) = -p(i, 2)
      p(izw, 3) = (0.5*shpr-pabs*cthe(jt))*(-1)**(jt+1)
      p(izw, 4) = 0.5*shpr*(1.-z(jt))
      p(izw, 5) = -sqrt(max(0.,p(izw,3)**2+ptabs**2-p(izw,4)**2))
    End Do
    i = mint(83) + 9
    k(ipu5, 1) = 1
    k(ipu5, 2) = kfres
    k(ipu5, 3) = i
    p(ipu5, 5) = shr
    p(ipu5, 1) = -p(ipu3, 1) - p(ipu4, 1)
    p(ipu5, 2) = -p(ipu3, 2) - p(ipu4, 2)
    p(ipu5, 3) = -p(ipu3, 3) - p(ipu4, 3)
    p(ipu5, 4) = shpr - p(ipu3, 4) - p(ipu4, 4)
    k(i, 1) = 21
    k(i, 2) = kfres
    Do j = 1, 5
      p(i, j) = p(ipu5, j)
    End Do
    n = ipu5
    mint(23) = kfres
  Else If (idoc==12) Then
    phi(1) = paru(2)*rlu(0)
    phi(2) = phi(1) - phir
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      p(i, 5) = ulmass(k(i,2))
      If (0.5*shpr*z(jt)<=p(i,5)) p(i, 5) = 0.
      pabs = sqrt(max(0.,(0.5*shpr*z(jt))**2-p(i,5)**2))
      ptabs = pabs*sqrt(max(0.,1.-cthe(jt)**2))
      p(i, 1) = ptabs*cos(phi(jt))
      p(i, 2) = ptabs*sin(phi(jt))
      p(i, 3) = pabs*cthe(jt)*(-1)**(jt+1)
      p(i, 4) = 0.5*shpr*z(jt)
      izw = mint(83) + 6 + jt
      k(izw, 1) = 21
      If (mint(14+jt)==mint(20+jt)) Then
        k(izw, 2) = 23
      Else
        k(izw, 2) = isign(24, luchge(mint(14+jt))-luchge(mint(20+jt)))
      End If
      k(izw, 3) = izw - 2
      p(izw, 1) = -p(i, 1)
      p(izw, 2) = -p(i, 2)
      p(izw, 3) = (0.5*shpr-pabs*cthe(jt))*(-1)**(jt+1)
      p(izw, 4) = 0.5*shpr*(1.-z(jt))
      p(izw, 5) = -sqrt(max(0.,p(izw,3)**2+ptabs**2-p(izw,4)**2))
      ipu = mint(84) + 4 + jt
      k(ipu, 1) = 3
      k(ipu, 2) = kfpr(isub, jt)
      k(ipu, 3) = mint(83) + 8 + jt
      If (iabs(k(ipu,2))<=10 .Or. k(ipu,2)==21) Then
        p(ipu, 5) = ulmass(k(ipu,2))
      Else
        p(ipu, 5) = sqrt(vint(63+mod(js+jt,2)))
      End If
      mint(22+jt) = k(izw, 2)
    End Do
    If (isub==72) k(mint(84)+4+int(1.5+rlu(0)), 2) = -24
    i1 = mint(83) + 7
    i2 = mint(83) + 8
    bexcm = (p(i1,1)+p(i2,1))/(p(i1,4)+p(i2,4))
    beycm = (p(i1,2)+p(i2,2))/(p(i1,4)+p(i2,4))
    bezcm = (p(i1,3)+p(i2,3))/(p(i1,4)+p(i2,4))
    gamcm = (p(i1,4)+p(i2,4))/shr
    bepcm = bexcm*p(i1, 1) + beycm*p(i1, 2) + bezcm*p(i1, 3)
    px = p(i1, 1) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*bexcm
    py = p(i1, 2) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*beycm
    pz = p(i1, 3) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*bezcm
    thecm = ulangl(pz, sqrt(px**2+py**2))
    phicm = ulangl(px, py)
    sqlam = (sh-p(ipu5,5)**2-p(ipu6,5)**2)**2 - 4.*p(ipu5, 5)**2*p(ipu6, 5)**2
    pabs = sqrt(max(0.,sqlam/(4.*sh)))
    cthwz = vint(23)
    sthwz = sqrt(max(0.,1.-cthwz**2))
    phiwz = vint(24) - phicm
    p(ipu5, 1) = pabs*sthwz*cos(phiwz)
    p(ipu5, 2) = pabs*sthwz*sin(phiwz)
    p(ipu5, 3) = pabs*cthwz
    p(ipu5, 4) = sqrt(pabs**2+p(ipu5,5)**2)
    p(ipu6, 1) = -p(ipu5, 1)
    p(ipu6, 2) = -p(ipu5, 2)
    p(ipu6, 3) = -p(ipu5, 3)
    p(ipu6, 4) = sqrt(pabs**2+p(ipu6,5)**2)
    Call ludbrb(ipu5, ipu6, thecm, phicm, dble(bexcm), dble(beycm), dble(bezcm))
    Do jt = 1, 2
      i1 = mint(83) + 8 + jt
      i2 = mint(84) + 4 + jt
      k(i1, 1) = 21
      k(i1, 2) = k(i2, 2)
      Do j = 1, 5
        p(i1, j) = p(i2, j)
      End Do
    End Do
    n = ipu6
    mint(7) = mint(83) + 9
    mint(8) = mint(83) + 10
  End If
  If (idoc>=8) Then
    Do j = 1, 2
      jc = j
      If (kcs==-1) jc = 3 - j
      If (icol(kcc,1,jc)/=0 .And. k(ipu1,1)==14) k(ipu1, j+3) = k(ipu1, j+3) + mint(84) + icol(kcc, 1, jc)
      If (icol(kcc,2,jc)/=0 .And. k(ipu2,1)==14) k(ipu2, j+3) = k(ipu2, j+3) + mint(84) + icol(kcc, 2, jc)
      If (icol(kcc,3,jc)/=0 .And. k(ipu3,1)==3) k(ipu3, j+3) = mstu(5)*(mint(84)+icol(kcc,3,jc))
      If (icol(kcc,4,jc)/=0 .And. k(ipu4,1)==3) k(ipu4, j+3) = mstu(5)*(mint(84)+icol(kcc,4,jc))
    End Do
    Do i = 1, 2
      i1 = mint(83) + idoc - 2 + i
      i2 = mint(84) + 2 + i
      k(i1, 1) = 21
      k(i1, 2) = k(i2, 2)
      If (idoc<=9) k(i1, 3) = 0
      If (idoc>=11) k(i1, 3) = mint(83) + 2 + i
      Do j = 1, 5
        p(i1, j) = p(i2, j)
      End Do
    End Do
  End If
  mint(52) = n
  If (isub==95) Then
    k(ipu3, 1) = k(ipu3, 1) + 10
    k(ipu4, 1) = k(ipu4, 1) + 10
    Do j = 41, 66
      vint(j) = 0.
    End Do
    Do i = mint(83) + 5, mint(83) + 8
      Do j = 1, 5
        p(i, j) = 0.
      End Do
    End Do
  End If
  Return
End Subroutine pyscat
