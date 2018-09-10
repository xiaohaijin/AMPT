Subroutine nnkaon(irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
  Parameter (maxstr=150001, maxr=1)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /nn/nnn
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Dimension px(4), py(4), pz(4)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  dm3 = 0.938
  dm4 = 0.938
  n = 0
  If (nchrg<=-1 .Or. nchrg>=3) dm3 = 1.232
  If (nchrg==-2 .Or. nchrg==4) dm4 = 1.232
  iblock = 0
  Call fstate(iseed, srt, dm3, dm4, px, py, pz, iflag)
  If (iflag<0) Then
    ictrl = -1
    n = n + 1
    Return
  End If
  iblock = 12
  pxrota = px(1)
  pyrota = py(1)
  pzrota = pz(1)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(1) = pxrota
  py(1) = pyrota
  pz(1) = pzrota
  pxrota = px(2)
  pyrota = py(2)
  pzrota = pz(2)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(2) = pxrota
  py(2) = pyrota
  pz(2) = pzrota
  pxrota = px(3)
  pyrota = py(3)
  pzrota = pz(3)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(3) = pxrota
  py(3) = pyrota
  pz(3) = pzrota
  pxrota = px(4)
  pyrota = py(4)
  pzrota = pz(4)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(4) = pxrota
  py(4) = pyrota
  pz(4) = pzrota
  nnn = nnn + 2
  lpion(nnn, irun) = 23
  If (nchrg==-1 .Or. nchrg==-2) Then
  End If
  epion(nnn, irun) = aka
  lpion(nnn-1, irun) = 21
  epion(nnn-1, irun) = aka
  e1cm = sqrt(dm3**2+px(1)**2+py(1)**2+pz(1)**2)
  p1beta = px(1)*betax + py(1)*betay + pz(1)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px(1)
  pt2i1 = betay*transf + py(1)
  pt3i1 = betaz*transf + pz(1)
  eti1 = dm3
  lb1 = 2
  If (nchrg>=-2 .And. nchrg<=1) lb1 = 2
  If (nchrg==-2 .Or. nchrg==-1) Then
    lb1 = 6
  End If
  If (nchrg==1 .Or. nchrg==2) lb1 = 1
  If (nchrg==3 .Or. nchrg==4) lb1 = 9
  e2cm = sqrt(dm4**2+px(2)**2+py(2)**2+pz(2)**2)
  p2beta = px(2)*betax + py(2)*betay + pz(2)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + px(2)
  pt2i2 = betay*transf + py(2)
  pt3i2 = betaz*transf + pz(2)
  eti2 = dm4
  lb2 = 2
  If (nchrg>=-1 .Or. nchrg<=1) lb2 = 2
  If (nchrg==2 .Or. nchrg==3) lb2 = 1
  If (nchrg==4) lb2 = 9
  If (nchrg==-2) lb2 = 6
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2
  epcmk = sqrt(epion(nnn-1,irun)**2+px(3)**2+py(3)**2+pz(3)**2)
  betak = px(3)*betax + py(3)*betay + pz(3)*betaz
  transf = gamma*(gamma*betak/(gamma+1.)+epcmk)
  ppion(1, nnn-1, irun) = betax*transf + px(3)
  ppion(2, nnn-1, irun) = betay*transf + py(3)
  ppion(3, nnn-1, irun) = betaz*transf + pz(3)
  rpion(1, nnn-1, irun) = r(1, i1)
  rpion(2, nnn-1, irun) = r(2, i1)
  rpion(3, nnn-1, irun) = r(3, i1)
  dppion(nnn-1, irun) = dpertp(i1)*dpertp(i2)
  epcmak = sqrt(epion(nnn,irun)**2+px(4)**2+py(4)**2+pz(4)**2)
  betaak = px(4)*betax + py(4)*betay + pz(4)*betaz
  transf = gamma*(gamma*betaak/(gamma+1.)+epcmak)
  ppion(1, nnn, irun) = betax*transf + px(4)
  ppion(2, nnn, irun) = betay*transf + py(4)
  ppion(3, nnn, irun) = betaz*transf + pz(4)
  rpion(1, nnn, irun) = r(1, i2)
  rpion(2, nnn, irun) = r(2, i2)
  rpion(3, nnn, irun) = r(3, i2)
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  Return
End Subroutine nnkaon
