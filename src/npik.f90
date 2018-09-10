Subroutine npik(irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, ratiok, iblock)
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
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
  Dimension bb(3), p1(4), p2(4), p3(4), px(4), py(4), pz(4)
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ictrl = 1
  lb1 = lb(i1)
  lb2 = lb(i2)
  k1 = i1
  k2 = i2
  If (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13)) Then
    k1 = i2
    k2 = i1
  End If
  lb(k1) = 1 + int(2*ranart(nseed))
  lb(k2) = 23
  pkmax = sqrt((srt**2-(aka+0.938+aka)**2)*(srt**2-(aka+0.938-aka)**2))/2./srt
  pk = ranart(nseed)*pkmax
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p3(1) = pk*sss*cos(fai)
  p3(2) = pk*sss*sin(fai)
  p3(3) = pk*css
  eip = srt - sqrt(aka**2+pk**2)
  rmnp = sqrt(eip**2-pk**2)
  Do i = 1, 3
    bb(i) = -1.*p3(i)/eip
  End Do
  pznp = sqrt((rmnp**2-(aka+0.938)**2)*(rmnp**2-(0.938-aka)**2))/2./rmnp
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p1(1) = pznp*sss*cos(fai)
  p1(2) = pznp*sss*sin(fai)
  p1(3) = pznp*css
  p1(4) = sqrt(0.938**2+pznp**2)
  p2(4) = sqrt(aka**2+pznp**2)
  Do i = 1, 3
    p2(i) = -1.*p1(i)
  End Do
  ilo = 1
  Call lorntz(ilo, bb, p1, p2)
  pxrota = p1(1)
  pyrota = p1(2)
  pzrota = p1(3)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p1(1) = pxrota
  p1(2) = pyrota
  p1(3) = pzrota
  pxrota = p2(1)
  pyrota = p2(2)
  pzrota = p2(3)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p2(1) = pxrota
  p2(2) = pyrota
  p2(3) = pzrota
  pxrota = p3(1)
  pyrota = p3(2)
  pzrota = p3(3)
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p3(1) = pxrota
  p3(2) = pyrota
  p3(3) = pzrota
  nnn = nnn + 1
  lpion(nnn, irun) = 21
  epion(nnn, irun) = aka
  e1cm = sqrt(0.938**2+p1(1)**2+p1(2)**2+p1(3)**2)
  p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + p1(1)
  pt2i1 = betay*transf + p1(2)
  pt3i1 = betaz*transf + p1(3)
  eti1 = 0.938
  lb1 = lb(k1)
  e2cm = sqrt(aka**2+p3(1)**2+p3(2)**2+p3(3)**2)
  p2beta = p3(1)*betax + p3(2)*betay + p3(3)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + p3(1)
  pt2i2 = betay*transf + p3(2)
  pt3i2 = betaz*transf + p3(3)
  eti2 = aka
  lb2 = lb(k2)
  p(1, k1) = pt1i1
  p(2, k1) = pt2i1
  p(3, k1) = pt3i1
  e(k1) = eti1
  lb(k1) = lb1
  p(1, k2) = pt1i2
  p(2, k2) = pt2i2
  p(3, k2) = pt3i2
  e(k2) = eti2
  lb(k2) = lb2
  iblock = 101
  epcmk = sqrt(epion(nnn,irun)**2+p2(1)**2+p2(2)**2+p2(3)**2)
  betak = p2(1)*betax + p2(2)*betay + p2(3)*betaz
  transf = gamma*(gamma*betak/(gamma+1.)+epcmk)
  ppion(1, nnn, irun) = betax*transf + p2(1)
  ppion(2, nnn, irun) = betay*transf + p2(2)
  ppion(3, nnn, irun) = betaz*transf + p2(3)
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  k = i2
  If (lb(i1)==1 .Or. lb(i1)==2) k = i1
  rpion(1, nnn, irun) = r(1, k)
  rpion(2, nnn, irun) = r(2, k)
  rpion(3, nnn, irun) = r(3, k)
  Return
End Subroutine npik
