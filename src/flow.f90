Subroutine flow(nt)
  Parameter (pi=3.1415926, apion=0.13957, aka=0.498)
  Parameter (maxstr=150001, maxr=1, amu=0.9383, etam=0.5475)
  Dimension ypion(-80:80), ypr(-80:80), ykaon(-80:80)
  Dimension pxpion(-80:80), pxpro(-80:80), pxkaon(-80:80)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /rr/massr(0:maxr)
  Common /run/num
  Common /input1/masspr, massta, iseed, iavoid, dt
  Save
  ycut1 = -2.6
  ycut2 = 2.6
  dy = 0.2
  ly = nint((ycut2-ycut1)/dy)
  Do kk = -80, 80
    pxpion(kk) = 0
    pxpro(kk) = 0
    pxkaon(kk) = 0
  End Do
  Do j = -ly, ly
    ypion(j) = 0
    ykaon(j) = 0
    ypr(j) = 0
  End Do
  nkaon = 0
  npr = 0
  npion = 0
  is = 0
  Do nrun = 1, num
    is = is + massr(nrun-1)
    Do j = 1, massr(nrun)
      i = j + is
      e00 = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
      y00 = 0.5*alog((e00+p(3,i))/(e00-p(3,i)))
      If (abs(y00)>=ycut2) Goto 20
      iy = nint(y00/dy)
      If (abs(iy)>=80) Goto 20
      If (e(i)==0) Goto 20
      If (lb(i)>=25) Goto 20
      If ((lb(i)<=5) .And. (lb(i)>=3)) Goto 50
      If (lb(i)==1 .Or. lb(i)==2) Goto 200
      If (lb(i)>=6 .And. lb(i)<=17) Goto 200
      If (lb(i)==23) Goto 400
      Goto 20
      50 npion = npion + 1
      ypion(iy) = ypion(iy) + 1
      pxpion(iy) = pxpion(iy) + p(1, i)/e(i)
      Goto 20
      200 npr = npr + 1
      pxpro(iy) = pxpro(iy) + p(1, i)/e(i)
      ypr(iy) = ypr(iy) + 1.
      Goto 20
      400 nkaon = nkaon + 1
      ykaon(iy) = ykaon(iy) + 1.
      pxkaon(iy) = pxkaon(iy) + p(1, i)/e(i)
    20 End Do
  End Do
  Do npt = -10, 10
    If (ypr(npt)==0) Goto 101
    pxpro(npt) = -pxpro(npt)/ypr(npt)
    dnuc = pxpro(npt)/sqrt(ypr(npt))
    101 If (ypion(npt)==0) Goto 102
    pxpion(npt) = -pxpion(npt)/ypion(npt)
    dnucp = pxpion(npt)/sqrt(ypion(npt))
    102 If (ykaon(npt)==0) Goto 3
    pxkaon(npt) = -pxkaon(npt)/ykaon(npt)
    dnuck = pxkaon(npt)/sqrt(ykaon(npt))
  3 End Do
  Do m = -ly, ly
    dypr = 0
    If (ypr(m)/=0) dypr = sqrt(ypr(m))/float(nrun)/dy
    ypr(m) = ypr(m)/float(nrun)/dy
    dypion = 0
    If (ypion(m)/=0) dypion = sqrt(ypion(m))/float(nrun)/dy
    ypion(m) = ypion(m)/float(nrun)/dy
    dykaon = 0
    If (ykaon(m)/=0) dykaon = sqrt(ykaon(m))/float(nrun)/dy
    ykaon(m) = ykaon(m)/float(nrun)/dy
  End Do
  Return
End Subroutine flow
