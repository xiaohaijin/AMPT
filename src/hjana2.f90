Subroutine hjana2
  Parameter (ymax=1.0, ymin=-1.0)
  Parameter (dmt=0.05, dy=0.2)
  Parameter (dr=0.2, dt=0.2)
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Dimension dyp2(50), dmyp2(200), deyp2(50)
  Dimension dyg2(50), dmyg2(200), deyg2(50)
  Dimension snyp2(50), smyp2(200), seyp2(50)
  Dimension snyg2(50), smyg2(200), seyg2(50)
  Dimension dnrpj2(50), dnrtg2(50), dnrin2(50), dnrtt2(50)
  Dimension dtpj2(50), dttg2(50), dtin2(50), dttot2(50)
  Dimension dyg2c(50), dmyg2c(50), deyg2c(50)
  Dimension snrpj2(50), snrtg2(50), snrin2(50), snrtt2(50)
  Dimension stpj2(50), sttg2(50), stin2(50), sttot2(50)
  Dimension snyg2c(50), smyg2c(50), seyg2c(50)
  Double Precision ataui, zt1, zt2, zt3
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Common /anim/nevent, isoft, isflag, izpc
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 200
      dmyp2(i) = smyp2(i)
      dmyg2(i) = smyg2(i)
    End Do
    Do i = 1, 50
      dyp2(i) = snyp2(i)
      deyp2(i) = seyp2(i)
      dyg2(i) = snyg2(i)
      deyg2(i) = seyg2(i)
      dnrpj2(i) = snrpj2(i)
      dnrtg2(i) = snrtg2(i)
      dnrin2(i) = snrin2(i)
      dnrtt2(i) = snrtt2(i)
      dtpj2(i) = stpj2(i)
      dttg2(i) = sttg2(i)
      dtin2(i) = stin2(i)
      dttot2(i) = sttot2(i)
      dyg2c(i) = snyg2c(i)
      dmyg2c(i) = smyg2c(i)
      deyg2c(i) = seyg2c(i)
    End Do
    nsubp = nsubps
    nsubg = nsubgs
    nisg = nisgs
  Else
    Do i = 1, 200
      smyp2(i) = dmyp2(i)
      smyg2(i) = dmyg2(i)
    End Do
    Do i = 1, 50
      snyp2(i) = dyp2(i)
      seyp2(i) = deyp2(i)
      snyg2(i) = dyg2(i)
      seyg2(i) = deyg2(i)
      snrpj2(i) = dnrpj2(i)
      snrtg2(i) = dnrtg2(i)
      snrin2(i) = dnrin2(i)
      snrtt2(i) = dnrtt2(i)
      stpj2(i) = dtpj2(i)
      sttg2(i) = dttg2(i)
      stin2(i) = dtin2(i)
      sttot2(i) = dttot2(i)
      snyg2c(i) = dyg2c(i)
      smyg2c(i) = dmyg2c(i)
      seyg2c(i) = deyg2c(i)
    End Do
    nsubps = nsubp
    nsubgs = nsubg
    nisgs = nisg
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 510
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      ityp = kfpj(i, j)
      px = pjpx(i, j)
      py = pjpy(i, j)
      pz = pjpz(i, j)
      pe = pjpe(i, j)
      pm = pjpm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 100
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      100 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 200
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 200
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      200 Continue
    End Do
  End Do
  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      ityp = kftj(i, j)
      px = pjtx(i, j)
      py = pjty(i, j)
      pz = pjtz(i, j)
      pe = pjte(i, j)
      pm = pjtm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 300
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      300 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 400
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 400
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      400 Continue
    End Do
  End Do
  510 Continue
  Do i = 1, nsg
    nj = njsg(i)
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) nj = njsgs(i)
    Do j = 1, nj
      ityp = k2sg(i, j)
      px = pxsg(i, j)
      py = pysg(i, j)
      pz = pzsg(i, j)
      pe = pesg(i, j)
      pm = pmsg(i, j)
      If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        ityp = k2sgs(i, j)
        px = sngl(pxsgs(i,j))
        py = sngl(pysgs(i,j))
        pz = sngl(pzsgs(i,j))
        pe = sngl(pesgs(i,j))
        pm = sngl(pmsgs(i,j))
      End If
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 500
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      500 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 600
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 600
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      600 Continue
    End Do
  End Do
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 520
  Do i = 1, ihnt2(1)
    j = i
    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 601
    dnrpj2(ir) = dnrpj2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    601 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 602
    dtpj2(it) = dtpj2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    602 Continue
  End Do
  Do i = 1, ihnt2(3)
    j = i + ihnt2(1)
    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 603
    dnrtg2(ir) = dnrtg2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    603 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 604
    dttg2(it) = dttg2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    604 Continue
  End Do
  520 Continue
  Do i = 1, nsg
    j = i + ihnt2(1) + ihnt2(3)
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) j = i
    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 605
    dnrin2(ir) = dnrin2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    605 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 606
    dtin2(it) = dtin2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    606 Continue
  End Do
  Do i = 1, mul
    ityp = ityp5(i)
    px = sngl(px5(i))
    py = sngl(py5(i))
    pz = sngl(pz5(i))
    pe = sngl(e5(i))
    pm = sngl(xmass5(i))
    xmt = sqrt(px**2+py**2+pm**2)
    dxmt = xmt - pm
    If (xmt>0.) Then
      rap = asinh(pz/xmt)
    Else
      Print *, ' IN HJANA2 mt=0'
      rap = 1000000.0*sign(1., pz)
    End If
    iy = 1 + int(abs(rap)/dy)
    If (iy<1 .Or. iy>50) Goto 700
    dyg2c(iy) = dyg2c(iy) + 1.0
    deyg2c(iy) = deyg2c(iy) + xmt
    700 Continue
    If (rap>ymax .Or. rap<=ymin) Goto 800
    imt = 1 + int(dxmt/dmt)
    If (imt>50) Goto 800
    dmyg2c(imt) = dmyg2c(imt) + 1.0/xmt
    800 Continue
  End Do
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 530
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      nsubp = nsubp + 1
      If (kfpj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do
  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      nsubp = nsubp + 1
      If (kftj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do
  530 Continue
  Do i = 1, nsg
    nj = njsg(i)
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) nj = njsgs(i)
    Do j = 1, nj
      nsubp = nsubp + 1
      If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        If (k2sgs(i,j)==21) nsubg = nsubg + 1
      Else
        If (k2sg(i,j)==21) nsubg = nsubg + 1
      End If
    End Do
  End Do
  nisg = nisg + nsg
  If (iout==1) Then
    Print *, ' in HJANA2 '
    Print *, ' total number of partons = ', nsubp/iw
    Print *, ' total number of gluons = ', nsubg/iw
    Print *, ' number of independent strings = ', nisg/iw
  End If
  Call hjan2a
  Call hjan2b
  Return
End Subroutine hjana2
