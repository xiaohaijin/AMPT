Subroutine hjana1
  Parameter (ymax=1.0, ymin=-1.0)
  Parameter (dmt=0.05, dy=0.2)
  Parameter (dr=0.2)
  Parameter (maxptn=400001, maxstr=150001)
  Dimension dyp1(50), dmyp1(200), deyp1(50)
  Dimension dyg1(50), dmyg1(200), deyg1(50)
  Dimension snyp1(50), smyp1(200), seyp1(50)
  Dimension snyg1(50), smyg1(200), seyg1(50)
  Dimension dnrpj1(50), dnrtg1(50), dnrin1(50), dnrtt1(50)
  Dimension dyg1c(50), dmyg1c(50), deyg1c(50)
  Dimension snrpj1(50), snrtg1(50), snrin1(50), snrtt1(50)
  Dimension snyg1c(50), smyg1c(50), seyg1c(50)
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0
  Common /para1/mul
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 200
      dmyp1(i) = smyp1(i)
      dmyg1(i) = smyg1(i)
    End Do
    Do i = 1, 50
      dyp1(i) = snyp1(i)
      deyp1(i) = seyp1(i)
      dyg1(i) = snyg1(i)
      deyg1(i) = seyg1(i)
      dnrpj1(i) = snrpj1(i)
      dnrtg1(i) = snrtg1(i)
      dnrin1(i) = snrin1(i)
      dnrtt1(i) = snrtt1(i)
      dyg1c(i) = snyg1c(i)
      dmyg1c(i) = smyg1c(i)
      deyg1c(i) = seyg1c(i)
    End Do
    nsubp = nsubps
    nsubg = nsubgs
    nisg = nisgs
  Else
    Do i = 1, 200
      smyp1(i) = dmyp1(i)
      smyg1(i) = dmyg1(i)
    End Do
    Do i = 1, 50
      snyp1(i) = dyp1(i)
      seyp1(i) = deyp1(i)
      snyg1(i) = dyg1(i)
      seyg1(i) = deyg1(i)
      snrpj1(i) = dnrpj1(i)
      snrtg1(i) = dnrtg1(i)
      snrin1(i) = dnrin1(i)
      snrtt1(i) = dnrtt1(i)
      snyg1c(i) = dyg1c(i)
      smyg1c(i) = dmyg1c(i)
      seyg1c(i) = deyg1c(i)
    End Do
    nsubps = nsubp
    nsubgs = nsubg
    nisgs = nisg
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
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
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 100
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      100 Continue
      imt = 1 + int(dxmt/dmt)
      If (rap>ymax .Or. rap<=ymin) Goto 200
      If (imt>200) Goto 200
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
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
        Print *, ' IN HJANA1 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 300
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      300 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 400
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 400
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
      End If
      400 Continue
    End Do
  End Do
  Do i = 1, nsg
    Do j = 1, njsg(i)
      ityp = k2sg(i, j)
      px = pxsg(i, j)
      py = pysg(i, j)
      pz = pzsg(i, j)
      pe = pesg(i, j)
      pm = pmsg(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA1 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(abs(rap)/dy)
      If (iy<1 .Or. iy>50) Goto 500
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      500 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 600
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 600
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
      End If
      600 Continue
    End Do
  End Do
  Do i = 1, ihnt2(1)
    yr = sqrt(yp(1,i)**2+yp(2,i)**2)
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 601
    dnrpj1(ir) = dnrpj1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    601 Continue
  End Do
  Do i = 1, ihnt2(3)
    yr = sqrt(yt(1,i)**2+yt(2,i)**2)
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 602
    dnrtg1(ir) = dnrtg1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    602 Continue
  End Do
  Do i = 1, nsg
    y1 = 0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
    y2 = 0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
    yr = sqrt(y1**2+y2**2)
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 603
    dnrin1(ir) = dnrin1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    603 Continue
  End Do
  Do i = 1, mul
    ityp = ityp0(i)
    px = sngl(px0(i))
    py = sngl(py0(i))
    pz = sngl(pz0(i))
    pe = sngl(e0(i))
    pm = sngl(xmass0(i))
    xmt = sqrt(px**2+py**2+pm**2)
    dxmt = xmt - pm
    If (xmt>0.) Then
      rap = asinh(pz/xmt)
    Else
      Print *, ' IN HJANA1 mt=0'
      rap = 1000000.0*sign(1., pz)
    End If
    iy = 1 + int(abs(rap)/dy)
    If (iy<1 .Or. iy>50) Goto 700
    dyg1c(iy) = dyg1c(iy) + 1.0
    deyg1c(iy) = deyg1c(iy) + xmt
    700 Continue
    If (rap>ymax .Or. rap<=ymin) Goto 800
    imt = 1 + int(dxmt/dmt)
    If (imt>50) Goto 800
    dmyg1c(imt) = dmyg1c(imt) + 1.0/xmt
    800 Continue
  End Do
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
  Do i = 1, nsg
    Do j = 1, njsg(i)
      nsubp = nsubp + 1
      If (k2sg(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do
  nisg = nisg + nsg
  If (iout==1) Then
    Print *, ' in HJANA1 '
    Print *, ' total number of partons = ', nsubp/iw
    Print *, ' total number of gluons = ', nsubg/iw
    Print *, ' number of independent strings = ', nisg/iw
  End If
  Return
End Subroutine hjana1
