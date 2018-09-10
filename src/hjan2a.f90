Subroutine hjan2a
  Parameter (dgx=0.2, dgy=0.2, dt=0.2)
  Parameter (maxptn=400001, maxstr=150001)
  Dimension dgxp2a(50), dgyp2a(50), dtp2a(50)
  Dimension dgxg2a(50), dgyg2a(50), dtg2a(50)
  Dimension sgxp2a(50), sgyp2a(50), stp2a(50)
  Dimension sgxg2a(50), sgyg2a(50), stg2a(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dgxp2a(i) = sgxp2a(i)
      dgyp2a(i) = sgyp2a(i)
      dtp2a(i) = stp2a(i)
      dgxg2a(i) = sgxg2a(i)
      dgyg2a(i) = sgyg2a(i)
      dtg2a(i) = stg2a(i)
    End Do
  Else
    Do i = 1, 50
      sgxp2a(i) = dgxp2a(i)
      sgyp2a(i) = dgyp2a(i)
      stp2a(i) = dtp2a(i)
      sgxg2a(i) = dgxg2a(i)
      sgyg2a(i) = dgyg2a(i)
      stg2a(i) = dtg2a(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      If (kfpj(i,j)/=21) Then
        igx = 1 + int(abs(yp(1,i))/dgx)
        If (igx>50 .Or. igx<1) Goto 100
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        100 Continue
        igy = 1 + int(abs(yp(2,i))/dgy)
        If (igy>50 .Or. igy<1) Goto 200
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        200 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do
  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      If (kftj(i,j)/=21) Then
        igx = 1 + int(abs(yt(1,i))/dgx)
        If (igx>50 .Or. igx<1) Goto 300
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        300 Continue
        igy = 1 + int(abs(yt(2,i))/dgy)
        If (igy>50 .Or. igy<1) Goto 400
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        400 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do
  Do i = 1, nsg
    Do j = 1, njsg(i)
      If (k2sg(i,j)/=21) Then
        igx = 1 + int(abs(0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2))))/dgx)
        If (igx>50 .Or. igx<1) Goto 500
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        500 Continue
        igy = 1 + int(abs(0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2))))/dgy)
        If (igy>50 .Or. igy<1) Goto 600
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        600 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do
  Do i = 1, mul
    igx = 1 + int(abs(sngl(gx5(i)))/dgx)
    If (igx>50 .Or. igx<1) Goto 700
    dgxg2a(igx) = dgxg2a(igx) + 1.0
    dgxp2a(igx) = dgxp2a(igx) + 1.0
    700 Continue
    igy = 1 + int(abs(sngl(gy5(i)))/dgy)
    If (igy>50 .Or. igy<1) Goto 800
    dgyg2a(igy) = dgyg2a(igy) + 1.0
    dgyp2a(igy) = dgyp2a(igy) + 1.0
    800 Continue
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '3:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      it = 1
    Else
      it = 1 + int(sqrt(diff2)/dt)
    End If
    If (it>50 .Or. it<1) Goto 900
    dtg2a(it) = dtg2a(it) + 1.0
    dtp2a(it) = dtp2a(it) + 1.0
    900 Continue
  End Do
  Return
End Subroutine hjan2a
