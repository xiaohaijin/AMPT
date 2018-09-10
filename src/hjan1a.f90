Subroutine hjan1a
  Parameter (maxptn=400001)
  Parameter (dgx=0.2, dgy=0.2, dt=0.2)
  Dimension dgxg1a(50), dgyg1a(50), dtg1a(50)
  Dimension sgxg1a(50), sgyg1a(50), stg1a(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dgxg1a(i) = sgxg1a(i)
      dgyg1a(i) = sgyg1a(i)
      dtg1a(i) = stg1a(i)
    End Do
  Else
    Do i = 1, 50
      sgxg1a(i) = dgxg1a(i)
      sgyg1a(i) = dgyg1a(i)
      stg1a(i) = dtg1a(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
  Do i = 1, mul
    igx = 1 + int(sngl(abs(gx5(i)))/dgx)
    If (igx>50 .Or. igx<1) Goto 100
    dgxg1a(igx) = dgxg1a(igx) + 1.0
    100 Continue
    igy = 1 + int(sngl(abs(gy5(i)))/dgy)
    If (igy>50 .Or. igy<1) Goto 200
    dgyg1a(igy) = dgyg1a(igy) + 1.0
    200 Continue
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '1:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      it = 1
    Else
      it = 1 + int(sqrt(diff2)/dt)
    End If
    If (it>50 .Or. it<1) Goto 300
    dtg1a(it) = dtg1a(it) + 1.0
    300 Continue
  End Do
  Call hjan1b
  Return
End Subroutine hjan1a
