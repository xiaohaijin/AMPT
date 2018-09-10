Subroutine hjan2b
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Parameter (dr=0.2, dt=0.2)
  Dimension dnrg2b(50), dtg2b(-24:25)
  Dimension snrg2b(50), stg2b(-24:25)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Double Precision ataui, zt1, zt2, zt3
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  Common /srec1/nsp, nst, nsi
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dnrg2b(i) = snrg2b(i)
      dtg2b(i-25) = stg2b(i-25)
    End Do
  Else
    Do i = 1, 50
      snrg2b(i) = dnrg2b(i)
      stg2b(i-25) = dtg2b(i-25)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
  Do i = 1, mul
    j = lstrg1(i)
    gx0 = sngl(zt1(j))
    gy0 = sngl(zt2(j))
    r0 = sqrt((sngl(gx5(i))-gx0)**2+(sngl(gy5(i))-gy0)**2)
    ir = 1 + int(r0/dr)
    If (ir>50 .Or. ir<1) Goto 100
    dnrg2b(ir) = dnrg2b(ir) + 1.0
    100 Continue
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '4:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1E-6
    Else
      tau7 = sqrt(diff2)
    End If
    dtau = tau7 - sngl(ataui(j))
    it = 1 + int(dtau/dt)
    If (it>25 .Or. it<-24) Goto 200
    dtg2b(it) = dtg2b(it) + 1.0
    200 Continue
  End Do
  Return
End Subroutine hjan2b
