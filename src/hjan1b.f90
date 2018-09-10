Subroutine hjan1b
  Parameter (maxptn=400001, maxstr=150001)
  Parameter (dr=0.2, dt=0.2)
  Dimension dnrg1b(50), dtg1b(50)
  Dimension snrg1b(50), stg1b(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  Common /srec1/nsp, nst, nsi
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arevt/iaevt, iarun, miss
  Common /arout/iout
  Save
  Data iw/0/
  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dnrg1b(i) = snrg1b(i)
      dtg1b(i) = stg1b(i)
    End Do
  Else
    Do i = 1, 50
      snrg1b(i) = dnrg1b(i)
      stg1b(i) = dtg1b(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
  Do i = 1, mul
    j = lstrg1(i)
    If (j<=nsp) Then
      k = j
      gx0 = yp(1, j)
      gy0 = yp(2, j)
    Else If (j<=nsp+nst) Then
      k = j - nsp
      gx0 = yt(1, k)
      gy0 = yt(2, k)
    Else
      k = j - nsp - nst
      gx0 = 0.5*(yp(1,iasg(k,1))+yt(1,iasg(k,2)))
      gy0 = 0.5*(yp(2,iasg(k,1))+yt(2,iasg(k,2)))
    End If
    r0 = sqrt((sngl(gx5(i))-gx0)**2+(sngl(gy5(i))-gy0)**2)
    ir = 1 + int(r0/dr)
    If (ir>50 .Or. ir<1) Goto 100
    dnrg1b(ir) = dnrg1b(ir) + 1.0
    100 Continue
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '5:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1E-6
    Else
      tau7 = sqrt(diff2)
    End If
    it = 1 + int(tau7/dt)
    If (it>50 .Or. it<1) Goto 200
    dtg1b(it) = dtg1b(it) + 1.0
    200 Continue
  End Do
  Return
End Subroutine hjan1b
