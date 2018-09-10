Subroutine zpstrg
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Real yp, yt, pxsg, pysg, pzsg, pesg, pmsg, hipr1, hint1, bb
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  Common /srec1/nsp, nst, nsi
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /anim/nevent, isoft, isflag, izpc
  Common /strg/np(maxstr)
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  Save
  If (isoft==5) Then
    Do i = 1, mul
      ityp5(i) = idfrz(i)
      gx5(i) = gxfrz(i)
      gy5(i) = gyfrz(i)
      gz5(i) = gzfrz(i)
      ft5(i) = ftfrz(i)
      px5(i) = pxfrz(i)
      py5(i) = pyfrz(i)
      pz5(i) = pzfrz(i)
      e5(i) = efrz(i)
      xmass5(i) = xmfrz(i)
    End Do
  End If
  Do i = 1, maxstr
    ataui(i) = 0D0
    zt1(i) = 0D0
    zt2(i) = 0D0
    zt3(i) = 0D0
    np(i) = 0
  End Do
  Do i = 1, mul
    istrg = lstrg1(i)
    diff2 = ft5(i)**2 - gz5(i)**2
    If (diff2<0D0) Then
      Write (6, *) '2:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1D-6
    Else
      tau7 = dsqrt(diff2)
    End If
    ataui(istrg) = ataui(istrg) + tau7
    zt1(istrg) = zt1(istrg) + gx5(i)
    zt2(istrg) = zt2(istrg) + gy5(i)
    zt3(istrg) = zt3(istrg) + gz5(i)
    np(istrg) = np(istrg) + 1
  End Do
  nstr = nsp + nst + nsi
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    Do i = 1, nstr
      If (np(i)/=0) Then
        ataui(i) = ataui(i)/np(i)
        zt1(i) = zt1(i)/np(i)
        zt2(i) = zt2(i)/np(i)
        zt3(i) = zt3(i)/np(i)
      End If
    End Do
    Return
  End If
  Do i = 1, nstr
    If (np(i)/=0) Then
      ataui(i) = ataui(i)/np(i)
      zt1(i) = zt1(i)/np(i)
      zt2(i) = zt2(i)/np(i)
      zt3(i) = zt3(i)/np(i)
    Else
      If (i<=nsp) Then
        j = i
        zt1(i) = dble(yp(1,j))
        zt2(i) = dble(yp(2,j))
        zt3(i) = 0D0
      Else If (i>nsp .And. i<=nsp+nst) Then
        j = i - nsp
        zt1(i) = dble(yt(1,j))
        zt2(i) = dble(yt(2,j))
        zt3(i) = 0D0
      Else
        j = i - nsp - nst
        zt1(i) = 0.5D0*dble((yp(1,iasg(j,1))+yt(1,iasg(j,2))))
        zt2(i) = 0.5D0*dble((yp(2,iasg(j,1))+yt(2,iasg(j,2))))
        zt3(i) = 0D0
      End If
    End If
  End Do
  bb = hint1(19)
  Do i = 1, nstr
    If (np(i)/=0) Then
      shift = 0D0
    Else
      shift = 0.5D0*dble(bb)
    End If
    If (i<=nsp) Then
      zt1(i) = zt1(i) + shift
    Else If (i>nsp .And. i<=nsp+nst) Then
      zt1(i) = zt1(i) - shift
    End If
  End Do
  Return
End Subroutine zpstrg
