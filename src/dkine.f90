Subroutine dkine(irun, i, nnn, nlab, iseed, wid, nt)
  Parameter (hbarc=0.19733)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, etam=0.5475, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  External iarflv, invflv
  Save
  px = p(1, i)
  py = p(2, i)
  pz = p(3, i)
  rx = r(1, i)
  ry = r(2, i)
  rz = r(3, i)
  dm = e(i)
  edelta = sqrt(dm**2+px**2+py**2+pz**2)
  pm = epion(nnn, irun)
  am = amp
  If (nlab==2) am = amn
  q2 = ((dm**2-am**2+pm**2)/(2.*dm))**2 - pm**2
  If (q2<=0.) q2 = 1.E-09
  q = sqrt(q2)
11 qx = 1. - 2.*ranart(nseed)
  qy = 1. - 2.*ranart(nseed)
  qz = 1. - 2.*ranart(nseed)
  qs = qx**2 + qy**2 + qz**2
  If (qs>1.) Goto 11
  pxp = q*qx/sqrt(qs)
  pyp = q*qy/sqrt(qs)
  pzp = q*qz/sqrt(qs)
  ep = sqrt(q**2+pm**2)
  pxn = -pxp
  pyn = -pyp
  pzn = -pzp
  en = sqrt(q**2+am**2)
  gd = edelta/dm
  fgd = gd/(1.+gd)
  bdx = px/edelta
  bdy = py/edelta
  bdz = pz/edelta
  bpp = bdx*pxp + bdy*pyp + bdz*pzp
  bpn = bdx*pxn + bdy*pyn + bdz*pzn
  p(1, i) = pxn + bdx*gd*(fgd*bpn+en)
  p(2, i) = pyn + bdy*gd*(fgd*bpn+en)
  p(3, i) = pzn + bdz*gd*(fgd*bpn+en)
  e(i) = am
  ppion(1, nnn, irun) = pxp + bdx*gd*(fgd*bpp+ep)
  ppion(2, nnn, irun) = pyp + bdy*gd*(fgd*bpp+ep)
  ppion(3, nnn, irun) = pzp + bdz*gd*(fgd*bpp+ep)
  dppion(nnn, irun) = dpertp(i)
  rpion(1, nnn, irun) = r(1, i)
  rpion(2, nnn, irun) = r(2, i)
  rpion(3, nnn, irun) = r(3, i)
  devio = sqrt(epion(nnn,irun)**2+ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2) - e1
  If (nt==ntmax) Then
     tau0 = hbarc/wid
     taudcy = tau0*(-1.)*alog(1.-ranart(nseed))
     taudcy = taudcy*e1/em1
     tfnl = tfnl + taudcy
     xfnl = xfnl + px1/e1*taudcy
     yfnl = yfnl + py1/e1*taudcy
     zfnl = zfnl + pz1/e1*taudcy
     r(1, i) = xfnl
     r(2, i) = yfnl
     r(3, i) = zfnl
     tfdcy(i) = tfnl
     rpion(1, nnn, irun) = xfnl
     rpion(2, nnn, irun) = yfnl
     rpion(3, nnn, irun) = zfnl
     tfdpi(nnn, irun) = tfnl
  End If
  Return
200 Format (A30, 2(1X,E10.4))
210 Format (I6, 5(1X,F8.3))
220 Format (A2, I5, 5(1X,F8.3))
End Subroutine dkine
