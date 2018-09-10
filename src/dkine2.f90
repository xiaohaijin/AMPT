Subroutine dkine2(irun, i, nnn, nlab, iseed, wid, nt)
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
  External iarflv, invflv
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  px = p(1, i)
  py = p(2, i)
  pz = p(3, i)
  rx = r(1, i)
  ry = r(2, i)
  rz = r(3, i)
  dm = e(i)
  edelta = sqrt(dm**2+px**2+py**2+pz**2)
  pm1 = epion(nnn, irun)
  pm2 = epion(nnn+1, irun)
  am = amn
  If (nlab==1) am = amp
  pmax2 = (dm**2-(am+pm1+pm2)**2)*(dm**2-(am-pm1-pm2)**2)/4/dm**2
  scheck = pmax2
  If (scheck<0) Then
     Write (99, *) 'scheck15: ', scheck
     scheck = 0.
  End If
  pmax = sqrt(scheck)
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1-css**2)
  fai = 2*pi*ranart(nseed)
  px0 = pmax*sss*cos(fai)
  py0 = pmax*sss*sin(fai)
  pz0 = pmax*css
  ep0 = sqrt(px0**2+py0**2+pz0**2+am**2)
  betax = -px0/(dm-ep0)
  betay = -py0/(dm-ep0)
  betaz = -pz0/(dm-ep0)
  scheck = 1 - betax**2 - betay**2 - betaz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck16: ', scheck
     Stop
  End If
  gd1 = 1./sqrt(scheck)
  fgd1 = gd1/(1+gd1)
  q2 = ((dm-ep0)/(2.*gd1))**2 - pm1**2
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
  ep = sqrt(q**2+pm1**2)
  pxn = -pxp
  pyn = -pyp
  pzn = -pzp
  en = sqrt(q**2+pm2**2)
  bpp1 = betax*pxp + betay*pyp + betaz*pzp
  bpn1 = betax*pxn + betay*pyn + betaz*pzn
  p1m = pxn + betax*gd1*(fgd1*bpn1+en)
  p2m = pyn + betay*gd1*(fgd1*bpn1+en)
  p3m = pzn + betaz*gd1*(fgd1*bpn1+en)
  epn = sqrt(p1m**2+p2m**2+p3m**2+pm2**2)
  p1p = pxp + betax*gd1*(fgd1*bpp1+ep)
  p2p = pyp + betay*gd1*(fgd1*bpp1+ep)
  p3p = pzp + betaz*gd1*(fgd1*bpp1+ep)
  epp = sqrt(p1p**2+p2p**2+p3p**2+pm1**2)
  gd = edelta/dm
  fgd = gd/(1.+gd)
  bdx = px/edelta
  bdy = py/edelta
  bdz = pz/edelta
  bp0 = bdx*px0 + bdy*py0 + bdz*pz0
  bpp = bdx*p1p + bdy*p2p + bdz*p3p
  bpn = bdx*p1m + bdy*p2m + bdz*p3m
  p(1, i) = px0 + bdx*gd*(fgd*bp0+ep0)
  p(2, i) = py0 + bdy*gd*(fgd*bp0+ep0)
  p(3, i) = pz0 + bdz*gd*(fgd*bp0+ep0)
  e(i) = am
  id(i) = 0
  enucl = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
  ppion(1, nnn, irun) = p1p + bdx*gd*(fgd*bpp+epp)
  ppion(2, nnn, irun) = p2p + bdy*gd*(fgd*bpp+epp)
  ppion(3, nnn, irun) = p3p + bdz*gd*(fgd*bpp+epp)
  epion1 = sqrt(ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2+epion(nnn,irun)**2)
  rpion(1, nnn, irun) = r(1, i)
  rpion(2, nnn, irun) = r(2, i)
  rpion(3, nnn, irun) = r(3, i)
  ppion(1, nnn+1, irun) = p1m + bdx*gd*(fgd*bpn+epn)
  ppion(2, nnn+1, irun) = p2m + bdy*gd*(fgd*bpn+epn)
  ppion(3, nnn+1, irun) = p3m + bdz*gd*(fgd*bpn+epn)
  dppion(nnn, irun) = dpertp(i)
  dppion(nnn+1, irun) = dpertp(i)
  epion2 = sqrt(ppion(1,nnn+1,irun)**2+ppion(2,nnn+1,irun)**2+ppion(3,nnn+1,irun)**2+epion(nnn+1,irun)**2)
  rpion(1, nnn+1, irun) = r(1, i)
  rpion(2, nnn+1, irun) = r(2, i)
  rpion(3, nnn+1, irun) = r(3, i)
  devio = sqrt(epion(nnn,irun)**2+ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2) + sqrt(epion(nnn+1,irun)**2+ppion(1,nnn+1,irun)**2+ppion(2,nnn+1,irun)**2+ppion(3,nnn+1,irun)**2) - e1
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
     rpion(1, nnn+1, irun) = xfnl
     rpion(2, nnn+1, irun) = yfnl
     rpion(3, nnn+1, irun) = zfnl
     tfdpi(nnn+1, irun) = tfnl
  End If
  Return
200 Format (A30, 2(1X,E10.4))
210 Format (I6, 5(1X,F8.3))
220 Format (A2, I5, 5(1X,F8.3))
End Subroutine dkine2
