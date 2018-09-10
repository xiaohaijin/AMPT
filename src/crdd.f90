Subroutine crdd(irun, px, py, pz, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (xmd=1.8756, npdmax=10000)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /input/nstar, ndirct, dir
  Common /nn/nnn
  Common /bg/betax, betay, betaz, gamma
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Dimension ppd(3, npdmax), lbpd(npdmax)
  Save
  n12 = 0
  m12 = 0
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  c2 = pz/pr
  If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .And. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Then
        Goto 108
     Else
        Return
     End If
  End If
  If (x1<=signn/sig) Then
     as = (3.65*(srt-1.8766))**6
     a = 6.0*as/(1.0+as)
     ta = -2.0*pr**2
     x = ranart(nseed)
     t1 = sngl(dlog(dble(1.-x)*dexp(dble(a)*dble(ta))+dble(x)))/a
     c1 = 1.0 - t1/ta
     t1 = 2.0*pi*ranart(nseed)
     iblock = 20
     Goto 107
  Else
     If (srt<2.15) Return
     Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
     akp = 0.498
     ak0 = 0.498
     ana = 0.938
     ada = 1.232
     al = 1.1157
     as = 1.1197
     xsk1 = 0
     xsk2 = 0
     xsk3 = 0
     xsk4 = 0
     xsk5 = 0
     t1nlk = ana + al + akp
     If (srt<=t1nlk) Goto 222
     xsk1 = 1.5*pplpk(srt)
     t1dlk = ada + al + akp
     t2dlk = ada + al - akp
     If (srt<=t1dlk) Goto 222
     es = srt
     pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
     pmdlk = sqrt(pmdlk2)
     xsk3 = 1.5*pplpk(srt)
     t1nsk = ana + as + akp
     t2nsk = ana + as - akp
     If (srt<=t1nsk) Goto 222
     pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
     pmnsk = sqrt(pmnsk2)
     xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
     t1dsk = ada + as + akp
     t2dsk = ada + as - akp
     If (srt<=t1dsk) Goto 222
     pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
     pmdsk = sqrt(pmdsk2)
     xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
     If (srt<=(2.*amn+aphi)) Goto 222
     xsk5 = 0.0001
222  sigk = xsk1 + xsk2 + xsk3 + xsk4
     xsk1 = 2.0*xsk1
     xsk2 = 2.0*xsk2
     xsk3 = 2.0*xsk3
     xsk4 = 2.0*xsk4
     sigk = 2.0*sigk + xsk5
     s2d = reab2d(i1, i2, srt)
     s2d = 0.
     If (((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=12)) .Or. ((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=6)) .Or. ((iabs(lb(i2))>=12) .And. (iabs(lb(i1))>=6))) Then
        signd = sigk + s2d
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signd+signn+sdprod)/sig) Return
        If ((sigk+sdprod)/sig>=ranart(nseed)) Goto 306
        Goto 1012
     End If
     idd = iabs(lb(i1)*lb(i2))
     If ((idd==63) .Or. (idd==64) .Or. (idd==48) .Or. (idd==49) .Or. (idd==11*11) .Or. (idd==10*10) .Or. (idd==88) .Or. (idd==66) .Or. (idd==90) .Or. (idd==70)) Then
        signd = x1535 + sigk + s2d
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+signd+sdprod)/sig) Return
        If (sigk/signd>ranart(nseed)) Goto 306
        If (s2d/(x1535+s2d)>ranart(nseed)) Goto 1012
        If (idd==63) n12 = 17
        If (idd==64) n12 = 20
        If (idd==48) n12 = 23
        If (idd==49) n12 = 24
        If (idd==121) n12 = 25
        If (idd==100) n12 = 26
        If (idd==88) n12 = 29
        If (idd==66) n12 = 31
        If (idd==90) n12 = 32
        If (idd==70) n12 = 35
        Goto 1011
     End If
     If ((idd==110) .Or. (idd==77) .Or. (idd==80)) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+x1535+sigk+s2d+sdprod)/sig) Return
        If (sigk/(x1535+sigk+s2d)>ranart(nseed)) Goto 306
        If (s2d/(x1535+s2d)>ranart(nseed)) Goto 1012
        If (idd==77) n12 = 30
        If ((idd==77) .And. (ranart(nseed)<=0.5)) n12 = 36
        If (idd==80) n12 = 34
        If ((idd==80) .And. (ranart(nseed)<=0.5)) n12 = 35
        If (idd==110) n12 = 27
        If ((idd==110) .And. (ranart(nseed)<=0.5)) n12 = 28
        Goto 1011
     End If
     If ((idd==54) .Or. (idd==56)) Then
        sig2 = (3./4.)*sigma(srt, 2, 0, 1)
        signd = 2.*(sig2+x1535) + sigk + s2d
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+signd+sdprod)/sig) Return
        If (sigk/signd>ranart(nseed)) Goto 306
        If (s2d/(2.*(sig2+x1535)+s2d)>ranart(nseed)) Goto 1012
        If (ranart(nseed)<x1535/(sig2+x1535)) Then
           If (idd==54) n12 = 18
           If ((idd==54) .And. (ranart(nseed)<=0.5)) n12 = 19
           If (idd==56) n12 = 21
           If ((idd==56) .And. (ranart(nseed)<=0.5)) n12 = 22
        Else
           If (idd==54) n12 = 13
           If ((idd==54) .And. (ranart(nseed)<=0.5)) n12 = 14
           If (idd==56) n12 = 15
           If ((idd==56) .And. (ranart(nseed)<=0.5)) n12 = 16
        End If
     End If
1011 Continue
     iblock = 5
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If ((n12>=13) .And. (n12<=16)) Then
        If (dmax<1.44) Then
           fm = fns(dmax, srt, 0.)
        Else
           xdmass = 1.44
           fm = fns(xdmass, srt, 1.)
        End If
        If (fm==0.) fm = 1.E-09
        ntry2 = 0
11      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry2 = ntry2 + 1
        If ((ranart(nseed)>fns(dm,srt,1.)/fm) .And. (ntry2<=10)) Goto 11
        If (dm>2.14) Goto 11
        Goto 13
     End If
     If ((n12>=17) .And. (n12<=36)) Then
        If (dmax<1.535) Then
           fm = fd5(dmax, srt, 0.)
        Else
           xdmass = 1.535
           fm = fd5(xdmass, srt, 1.)
        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
12      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fd5(dm,srt,1.)/fm) .And. (ntry1<=10)) Goto 12
        If (dm>1.84) Goto 12
     End If
13   Continue
     If (n12==13) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 11
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 11
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==14) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 10
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 10
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==15) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 11
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 11
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==16) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 10
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 10
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==17) Then
        lb(i2) = 13
        e(i2) = dm
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (n12==18) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==19) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==20) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==21) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==22) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==23) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==24) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (n12==25) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (n12==26) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (n12==27) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==28) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==27) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==29) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==30) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==31) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==32) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==33) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==34) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     If (n12==35) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (n12==36) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
1012 Continue
     iblock = 55
     lb1 = lb(i1)
     lb2 = lb(i2)
     ich = iabs(lb1*lb2)
     If (ich==9*6) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (ich==8*7) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (ich==9*7) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (ich==8*8) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (ich==8*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (ich==6*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (ich==11*11 .Or. ich==13*13 .Or. ich==11*13) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (ich==10*10 .Or. ich==12*12 .Or. ich==10*12) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (ich==10*11 .Or. ich==12*13 .Or. ich==10*13 .Or. ich==11*12) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (ich==11*8 .Or. ich==13*8) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (ich==11*7 .Or. ich==13*7) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     If (ich==11*6 .Or. ich==13*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (ich==10*9 .Or. ich==12*9) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     If (ich==10*7 .Or. ich==12*7) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     If (ich==10*8 .Or. ich==12*8) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     lb(i1) = 1
     e(i1) = amp
     lb(i2) = 2
     e(i2) = amn
200  em1 = e(i1)
     em2 = e(i2)
     pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
     If (pr2<=0.) pr2 = 1.E-09
     pr = sqrt(pr2)/(2.*srt)
     If (srt<=2.14) c1 = 1.0 - 2.0*ranart(nseed)
     If (srt>2.14 .And. srt<=2.4) c1 = ang(srt, iseed)
     If (srt>2.4) Then
        xptr = 0.33*pr
        cc1 = ptr(xptr, iseed)
        scheck = pr**2 - cc1**2
        If (scheck<0) Then
           Write (99, *) 'scheck7: ', scheck
           scheck = 0.
        End If
        c1 = sqrt(scheck)/pr
     End If
     t1 = 2.0*pi*ranart(nseed)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(i1) = -lb(i1)
        lb(i2) = -lb(i2)
     End If
  End If
107 scheck = 1.0 - c1**2
  If (scheck<0) Then
     Write (99, *) 'scheck8: ', scheck
     scheck = 0.
  End If
  s1 = sqrt(scheck)
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck9: ', scheck
     scheck = 0.
  End If
  s2 = sqrt(scheck)
  ct1 = cos(t1)
  st1 = sin(t1)
  ct2 = cos(t2)
  st2 = sin(t2)
  pz = pr*(c1*c2-s1*s2*ct1)
  ss = c2*s1*ct1 + s2*c1
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  Return
306 Continue
  If (xsk5/sigk>ranart(nseed)) Then
     pz1 = p(3, i1)
     pz2 = p(3, i2)
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 1 + int(2*ranart(nseed))
     nnn = nnn + 1
     lpion(nnn, irun) = 29
     epion(nnn, irun) = aphi
     iblock = 222
     Goto 208
  End If
  iblock = 10
  If (ianti==1) iblock = -10
  pz1 = p(3, i1)
  pz2 = p(3, i2)
  nnn = nnn + 1
  lpion(nnn, irun) = 23
  epion(nnn, irun) = aka
  If (srt<=2.63) Then
     ic = 1
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 14
     Goto 208
  End If
  If (srt<=2.74 .And. srt>2.63) Then
     If (xsk1/(xsk1+xsk2)>ranart(nseed)) Then
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
     Else
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 15 + int(3*ranart(nseed))
        ic = 2
     End If
     Goto 208
  End If
  If (srt<=2.77 .And. srt>2.74) Then
     If (xsk1/(xsk1+xsk2+xsk3)>ranart(nseed)) Then
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk2/(xsk2+xsk3)>ranart(nseed)) Then
           ic = 2
           lb(i1) = 1 + int(2*ranart(nseed))
           lb(i2) = 15 + int(3*ranart(nseed))
        Else
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
        End If
        Goto 208
     End If
  End If
  If (srt>2.77) Then
     If (xsk1/(xsk1+xsk2+xsk3+xsk4)>ranart(nseed)) Then
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk3/(xsk2+xsk3+xsk4)>ranart(nseed)) Then
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
           Goto 208
        Else
           If (xsk2/(xsk2+xsk4)>ranart(nseed)) Then
              lb(i1) = 1 + int(2*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
              ic = 2
           Else
              ic = 4
              lb(i1) = 6 + int(4*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
           End If
           Goto 208
        End If
     End If
  End If
208 Continue
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==23) lpion(nnn, irun) = 21
  End If
  lbi1 = lb(i1)
  lbi2 = lb(i2)
  ntry1 = 0
129 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 129
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  epcm = sqrt(aka**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lbi1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lbi2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  lb1 = lb(i1)
  lb2 = lb(i2)
  am1 = em1
  am2 = em2
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  Return
108 Continue
  If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
     ndloop = npertd
  Else If (idpert==2 .And. npertd>=1) Then
     ndloop = npertd + 1
  Else
     ndloop = 1
  End If
  dprob1 = sdprod/sig/float(npertd)
  Do idloop = 1, ndloop
     Call bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
     Call rotate(px, py, pz, pxd, pyd, pzd)
     xmass = xmd
     e1dcm = sqrt(xmass**2+pxd**2+pyd**2+pzd**2)
     p1dbeta = pxd*betax + pyd*betay + pzd*betaz
     transf = gamma*(gamma*p1dbeta/(gamma+1.)+e1dcm)
     pxi1 = betax*transf + pxd
     pyi1 = betay*transf + pyd
     pzi1 = betaz*transf + pzd
     If (ianti==0) Then
        lbd = 42
     Else
        lbd = -42
     End If
     If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
        nnn = nnn + 1
        ppion(1, nnn, irun) = pxi1
        ppion(2, nnn, irun) = pyi1
        ppion(3, nnn, irun) = pzi1
        epion(nnn, irun) = xmd
        lpion(nnn, irun) = lbd
        rpion(1, nnn, irun) = r(1, i1)
        rpion(2, nnn, irun) = r(2, i1)
        rpion(3, nnn, irun) = r(3, i1)
        dppion(nnn, irun) = sdprod/sig/float(npertd)
     Else If (idpert==2 .And. idloop<=npertd) Then
        ppd(1, idloop) = pxi1
        ppd(2, idloop) = pyi1
        ppd(3, idloop) = pzi1
        lbpd(idloop) = lbd
     Else
        e(i1) = xmm
        e2picm = sqrt(xmm**2+pxd**2+pyd**2+pzd**2)
        p2pibeta = -pxd*betax - pyd*betay - pzd*betaz
        transf = gamma*(gamma*p2pibeta/(gamma+1.)+e2picm)
        pxi2 = betax*transf - pxd
        pyi2 = betay*transf - pyd
        pzi2 = betaz*transf - pzd
        p(1, i1) = pxi2
        p(2, i1) = pyi2
        p(3, i1) = pzi2
        lb(i1) = lbm
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        id(i1) = 2
        id1 = id(i1)
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        lb1 = lb(i1)
        p(1, i2) = pxi1
        p(2, i2) = pyi1
        p(3, i2) = pzi1
        lb(i2) = lbd
        lb2 = lb(i2)
        e(i2) = xmd
        eti2 = e(i2)
        id(i2) = 2
        If (idpert==2 .And. idloop==ndloop) Then
           Do ipertd = 1, npertd
              nnn = nnn + 1
              ppion(1, nnn, irun) = ppd(1, ipertd)
              ppion(2, nnn, irun) = ppd(2, ipertd)
              ppion(3, nnn, irun) = ppd(3, ipertd)
              epion(nnn, irun) = xmd
              lpion(nnn, irun) = lbpd(ipertd)
              rpion(1, nnn, irun) = r(1, i1)
              rpion(2, nnn, irun) = r(2, i1)
              rpion(3, nnn, irun) = r(3, i1)
              dppion(nnn, irun) = 1./float(npertd)
           End Do
        End If
     End If
  End Do
  iblock = 501
  Return
End Subroutine crdd
