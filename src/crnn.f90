Subroutine crnn(irun, px, py, pz, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383, aphi=1.020)
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
  Common /table/xarray(0:1000), earray(0:1000)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /rndf77/nseed
  Common /dpi/em2, lb2
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
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
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))==1 .Or. iabs(lb(i1))==2) .And. (iabs(lb(i2))==1 .Or. iabs(lb(i2))==2)) Then
        Goto 108
     Else
        Return
     End If
  End If
  If (x1<=(signn/sig)) Then
     as = (3.65*(srt-1.8766))**6
     a = 6.0*as/(1.0+as)
     ta = -2.0*pr**2
     x = ranart(nseed)
     t1 = sngl(dlog(dble(1.-x)*dexp(dble(a)*dble(ta))+dble(x)))/a
     c1 = 1.0 - t1/ta
     t1 = 2.0*pi*ranart(nseed)
     iblock = 1
     Goto 107
  Else
     If (srt<2.012) Return
     Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
     sig3 = 3.*(x3pi(srt)+x33pi(srt))
     sig4 = 4.*x2pi(srt)
     s4pi = x4pi(srt)
     srho = xrho(srt)
     somega = omega(srt)
     akp = 0.498
     ak0 = 0.498
     ana = 0.94
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
     lb1 = iabs(lb(i1))
     lb2 = iabs(lb(i2))
     If ((lb(i1)*lb(i2)==1) .Or. ((lb1<=17 .And. lb1>=14) .And. (lb2<=17 .And. lb2>=14)) .Or. ((lb1<=2) .And. (lb2<=17 .And. lb2>=14)) .Or. ((lb2<=2) .And. (lb1<=17 .And. lb1>=14))) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sig2 = 1.5*sigma(srt, 1, 1, 1)
        signd = sig1 + sig2 + sig3 + sig4 + x1535 + sigk + s4pi + srho + somega
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(sigk+x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 306
        If (ranart(nseed)<=s4pi/(x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 307
        If (ranart(nseed)<=srho/(x1535+sig4+sig2+sig1+srho+somega)) Goto 308
        If (ranart(nseed)<=somega/(x1535+sig4+sig2+sig1+somega)) Goto 309
        If (ranart(nseed)<=x1535/(sig1+sig2+sig4+x1535)) Then
           n12 = 9
        Else
           If (ranart(nseed)<=sig4/(sig1+sig2+sig4)) Then
              n12 = 66
              Goto 1012
           Else
              n12 = 3
              If (ranart(nseed)>sig1/(sig1+sig2)) n12 = 4
           End If
        End If
        Goto 1011
     End If
     If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sig2 = 1.5*sigma(srt, 1, 1, 1)
        signd = sig1 + sig2 + x1535 + sig3 + sig4 + sigk + s4pi + srho + somega
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(sigk+x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 306
        If (ranart(nseed)<=s4pi/(x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 307
        If (ranart(nseed)<=srho/(x1535+sig4+sig2+sig1+srho+somega)) Goto 308
        If (ranart(nseed)<=somega/(x1535+sig4+sig2+sig1+somega)) Goto 309
        If (ranart(nseed)<=x1535/(x1535+sig1+sig2+sig4)) Then
           n12 = 10
        Else
           If (ranart(nseed)<=sig4/(sig1+sig2+sig4)) Then
              n12 = 67
              Goto 1013
           Else
              n12 = 6
              If (ranart(nseed)>sig1/(sig1+sig2)) n12 = 5
           End If
        End If
        Goto 1011
     End If
     If (lb(i1)*lb(i2)==2) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        If (nstar==1) Then
           sig2 = (3./4.)*sigma(srt, 2, 0, 1)
        Else
           sig2 = 0.
        End If
        signd = 2.*(sig1+sig2+x1535) + sig3 + sig4 + sigk + s4pi + srho + somega
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(signd-sig3)) Goto 306
        If (ranart(nseed)<=s4pi/(signd-sig3-sigk)) Goto 307
        If (ranart(nseed)<=srho/(signd-sig3-sigk-s4pi)) Goto 308
        If (ranart(nseed)<=somega/(signd-sig3-sigk-s4pi-srho)) Goto 309
        If (ranart(nseed)<x1535/(sig1+sig2+x1535+0.5*sig4)) Then
           n12 = 11
           If (ranart(nseed)<=0.5) n12 = 12
        Else
           If (ranart(nseed)<=sig4/(sig4+2.*(sig1+sig2))) Then
              n12 = 68
              Goto 1014
           Else
              If (ranart(nseed)<=sig1/(sig1+sig2)) Then
                 n12 = 2
                 If (ranart(nseed)>=0.5) n12 = 1
              Else
                 n12 = 8
                 If (ranart(nseed)>=0.5) n12 = 7
              End If
           End If
        End If
     End If
1011 iblock = 2
     Continue
     dmax = srt - avmass - 0.005
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If (n12<7) Then
        If (dmax<1.232) Then
           fm = fde(dmax, srt, 0.)
        Else
           xdmass = 1.232
           fm = fde(xdmass, srt, 1.)
        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
10      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fde(dm,srt,1.)/fm) .And. (ntry1<=30)) Goto 10
        If (dm>1.47) Goto 10
        Goto 13
     End If
     If ((n12==7) .Or. (n12==8)) Then
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
     If (n12>=17) Then
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
        Goto 13
     End If
1012 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==66) Then
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           lb(i1) = 9
           lb(i2) = 7
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           lb(i1) = 8
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           lb(i1) = 9
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
        If (xfinal>0.75) Then
           lb(i1) = 8
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
     End If
1013 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==67) Then
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           lb(i1) = 7
           lb(i2) = 7
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           lb(i1) = 6
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           lb(i1) = 7
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
        If (xfinal>0.75) Then
           lb(i1) = 8
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
     End If
1014 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==68) Then
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           lb(i1) = 7
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           lb(i1) = 9
           lb(i2) = 6
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           lb(i1) = 7
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
        If (xfinal>0.75) Then
           lb(i1) = 8
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
        End If
     End If
13   Continue
     If (n12==1) Then
        If (iabs(lb(i1))==1) Then
           lb(i2) = 2
           lb(i1) = 8
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 8
           e(i2) = dm
        End If
        Goto 200
     End If
     If (n12==2) Then
        If (iabs(lb(i1))==2) Then
           lb(i2) = 1
           lb(i1) = 7
           e(i1) = dm
        Else
           lb(i1) = 1
           lb(i2) = 7
           e(i2) = dm
        End If
        Goto 200
     End If
     If (n12==3) Then
        lb(i1) = 9
        e(i1) = dm
        lb(i2) = 2
        e(i2) = amn
        Goto 200
     End If
     If (n12==4) Then
        lb(i2) = 1
        lb(i1) = 8
        e(i1) = dm
        Goto 200
     End If
     If (n12==5) Then
        lb(i2) = 2
        lb(i1) = 7
        e(i1) = dm
        Goto 200
     End If
     If (n12==6) Then
        lb(i1) = 6
        e(i1) = dm
        lb(i2) = 1
        e(i2) = amp
        Goto 200
     End If
     If (n12==7) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 200
     End If
     If (n12==8) Then
        If (iabs(lb(i1))==1) Then
           lb(i2) = 2
           lb(i1) = 11
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 11
           e(i2) = dm
        End If
        Goto 200
     End If
     If (n12==9) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           lb(i1) = 13
           e(i1) = dm
        Else
           lb(i1) = 1
           lb(i2) = 13
           e(i2) = dm
        End If
        Goto 200
     End If
     If (n12==10) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 2
           lb(i1) = 12
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 12
           e(i2) = dm
        End If
        Goto 200
     End If
     If (n12==11) Then
        If (iabs(lb(i1))==2) Then
           lb(i1) = 2
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 200
     End If
     If (n12==12) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           lb(i1) = 12
           e(i1) = dm
        End If
     End If
  End If
200 em1 = e(i1)
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
        Write (99, *) 'scheck2: ', scheck
        scheck = 0.
     End If
     c1 = sqrt(scheck)/pr
  End If
  t1 = 2.0*pi*ranart(nseed)
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If
  Goto 107
106 Continue
  ntry1 = 0
123 Call ddp2(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=40)) Goto 123
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.2) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 9
        lb(i2) = 7
        Goto 205
     End If
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 8
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir<=0.6) .And. (xdir>0.4)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir<=0.8) .And. (xdir>0.6)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 6
        Goto 205
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
  End If
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.2) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 6
        lb(i2) = 7
        Goto 205
     End If
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 6
        lb(i2) = 9
        Goto 205
     End If
     If ((xdir>0.4) .And. (xdir<=0.6)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir>0.6) .And. (xdir<=0.8)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 7
        lb(i2) = 7
        Goto 205
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
  End If
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.17) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 6
        lb(i2) = 9
        Goto 205
     End If
     If ((xdir<=0.34) .And. (xdir>0.17)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 9
        Goto 205
     End If
     If ((xdir>0.34) .And. (xdir<=0.51)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir>0.51) .And. (xdir<=0.68)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 8
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir>0.68) .And. (xdir<=0.85)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
     If (xdir>0.85) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 7
     End If
  End If
205 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==3) Then
        lpion(nnn, irun) = 5
     Else If (lpion(nnn,irun)==5) Then
        lpion(nnn, irun) = 3
     End If
  End If
  lb1 = lb(i1)
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 4
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  Goto 90005
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
  Goto 90005
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
  iblock = 9
  If (ianti==1) iblock = -9
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
  ntry1 = 0
127 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 127
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
  lbi1 = lb(i1)
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lbi2 = lb(i2)
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
  Goto 90005
307 Continue
  ntry1 = 0
125 Call ddrho(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, amrho, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 125
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  arho = amrho
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.2) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 7
        Goto 2051
     End If
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 8
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir<=0.6) .And. (xdir>0.4)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir<=0.8) .And. (xdir>0.6)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 6
        Goto 2051
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
  End If
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.2) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 7
        Goto 2051
     End If
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 9
        Goto 2051
     End If
     If ((xdir>0.4) .And. (xdir<=0.6)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir>0.6) .And. (xdir<=0.8)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 7
        Goto 2051
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
  End If
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.17) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 9
        Goto 2051
     End If
     If ((xdir<=0.34) .And. (xdir>0.17)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 9
        Goto 2051
     End If
     If ((xdir>0.34) .And. (xdir<=0.51)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir>0.51) .And. (xdir<=0.68)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 8
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir>0.68) .And. (xdir<=0.85)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
     If (xdir>0.85) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 7
     End If
  End If
2051 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==25) Then
        lpion(nnn, irun) = 27
     Else If (lpion(nnn,irun)==27) Then
        lpion(nnn, irun) = 25
     End If
  End If
  lb1 = lb(i1)
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 44
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  Goto 90005
308 Continue
  ntry1 = 0
126 Call pprho(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, amrho, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 126
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  arho = amrho
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.5) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 1
        Goto 2052
     Else
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
     End If
  End If
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.5) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 2
        lb(i2) = 2
        Goto 2052
     Else
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
     End If
  End If
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.33) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
     Else If ((xdir<=0.67) .And. (xdir>0.34)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 1
        Goto 2052
     Else
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 2
        lb(i2) = 2
        Goto 2052
     End If
  End If
2052 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==25) Then
        lpion(nnn, irun) = 27
     Else If (lpion(nnn,irun)==27) Then
        lpion(nnn, irun) = 25
     End If
  End If
  lb1 = lb(i1)
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 45
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  Goto 90005
309 Continue
  ntry1 = 0
138 Call ppomga(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 138
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  aomega = 0.782
  If (lb(i1)*lb(i2)==1) Then
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 1
     lb(i2) = 1
     Goto 2053
  End If
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 2
     lb(i2) = 2
     Goto 2053
  End If
  If (lb(i1)*lb(i2)==2) Then
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 1
     lb(i2) = 2
     Goto 2053
  End If
2053 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If
  lb1 = lb(i1)
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 46
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  Goto 90005
90005 Continue
  Return
107 If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If
  s1 = 1.0 - c1**2
  If (s1<=0) s1 = 0
  s1 = sqrt(s1)
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck3: ', scheck
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
End Subroutine crnn
