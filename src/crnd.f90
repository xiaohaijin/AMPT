Subroutine crnd(irun, px, py, pz, srt, i1, i2, iblock, signn, sig, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
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
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))==1 .Or. iabs(lb(i1))==2) .And. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Then
        Goto 108
     Else If ((iabs(lb(i2))==1 .Or. iabs(lb(i2))==2) .And. (iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13)) Then
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
     iblock = 1
     Goto 107
  Else
     If (srt<2.04) Return
     If (((iabs(lb(i1))==2 .Or. iabs(lb(i2))==2) .And. (lb(i1)*lb(i2))==20) .Or. (lb(i1)*lb(i2))==13) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
     End If
     prf = sqrt(0.25*srt**2-avmass**2)
     If (em1>1.) Then
        deltam = em1
     Else
        deltam = em2
     End If
     renom = deltam*prf**2/denom(srt, 1.)/pr
     renomn = deltam*prf**2/denom(srt, 2.)/pr
     renom1 = deltam*prf**2/denom(srt, -1.)/pr
     If ((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) renom = 0.
     If ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) renom = 0.
     If ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) renom = 0.
     If ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9)) renom = 0.
     Call m1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
     x1440 = (3./4.)*sigma(srt, 2, 0, 1)
     If (((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) .Or. ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) .Or. ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) .Or. ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9))) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If ((sigk+signn+sdprod)/sig>=x1) Goto 306
     End If
     If (lb(i1)*lb(i2)==18 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 3
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 37
           Else
              Return
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==6 .And. ((iabs(lb(i1))==1) .Or. (iabs(lb(i2))==1))) Then
        signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 6
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 47
           Else
              Return
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==8 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
        signd = 1.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 4
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 39
           Else
              m12 = 40
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==14 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = 1.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 5
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 48
           Else
              m12 = 49
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==16 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        sigdn = 0.5*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+2.*x1440+2.*x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+2*x1440+2*x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+2.*x1440+2.*x1535)) Then
           m12 = 1
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 41
              If (ranart(nseed)<=0.5) m12 = 43
           Else
              m12 = 42
              If (ranart(nseed)<=0.5) m12 = 44
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==7) Then
        signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        sigdn = 0.5*signd*renom
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+2.*x1440+2.*x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+2*x1440+2*x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+2.*x1440+2.*x1535)) Then
           m12 = 2
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 50
              If (ranart(nseed)<=0.5) m12 = 51
           Else
              m12 = 52
              If (ranart(nseed)<=0.5) m12 = 53
           End If
           Goto 204
        End If
     End If
     If (lb(i1)*lb(i2)==10 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
        signd = (3./4.)*sigma(srt, 2, 0, 1)
        sigdn = signd*renomn
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1535)) Then
           m12 = 7
           Goto 206
        Else
           m12 = 54
           If (ranart(nseed)<=0.5) m12 = 55
        End If
        Goto 204
     End If
     If (lb(i1)*lb(i2)==22 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = (3./4.)*sigma(srt, 2, 0, 1)
        sigdn = signd*renomn
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+x1535+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1535)) Then
           m12 = 8
           Goto 206
        Else
           m12 = 56
           If (ranart(nseed)<=0.5) m12 = 57
        End If
        Goto 204
     End If
     If ((iabs(lb(i1))==12) .Or. (iabs(lb(i1))==13) .Or. (iabs(lb(i2))==12) .Or. (iabs(lb(i2))==13)) Then
        signd = x1535
        sigdn = signd*renom1
        If (x1<=((signn+sdprod)/sig)) Goto 108
        If (x1>(signn+sigdn+sigk+sdprod)/sig) Return
        If (sigk/(sigk+sigdn)>ranart(nseed)) Goto 306
        If (lb(i1)*lb(i2)==24) m12 = 10
        If (lb(i1)*lb(i2)==12) m12 = 12
        If (lb(i1)*lb(i2)==26) m12 = 11
        If (lb(i1)*lb(i2)==13) m12 = 9
        Goto 206
     End If
204  Continue
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If ((m12==37) .Or. (m12==39) .Or. (m12==41) .Or. (m12==43) .Or. (m12==46) .Or. (m12==48) .Or. (m12==50) .Or. (m12==51)) Then
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
     Else
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
     prf = 0.
     pf2 = ((srt**2-dm**2+avmass**2)/(2.*srt))**2 - avmass**2
     If (pf2>0.) prf = sqrt(pf2)
     If (m12==37) Then
        If (iabs(lb(i1))==9) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==38) Then
        If (iabs(lb(i1))==9) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==39) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==40) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==41) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==42) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==43) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==44) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==46) Then
        If (iabs(lb(i1))==6) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==47) Then
        If (iabs(lb(i1))==6) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==48) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==49) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==50) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==51) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==52) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==53) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==54) Then
        If (iabs(lb(i1))==10) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==55) Then
        If (iabs(lb(i1))==10) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==56) Then
        If (iabs(lb(i1))==11) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     If (m12==57) Then
        If (iabs(lb(i1))==11) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
     End If
     Goto 207
206  If (m12==1) Then
        If (iabs(lb(i1))==8) Then
           lb(i2) = 2
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 2
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 207
     End If
     If (m12==2) Then
        If (iabs(lb(i1))==7) Then
           lb(i2) = 1
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 207
     End If
     If (m12==3) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     If (m12==4) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     If (m12==5) Then
        lb(i1) = 2
        lb(i2) = 2
        e(i1) = amn
        e(i2) = amn
        Goto 207
     End If
     If (m12==6) Then
        lb(i1) = 2
        lb(i2) = 2
        e(i1) = amn
        e(i2) = amn
        Goto 207
     End If
     If (m12==7) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 2
           e(i1) = amp
           e(i2) = amn
        Else
           lb(i1) = 2
           lb(i2) = 1
           e(i1) = amn
           e(i2) = amp
        End If
        Goto 207
     End If
     If (m12==8) Then
        If (iabs(lb(i1))==2) Then
           lb(i1) = 2
           lb(i2) = 1
           e(i1) = amn
           e(i2) = amp
        Else
           lb(i1) = 1
           lb(i2) = 2
           e(i1) = amp
           e(i2) = amn
        End If
        Goto 207
     End If
     If (m12==9) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     If (m12==12) Then
        lb(i1) = 2
        lb(i2) = 1
        e(i1) = amn
        e(i2) = amp
        Goto 207
     End If
     If (m12==11) Then
        lb(i1) = 2
        lb(i2) = 1
        e(i1) = amn
        e(i2) = amp
        Goto 207
     End If
     If (m12==12) Then
        lb(i1) = 1
        lb(i2) = 2
        e(i1) = amp
        e(i2) = amn
     End If
207  pr = prf
     c1 = 1.0 - 2.0*ranart(nseed)
     If (srt<=2.14) c1 = 1.0 - 2.0*ranart(nseed)
     If (srt>2.14 .And. srt<=2.4) c1 = ang(srt, iseed)
     If (srt>2.4) Then
        xptr = 0.33*pr
        cc1 = ptr(xptr, iseed)
        scheck = pr**2 - cc1**2
        If (scheck<0) Then
           Write (99, *) 'scheck4: ', scheck
           scheck = 0.
        End If
        c1 = sqrt(scheck)/pr
     End If
     t1 = 2.0*pi*ranart(nseed)
     iblock = 3
  End If
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If
107 If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If
  scheck = 1.0 - c1**2
  If (scheck<0) Then
     Write (99, *) 'scheck5: ', scheck
     scheck = 0.
  End If
  s1 = sqrt(scheck)
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck6: ', scheck
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
  iblock = 11
  If (ianti==1) iblock = -11
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
128 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 128
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
  If (lpion(nnn,irun)/=29) iblock = 11
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
End Subroutine crnd
