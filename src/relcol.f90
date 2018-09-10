Subroutine relcol(lcoll, lbloc, lcnne, ldd, lpp, lppk, lpn, lpd, lrho, lomega, lkn, lnnk, lddk, lndk, lcnnd, lcndn, ldirt, ldecay, lres, ldou, lddrho, lnnrho, lnnom, nt, ntmax, sp, akaon, sk)
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, aks=0.895)
  Parameter (aa1=1.26, aphi=1.02, ap1=0.13496)
  Parameter (maxx=20, maxz=24)
  Parameter (rrkk=0.6, prkk=0.3, srhoks=5., esbin=0.04)
  Dimension massrn(0:maxr), rt(3, maxstr), pt(3, maxstr), et(maxstr)
  Dimension lt(maxstr), prot(maxstr)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Common /hh/proper(maxstr)
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /input/nstar, ndirct, dir
  Common /nn/nnn
  Common /rr/massr(0:maxr)
  Common /ss/inout(20)
  Common /bg/betax, betay, betaz, gamma
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /pe/propi(maxstr, maxr)
  Common /kkk/tkaon(7), ekaon(7, 0:2000)
  Common /kaon/ak(3, 50, 36), speck(50, 36, 7), mf
  Common /table/xarray(0:1000), earray(0:1000)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  Common /lastt/itimeh, bimp
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  Common /resdcy/nsav, iksdcy
  Common /rndf77/nseed
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Dimension ftpisv(maxstr, maxr), fttemp(maxstr)
  Common /dpi/em2, lb2
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Dimension dptemp(maxstr)
  Common /para8/idpert, npertd, idxsec
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./
  Call inidcy
  resona = 5.
  nodelt = 0
  sumsrt = 0.
  lcoll = 0
  lbloc = 0
  lcnne = 0
  ldd = 0
  lpp = 0
  lpd = 0
  lpdr = 0
  lrho = 0
  lrhor = 0
  lomega = 0
  lomgar = 0
  lpn = 0
  lkn = 0
  lnnk = 0
  lddk = 0
  lndk = 0
  lppk = 0
  lcnnd = 0
  lcndn = 0
  ldirt = 0
  ldecay = 0
  lres = 0
  ldou = 0
  lddrho = 0
  lnnrho = 0
  lnnom = 0
  msum = 0
  massrn(0) = 0
  Do il = 1, 5
     tkaon(il) = 0
     Do is = 1, 2000
        ekaon(il, is) = 0
     End Do
  End Do
  Do i = 1, num
     Do j = 1, maxstr
        propi(j, i) = 1.
     End Do
  End Do
  Do i = 1, maxstr
     fttemp(i) = 0.
     Do irun = 1, maxr
        ftpisv(i, irun) = 0.
     End Do
  End Do
  sp = 0
  akaon = 0
  sk = 0
  mass = 0
  Do irun = 1, num
     nnn = 0
     msum = msum + massr(irun-1)
     j10 = 2
     If (nt==ntmax) j10 = 1
     Do j1 = j10, massr(irun)
        i1 = j1 + msum
        If (e(i1)==0.) Goto 798
        If (lb(i1)<-45 .Or. lb(i1)>45) Goto 798
        x1 = r(1, i1)
        y1 = r(2, i1)
        z1 = r(3, i1)
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        am1 = em1
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        id1 = id(i1)
        lb1 = lb(i1)
        If (nt==ntmax .And. (lb1==21 .Or. lb1==23)) Then
           pk0 = ranart(nseed)
           If (pk0<0.25) Then
              lb(i1) = 22
           Else If (pk0<0.50) Then
              lb(i1) = 24
           End If
           lb1 = lb(i1)
        End If
        If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13) .Or. (iksdcy==1 .And. lb1==24) .Or. iabs(lb1)==16 .Or. (ipi0dcy==1 .And. nt==ntmax .And. lb1==4)) Then
           Continue
        Else
           Goto 1
        End If
        If (lb1>=25 .And. lb1<=27) Then
           wid = 0.151
        Else If (lb1==28) Then
           wid = 0.00841
        Else If (lb1==29) Then
           wid = 0.00443
        Else If (iabs(lb1)==30) Then
           wid = 0.051
        Else If (lb1==0) Then
           wid = 1.18E-6
        Else If (iksdcy==1 .And. lb1==24) Then
           wid = 7.36E-15
        Else If (iabs(lb1)==16) Then
           wid = 8.87E-6
        Else If (lb1==32) Then
           Call wida1(em1, rhomp, wid, iseed)
        Else If (iabs(lb1)>=6 .And. iabs(lb1)<=9) Then
           wid = width(em1)
        Else If ((iabs(lb1)==10) .Or. (iabs(lb1)==11)) Then
           wid = w1440(em1)
        Else If ((iabs(lb1)==12) .Or. (iabs(lb1)==13)) Then
           wid = w1535(em1)
        Else If (ipi0dcy==1 .And. nt==ntmax .And. lb1==4) Then
           wid = 7.85E-9
        End If
        If (nt==ntmax) Then
           pdecay = 1.1
           If (iphidcy==0 .And. iabs(lb1)==29) pdecay = 0.
        Else
           t0 = 0.19733/wid
           gfactr = e1/em1
           t0 = t0*gfactr
           If (t0>0.) Then
              pdecay = 1. - exp(-dt/t0)
           Else
              pdecay = 0.
           End If
        End If
        xdecay = ranart(nseed)
        If (xdecay<pdecay) Then
           idecay = irun
           tfnl = nt*dt
           If (nt==ntmax .And. ftsv(i1)>((ntmax-1)*dt)) tfnl = ftsv(i1)
           xfnl = x1
           yfnl = y1
           zfnl = z1
           If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=9) .Or. (iksdcy==1 .And. lb1==24) .Or. iabs(lb1)==16 .Or. (ipi0dcy==1 .And. nt==ntmax .And. lb1==4)) Then
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              dpertp(i1) = dp1n
              If (nt==ntmax) Then
                 r(1, i1) = xfnl
                 r(2, i1) = yfnl
                 r(3, i1) = zfnl
                 tfdcy(i1) = tfnl
              End If
              If (iabs(lb1)>=6 .And. iabs(lb1)<=9) Then
                 ldecay = ldecay + 1
              End If
           Else If (iabs(lb1)==10 .Or. iabs(lb1)==11) Then
              nnn = nnn + 1
              ldecay = ldecay + 1
              pnstar = 1.
              If (e(i1)>1.22) pnstar = 0.6
              If (ranart(nseed)<=pnstar) Then
                 Call decay(idecay, i1, nnn, iseed, wid, nt)
              Else
                 Call decay2(idecay, i1, nnn, iseed, wid, nt)
                 nnn = nnn + 1
              End If
           Else If (iabs(lb1)==12 .Or. iabs(lb1)==13) Then
              nnn = nnn + 1
              Call decay(idecay, i1, nnn, iseed, wid, nt)
              ldecay = ldecay + 1
           End If
           If (nt==ntmax) Then
              If (lb(i1)==25 .Or. lb(i1)==26 .Or. lb(i1)==27) Then
                 wid = 0.151
              Else If (lb(i1)==0) Then
                 wid = 1.18E-6
              Else If (lb(i1)==24 .And. iksdcy==1) Then
                 wid = 7.36E-15
              Else If (ipi0dcy==1 .And. lb(i1)==4) Then
                 wid = 7.85E-9
              Else
                 Goto 9000
              End If
              lb1 = lb(i1)
              px1 = p(1, i1)
              py1 = p(2, i1)
              pz1 = p(3, i1)
              em1 = e(i1)
              e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              r(1, i1) = xfnl
              r(2, i1) = yfnl
              r(3, i1) = zfnl
              tfdcy(i1) = tfnl
              dpertp(i1) = dp1n
           End If
           If (nt==ntmax .And. ipi0dcy==1 .And. lb(i1)==4) Then
              wid = 7.85E-9
              lb1 = lb(i1)
              px1 = p(1, i1)
              py1 = p(2, i1)
              pz1 = p(3, i1)
              em1 = e(i1)
              e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              r(1, i1) = xfnl
              r(2, i1) = yfnl
              r(3, i1) = zfnl
              tfdcy(i1) = tfnl
              dpertp(i1) = dp1n
           End If
9000       Goto 798
        End If
1       If (nt==ntmax) Goto 798
        x1 = r(1, i1)
        y1 = r(2, i1)
        z1 = r(3, i1)
        Do j2 = 1, j1 - 1
           i2 = j2 + msum
           If (e(i2)==0.) Goto 600
           If (e(i1)==0.) Goto 800
           If (lb(i2)<-45 .Or. lb(i2)>45) Goto 600
           x2 = r(1, i2)
           y2 = r(2, i2)
           z2 = r(3, i2)
           dr0max = 5.
           ilb1 = iabs(lb(i1))
           ilb2 = iabs(lb(i2))
           If (ilb1==42 .Or. ilb2==42) Then
              If ((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) Then
                 If ((lb(i1)*lb(i2))>0) dr0max = 10.
              End If
           End If
           If (((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)>dr0max**2) Goto 600
           If (id(i1)*id(i2)==iavoid) Goto 400
           id1 = id(i1)
           id2 = id(i2)
           ix1 = nint(x1/dx)
           iy1 = nint(y1/dy)
           iz1 = nint(z1/dz)
           px1 = p(1, i1)
           py1 = p(2, i1)
           pz1 = p(3, i1)
           em1 = e(i1)
           am1 = em1
           lb1 = lb(i1)
           e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
           ipx1 = nint(px1/dpx)
           ipy1 = nint(py1/dpy)
           ipz1 = nint(pz1/dpz)
           lb2 = lb(i2)
           px2 = p(1, i2)
           py2 = p(2, i2)
           pz2 = p(3, i2)
           em2 = e(i2)
           am2 = em2
           lb1i = lb(i1)
           lb2i = lb(i2)
           px1i = p(1, i1)
           py1i = p(2, i1)
           pz1i = p(3, i1)
           em1i = e(i1)
           px2i = p(1, i2)
           py2i = p(2, i2)
           pz2i = p(3, i2)
           em2i = e(i2)
           eini = sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2) + sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
           pxini = p(1, i1) + p(1, i2)
           pyini = p(2, i1) + p(2, i2)
           pzini = p(3, i1) + p(3, i2)
           nnnini = nnn
           iblock = 0
           deltr0 = 3.
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb1)>=30 .And. iabs(lb1)<=45)) deltr0 = 5.0
           If ((iabs(lb2)>=14 .And. iabs(lb2)<=17) .Or. (iabs(lb2)>=30 .And. iabs(lb2)<=45)) deltr0 = 5.0
           If (lb1==28 .And. lb2==28) deltr0 = 4.84
           If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
              e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
              spipi = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
              If (spipi>=(4*0.77**2)) deltr0 = 3.5
           End If
           If (lb1==23 .And. (lb2>=14 .And. lb2<=17)) Goto 3699
           If (lb2==23 .And. (lb1>=14 .And. lb1<=17)) Goto 3699
           If (lb1==21 .And. lb2==23) Goto 3699
           If (lb2==21 .And. lb1==23) Goto 3699
           If (lb1==30 .And. lb2==21) Goto 3699
           If (lb2==30 .And. lb1==21) Goto 3699
           If (lb1==-30 .And. lb2==23) Goto 3699
           If (lb2==-30 .And. lb1==23) Goto 3699
           If (lb1==-30 .And. lb2==30) Goto 3699
           If (lb2==-30 .And. lb1==30) Goto 3699
           If (lb1==21 .Or. lb1==23) Then
              If (lb2==0 .Or. (lb2>=25 .And. lb2<=28)) Then
                 Goto 3699
              End If
           Else If (lb2==21 .Or. lb2==23) Then
              If (lb1==0 .Or. (lb1>=25 .And. lb1<=28)) Then
                 Goto 3699
              End If
           End If
           If (iabs(lb1)==30 .And. (lb2==0 .Or. (lb2>=25 .And. lb2<=28) .Or. (lb2>=3 .And. lb2<=5))) Then
              Goto 3699
           Else If (iabs(lb2)==30 .And. (lb1==0 .Or. (lb1>=25 .And. lb1<=28) .Or. (lb1>=3 .And. lb1<=5))) Then
              Goto 3699
           Else If (iabs(lb1)==30 .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Then
              Goto 3699
           End If
           If (iabs(lb2)==30 .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Then
              Goto 3699
           End If
           If ((lb1==23 .Or. lb1==21) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Then
              Goto 3699
           Else If ((lb2==23 .Or. lb2==21) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Then
              Goto 3699
           End If
           rppmax = 3.57
           If ((lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) Then
              deltr0 = rppmax
              Goto 2699
           Else If ((lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6)) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13))) Then
              deltr0 = rppmax
              Goto 2699
           End If
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 3699
           If (iabs(lb1)==42 .Or. iabs(lb2)==42) Then
              ilb1 = iabs(lb1)
              ilb2 = iabs(lb2)
              If ((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) Then
                 If ((lb1*lb2)>0) deltr0 = 9.5
              End If
           End If
           If ((iabs(lb1)>=40 .And. iabs(lb1)<=45) .Or. (iabs(lb2)>=40 .And. iabs(lb2)<=45)) Goto 3699
           If ((lb1==29 .And. ((lb2>=1 .And. lb2<=13) .Or. (lb2>=21 .And. lb2<=28) .Or. iabs(lb2)==30)) .Or. (lb2==29 .And. ((lb1>=1 .And. lb1<=13) .Or. (lb1>=21 .And. lb1<=28) .Or. iabs(lb1)==30))) Then
              deltr0 = 3.0
              Goto 3699
           End If
           If (iabs(lb1)==30 .Or. iabs(lb2)==30) Goto 400
           If (lb1==23 .And. (lb2<1 .Or. lb2>17)) Goto 400
           If (lb2==23 .And. (lb1<1 .Or. lb1>17)) Goto 400
           If (((lb1<=-1 .And. lb1>=-13) .And. (lb2==0 .Or. (lb2>=3 .And. lb2<=5) .Or. (lb2>=25 .And. lb2<=28))) .Or. ((lb2<=-1 .And. lb2>=-13) .And. (lb1==0 .Or. (lb1>=3 .And. lb1<=5) .Or. (lb1>=25 .And. lb1<=28)))) Then
           Else If (((lb1==-1 .Or. lb1==-2) .And. (lb2<-5 .And. lb2>=-13)) .Or. ((lb2==-1 .Or. lb2==-2) .And. (lb1<-5 .And. lb1>=-13))) Then
           Else If ((lb1==-1 .Or. lb1==-2) .And. (lb2==-1 .Or. lb2==-2)) Then
           Else If ((lb1<-5 .And. lb1>=-13) .And. (lb2<-5 .And. lb2>=-13)) Then
           End If
2699       Continue
           If (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=17)) Then
              If (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=17)) Then
                 deltr0 = 2.
              End If
           End If
3699       rsqare = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
           If (rsqare>deltr0**2) Goto 400
           ix2 = nint(x2/dx)
           iy2 = nint(y2/dy)
           iz2 = nint(z2/dz)
           ipx2 = nint(px2/dpx)
           ipy2 = nint(py2/dpy)
           ipz2 = nint(pz2/dpz)
           Call cms(i1, i2, pcx, pcy, pcz, srt)
           drmax = dr0max
           Call distc0(drmax, deltr0, dt, ifirst, pcx, pcy, pcz, x1, y1, z1, px1, py1, pz1, em1, x2, y2, z2, px2, py2, pz2, em2)
           If (ifirst==-1) Goto 400
           iss = nint(srt/esbin)
           If (iss>2000) iss = 2000
           If (iabs(lb1)==42 .Or. iabs(lb2)==42) Then
              ilb1 = iabs(lb1)
              ilb2 = iabs(lb2)
              If (lb1==0 .Or. (lb1>=3 .And. lb1<=5) .Or. (lb1>=25 .And. lb1<=28) .Or. lb2==0 .Or. (lb2>=3 .And. lb2<=5) .Or. (lb2>=25 .And. lb2<=28)) Then
                 Goto 505
              Else If (((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) .And. (lb1*lb2)>0) Then
                 Goto 506
              Else
                 Goto 400
              End If
           End If
           If (((lb1==23 .Or. lb1==30) .And. (lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6))) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)))) Then
              bmass = 0.938
              If (srt<=(bmass+aka)) Then
                 pkaon = 0.
              Else
                 pkaon = sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
              End If
              sigela = 0.5*(akpel(pkaon)+aknel(pkaon))
              sigsgm = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
              sig = sigela + sigsgm + akplam(pkaon)
              If (sig>1.E-7) Then
                 icase = 3
                 brel = sigela/sig
                 brsgm = sigsgm/sig
                 brsig = sig
                 nchrg = 1
                 Goto 3555
              End If
              Goto 400
           End If
           If (((lb1>=-17 .And. lb1<=-14) .And. (lb2>=3 .And. lb2<=5)) .Or. ((lb2>=-17 .And. lb2<=-14) .And. (lb1>=3 .And. lb1<=5))) Then
              nchrg = -100
              If ((lb1==-15 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==-15 .And. (lb1==5 .Or. lb1==27))) Then
                 nchrg = -2
                 bmass = 1.232
                 Goto 110
              End If
              If ((lb1==-15 .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. (lb2==-15 .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==5 .Or. lb2==27)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==5 .Or. lb1==27))) Then
                 nchrg = -1
                 bmass = 0.938
                 Goto 110
              End If
              If ((lb1==-15 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==-15 .And. (lb1==3 .Or. lb1==25)) .Or. (lb1==-17 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==-17 .And. (lb1==5 .Or. lb1==27)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28))) Then
                 nchrg = 0
                 bmass = 0.938
                 Goto 110
              End If
              If ((lb1==-17 .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. (lb2==-17 .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==3 .Or. lb2==25)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==3 .Or. lb1==25))) Then
                 nchrg = 1
                 bmass = 1.232
              End If
110           sig = 0.
              If (nchrg/=-100 .And. srt>=(aka+bmass)) Then
                 icase = 4
                 pkaon = sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
                 If (lb1==-14 .Or. lb2==-14) Then
                    If (nchrg>=0) sigma0 = akplam(pkaon)
                    If (nchrg<0) sigma0 = aknlam(pkaon)
                 Else
                    If (nchrg>=0) sigma0 = akpsgm(pkaon)
                    If (nchrg<0) sigma0 = aknsgm(pkaon)
                    sigma0 = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
                 End If
                 sig = (srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
                 If (nchrg==-2 .Or. nchrg==2) sig = 2.*sig
                 If (lb1==-14 .Or. lb2==-14) Then
                    sig = 4.0/3.0*sig
                 Else If (nchrg==-2 .Or. nchrg==2) Then
                    sig = 8.0/9.0*sig
                 Else
                    sig = 4.0/9.0*sig
                 End If
              End If
              icase = 4
              sigela = 10.
              sig = sig + sigela
              brel = sigela/sig
              brsgm = 0.
              brsig = sig
              Goto 3555
           End If
           If (((lb1==21 .Or. lb1==-30) .And. (lb2>=14 .And. lb2<=17)) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1>=14 .And. lb1<=17))) Then
              kp = 0
              Goto 3455
           End If
           If (((lb1==23 .Or. lb1==30) .And. (lb2<=-14 .And. lb2>=-17)) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1<=-14 .And. lb1>=-17))) Then
              kp = 1
              Goto 3455
           End If
           If (((lb1==21 .Or. lb1==-30) .And. (lb2==40 .Or. lb2==41)) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1==40 .Or. lb1==41))) Then
              kp = 0
              Goto 3455
           End If
           If (((lb1==23 .Or. lb1==30) .And. (lb2==-40 .Or. lb2==-41)) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==-40 .Or. lb1==-41))) Then
              kp = 1
              Goto 3455
           End If
           kp = 3
           If ((((lb1>=3 .And. lb1<=5) .Or. lb1==0) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) .Or. (((lb2>=3 .And. lb2<=5) .Or. lb2==0) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41))) Goto 3455
           If (((lb1>=3 .And. lb1<=5) .And. iabs(lb2)==45) .Or. ((lb2>=3 .And. lb2<=5) .And. iabs(lb1)==45)) Goto 3455
           If (lb1==23 .And. (lb2>=14 .And. lb2<=17)) Goto 5699
           If (lb2==23 .And. (lb1>=14 .And. lb1<=17)) Goto 5699
           If (lb1==21 .And. (lb2>=-17 .And. lb2<=-14)) Goto 5699
           If (lb2==21 .And. (lb1>=-17 .And. lb1<=-14)) Goto 5699
           If ((((lb1==1 .Or. lb1==2) .Or. (lb1>=6 .And. lb1<=13)) .And. (lb2>=-17 .And. lb2<=-14)) .Or. (((lb2==1 .Or. lb2==2) .Or. (lb2>=6 .And. lb2<=13)) .And. (lb1>=-17 .And. lb1<=-14))) Goto 5999
           If ((((lb1==-1 .Or. lb1==-2) .Or. (lb1<=-6 .And. lb1>=-13)) .And. (lb2>=14 .And. lb2<=17)) .Or. (((lb2==-1 .Or. lb2==-2) .Or. (lb2<=-6 .And. lb2>=-13)) .And. (lb1>=14 .And. lb1<=17))) Goto 5999
           If (lb1==21 .And. lb2==23) Goto 8699
           If (lb2==21 .And. lb1==23) Goto 8699
           If (lb1==30 .And. lb2==21) Goto 8699
           If (lb2==30 .And. lb1==21) Goto 8699
           If (lb1==-30 .And. lb2==23) Goto 8699
           If (lb2==-30 .And. lb1==23) Goto 8699
           If (lb1==-30 .And. lb2==30) Goto 8699
           If (lb2==-30 .And. lb1==30) Goto 8699
           If (((lb1==23 .Or. lb1==21 .Or. iabs(lb1)==30) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==23 .Or. lb2==21 .Or. iabs(lb2)==30) .And. (lb1>=25 .And. lb1<=28))) Goto 8799
           If ((iabs(lb1)==30 .And. (lb2>=3 .And. lb2<=5)) .Or. (iabs(lb2)==30 .And. (lb1>=3 .And. lb1<=5))) Goto 8799
           If ((lb1==29 .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=9))) .Or. (lb2==29 .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=9)))) Goto 7222
           If ((lb1==29 .And. ((lb2>=3 .And. lb2<=5) .Or. (lb2>=21 .And. lb2<=28) .Or. iabs(lb2)==30)) .Or. (lb2==29 .And. ((lb1>=3 .And. lb1<=5) .Or. (lb1>=21 .And. lb1<=28) .Or. iabs(lb1)==30))) Then
              Goto 7444
           End If
           If (((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. iabs(lb1)>=40) .And. ((lb2>=25 .And. lb2<=29) .Or. lb2==0)) Goto 888
           If (((iabs(lb2)>=14 .And. iabs(lb2)<=17) .Or. iabs(lb2)>=40) .And. ((lb1>=25 .And. lb1<=29) .Or. lb1==0)) Goto 888
           If (((lb1==23 .Or. lb1==30) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13)))) Goto 888
           If (((lb1==21 .Or. lb1==-30) .And. (lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6))) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)))) Goto 888
           If (((lb1>=14 .And. lb1<=17) .And. (lb2>=6 .And. lb2<=13)) .Or. ((lb2>=14 .And. lb2<=17) .And. (lb1>=6 .And. lb1<=13))) Goto 7799
           If (((lb1<=-14 .And. lb1>=-17) .And. (lb2<=-6 .And. lb2>=-13)) .Or. ((lb2<=-14 .And. lb2>=-17) .And. (lb1<=-6 .And. lb1>=-13))) Goto 7799
           If (iabs(lb1)>=40 .Or. iabs(lb2)>=40 .Or. (lb1<=-14 .And. lb1>=-17) .Or. (lb2<=-14 .And. lb2>=-17)) Goto 400
           If ((lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) Then
              Goto 2799
           Else If ((lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6)) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13))) Then
              Goto 2799
           End If
           inewka = irun
           Call newka(icase, inewka, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, iblock)
           If (ictrl==1) Goto 400
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 400
           If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If (lb1==0 .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If (lb2==0 .And. (lb1>=3 .And. lb1<=5)) Goto 777
           If (lb1==0 .And. lb2==0) Goto 777
           If ((lb1>=25 .And. lb1<=28) .And. (lb2>=25 .And. lb2<=28)) Goto 777
           If ((lb1>=25 .And. lb1<=28) .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If ((lb2>=25 .And. lb2<=28) .And. (lb1>=3 .And. lb1<=5)) Goto 777
           If ((lb1>=25 .And. lb1<=28) .And. lb2==0) Goto 777
           If ((lb2>=25 .And. lb2<=28) .And. lb1==0) Goto 777
           If ((lb1==23 .Or. lb1==21) .And. (lb2>=3 .And. lb2<=5)) Goto 889
           If ((lb2==23 .Or. lb2==21) .And. (lb1>=3 .And. lb1<=5)) Goto 889
           If (iabs(lb1)==30 .Or. iabs(lb2)==30) Goto 400
           If (lb1==21 .Or. lb2==21) Goto 400
           If (lb1==23 .Or. lb2==23) Goto 400
           If ((lb1>=3 .And. lb1<=5) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 3
           If ((lb2>=3 .And. lb2<=5) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 3
           If ((lb1>=25 .And. lb1<=28) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 33
           If ((lb2>=25 .And. lb2<=28) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 33
           If (lb1==0 .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 547
           If (lb2==0 .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 547
           If ((lb1==1 .Or. lb1==2) .And. (lb2>5 .And. lb2<=13)) Goto 44
           If ((lb2==1 .Or. lb2==2) .And. (lb1>5 .And. lb1<=13)) Goto 44
           If ((lb1==-1 .Or. lb1==-2) .And. (lb2<-5 .And. lb2>=-13)) Goto 44
           If ((lb2==-1 .Or. lb2==-2) .And. (lb1<-5 .And. lb1>=-13)) Goto 44
           If ((lb1==1 .Or. lb1==2) .And. (lb2==1 .Or. lb2==2)) Goto 4
           If ((lb1==-1 .Or. lb1==-2) .And. (lb2==-1 .Or. lb2==-2)) Goto 4
           If ((lb1>5 .And. lb1<=13) .And. (lb2>5 .And. lb2<=13)) Goto 444
           If ((lb1<-5 .And. lb1>=-13) .And. (lb2<-5 .And. lb2>=-13)) Goto 444
           If ((lb1<3) .And. (lb2>=14 .And. lb2<=17)) Goto 400
           If ((lb2<3) .And. (lb1>=14 .And. lb1<=17)) Goto 400
           If ((lb1>=14 .And. lb1<=17) .And. (lb2>=14 .And. lb2<=17)) Goto 400
           Goto 400
547        If (lb1*lb2==0) Then
              ece = (em1+em2+0.02)**2
              xkaon0 = 0.
              If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
              If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
              xkaon0 = 2.0*xkaon0
              xkaon = xkaon0
              xeta = xn1535(i1, i2, 0)
              If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .Or. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) xeta = 0.
              If ((xeta+xkaon)<=1.E-06) Goto 400
              dse = sqrt((xeta+xkaon)/pi)
              deltre = dse + 0.1
              px1cm = pcx
              py1cm = pcy
              pz1cm = pcz
              Call distce(i1, i2, deltre, dse, dt, ece, srt, ic, pcx, pcy, pcz)
              If (ic==-1) Goto 400
              ekaon(4, iss) = ekaon(4, iss) + 1
              If (xkaon0/(xkaon+xeta)>ranart(nseed)) Then
                 Call cren(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
                 If (iblock==7) Then
                    lpn = lpn + 1
                 Else If (iblock==-7) Then
                 End If
                 em1 = e(i1)
                 em2 = e(i2)
                 Goto 440
              End If
              resona = 1.
              Goto 98
           End If
3          Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           xkaon0 = 0.
           If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
           If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
           xkaon0 = 2.0*xkaon0
           xphi = 0.
           If ((((lb1>=1 .And. lb1<=2) .Or. (lb1>=6 .And. lb1<=9)) .Or. ((lb2>=1 .And. lb2<=2) .Or. (lb2>=6 .And. lb2<=9))) .And. srt>1.958) Call pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
           If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .Or. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Goto 31
           ec = (em1+em2+0.02)**2
           xkaon = 0.
           If (srt>1.23) xkaon = (pionpp(srt)+pipp1(srt))/2.
           If ((lb1*lb2==5 .Or. ((lb1*lb2==6) .And. (lb1==3 .Or. lb2==3))) .Or. (lb1*lb2==-3 .Or. ((lb1*lb2==-10) .And. (lb1==5 .Or. lb2==5)))) Then
              xmax = 190.
              xmaxn = 0
              xmaxn1 = 0
              xdirct = dirct1(srt)
              Goto 678
           End If
           If ((lb1*lb2==3 .Or. ((lb1*lb2==10) .And. (lb1==5 .Or. lb2==5))) .Or. (lb1*lb2==-5 .Or. ((lb1*lb2==-6) .And. (lb1==3 .Or. lb2==3)))) Then
              xmax = 27.
              xmaxn = 2./3.*25.*0.6
              xmaxn1 = 2./3.*40.*0.5
              xdirct = dirct2(srt)
              Goto 678
           End If
           If ((lb1==4 .Or. lb2==4) .And. (iabs(lb1*lb2)==4 .Or. iabs(lb1*lb2)==8)) Then
              xmax = 50.
              xmaxn = 1./3.*25*0.6
              xmaxn1 = 1/3.*40.*0.5
              xdirct = dirct3(srt)
              Goto 678
           End If
678        xnpin1 = 0
           xnpin = 0
           xnpid = xnpi(i1, i2, 1, xmax)
           If (xmaxn1/=0) xnpin1 = xnpi(i1, i2, 2, xmaxn1)
           If (xmaxn/=0) xnpin = xnpi(i1, i2, 0, xmaxn)
           xres = xnpid + xnpin + xnpin1
           xnelas = xres + xdirct
           icheck = 1
           Goto 34
31         ec = (em1+em2+0.02)**2
           xreab = reab(i1, i2, srt, 1)
           If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) xreab = 0.
           xkaon = xkaon0 + xreab
           If ((iabs(lb1)>9 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>9 .And. iabs(lb2)<=13)) Then
              xnelas = 1.0
           Else
              xnelas = dpion(em1, em2, lb1, lb2, srt)
           End If
           icheck = 2
34         If ((xnelas+xkaon+xphi)<=0.000001) Goto 400
           ds = sqrt((xnelas+xkaon+xphi)/pi)
           deltar = ds + 0.1
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, pcx, pcy, pcz)
           If (ic==-1) Goto 400
           ekaon(4, iss) = ekaon(4, iss) + 1
           If (icheck==2) Then
              If (xnelas/(xnelas+xkaon+xphi)>=ranart(nseed)) Then
                 Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
                 Goto 440
              Else
                 Goto 96
              End If
           End If
           If ((xkaon+xphi)/(xkaon+xphi+xnelas)>ranart(nseed)) Goto 95
           If (xdirct/xnelas>=ranart(nseed)) Then
              Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
              Goto 440
           End If
           If ((lb1*lb2==5 .Or. ((lb1*lb2==6) .And. (lb1==3 .Or. lb2==3))) .Or. (lb1*lb2==-3 .Or. ((lb1*lb2==-10) .And. (lb1==5 .Or. lb2==5)))) Then
              Goto 99
           Else
              xx = (xnpin+xnpin1)/xres
              If (ranart(nseed)<xx) Then
                 xx0 = xnpin/(xnpin+xnpin1)
                 If (ranart(nseed)<xx0) Then
                    resona = 0.
                    Goto 97
                 Else
                    resona = 1.
                    Goto 98
                 End If
              Else
                 Goto 99
              End If
           End If
97         Continue
           If (resona==0.) Then
              i = i1
              If (em1<0.6) i = i2
              If ((lb1*lb2==10 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-6 .And. (lb1==3 .Or. lb2==3))) Then
                 lb(i) = 11
                 Goto 303
              End If
              If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 11
                 Goto 303
              End If
              If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 10
                 Goto 303
              End If
              If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
                 lb(i) = 10
              End If
303           Call dreson(i1, i2)
              If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
              lres = lres + 1
              Goto 101
           End If
98         If (resona==1.) Then
              i = i1
              If (em1<0.6) i = i2
              If ((lb1*lb2==10 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-6 .And. (lb1==3 .Or. lb2==3))) Then
                 lb(i) = 13
                 Goto 304
              End If
              If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 13
                 Goto 304
              End If
              If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 12
                 Goto 304
              End If
              If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
                 lb(i) = 12
                 Goto 304
              End If
              If (lb(i1)*lb(i2)==0) Then
                 If (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1) Then
                    lb(i) = 13
                    Goto 304
                 Else
                    lb(i) = 12
                 End If
              End If
304           Call dreson(i1, i2)
              If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
              lres = lres + 1
              Goto 101
           End If
99         lres = lres + 1
           i = i1
           If (em1<=0.6) i = i2
           If ((lb(i1)*lb(i2)==5) .Or. (lb(i1)*lb(i2)==-3)) Then
              lb(i) = 9
              Goto 305
           End If
           If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
              lb(i) = 8
              Goto 305
           End If
           If ((lb(i1)*lb(i2)==10 .And. (lb(i1)==5 .Or. lb(i2)==5)) .Or. (lb(i1)*lb(i2)==-6 .And. (lb(i1)==3 .Or. lb(i2)==3))) Then
              lb(i) = 8
              Goto 305
           End If
           If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
              lb(i) = 7
              Goto 305
           End If
           If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
              lb(i) = 7
              Goto 305
           End If
           If ((lb(i1)*lb(i2)==6 .And. (lb(i1)==3 .Or. lb(i2)==3)) .Or. (lb(i1)*lb(i2)==-10 .And. (lb(i1)==5 .Or. lb(i2)==5))) Then
              lb(i) = 6
           End If
305        Call dreson(i1, i2)
           If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
           Goto 101
889        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           spika = 60./(1.+4.*(srt-0.895)**2/(0.05)**2)
           Call crkpla(px1cm, py1cm, pz1cm, ec, srt, spika, emm1, emm2, lbp1, lbp2, i1, i2, icase, srhoks)
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If
           If (icase==1) Then
              Call ksreso(i1, i2)
              iblock = 171
              lres = lres + 1
              Goto 101
           Else If (icase==2) Then
              iblock = 71
           Else If (iabs(icase)==5) Then
              iblock = 88
           Else
              iblock = 222
           End If
           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           ntag = 0
           Goto 440
33         Continue
           em1 = e(i1)
           em2 = e(i2)
           xelstc = 0
           If ((lb1>=25 .And. lb1<=28) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2)) xelstc = erhon(em1, em2, lb1, lb2, srt)
           If ((lb2>=25 .And. lb2<=28) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2)) xelstc = erhon(em1, em2, lb1, lb2, srt)
           ec = (em1+em2+0.02)**2
           xkaon0 = 0
           If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
           If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
           If (xkaon0<0) xkaon0 = 0
           xkaon0 = 2.0*xkaon0
           xkaon = xkaon0
           ichann = 0
           xphi = 0.
           If (((((lb1>=1 .And. lb1<=2) .Or. (lb1>=6 .And. lb1<=9)) .And. (lb2>=25 .And. lb2<=27)) .Or. (((lb2>=1 .And. lb2<=2) .Or. (lb2>=6 .And. lb2<=9)) .And. (lb1>=25 .And. lb1<=27))) .And. srt>1.958) Call pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
           If ((iabs(lb1)>=6 .And. lb2>=25) .Or. (lb1>=25 .And. iabs(lb2)>=6)) Then
              ichann = 1
              ictrl = 2
              If (lb1==28 .Or. lb2==28) ictrl = 3
              xreab = reab(i1, i2, srt, ictrl)
              If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) xreab = 0.
              If (xreab<0) xreab = 1.E-06
              xkaon = xkaon0 + xreab
              xelstc = 1.0
           End If
           ds = sqrt((xkaon+xphi+xelstc)/pi)
           deltar = ds + 0.1
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, pcx, pcy, pcz)
           If (ic==-1) Goto 400
           ekaon(4, iss) = ekaon(4, iss) + 1
           If (xelstc/(xelstc+xkaon+xphi)>ranart(nseed)) Then
              Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
              Goto 440
           End If
           Call crrd(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           If (iblock==81) lrhor = lrhor + 1
           If (iblock==82) lomgar = lomgar + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
95         Continue
           Call crpn(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           If (iblock==77) lpd = lpd + 1
           If (iblock==78) lrho = lrho + 1
           If (iblock==79) lomega = lomega + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
96         Continue
           Call crpd(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           If (iblock==80) lpdr = lpdr + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
101        Continue
           If (e(i2)==0.) Goto 600
           If (e(i1)==0.) Goto 800
44         Continue
           cutoff = em1 + em2 + 0.02
           If (srt<=cutoff) Goto 400
           If (srt>2.245) Then
              signn = pp2(srt)
           Else
              signn = 35.0/(1.+(srt-cutoff)*100.0) + 20.0
           End If
           Call xnd(pcx, pcy, pcz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
           sig = signn + xinel
           ec = (em1+em2+0.02)**2
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           ipdflag = 0
           If (idpert==1) Then
              ipert1 = 1
              sigr0 = sig
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 363
              signn0 = 0.
              Call crnd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, signn0, sigr0, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
              ipdflag = 1
363           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           ds = sqrt(sig/(10.*pi))
           deltar = ds + 0.1
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           ekaon(3, iss) = ekaon(3, iss) + 1
           Goto 361
361        Continue
           Call crnd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, signn, sig, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
           If (iblock==0 .And. ipdflag==1) iblock = 501
           If (iblock==11) Then
              lndk = lndk + 1
              Goto 400
           Else If (iblock==-11 .Or. iblock==501) Then
              Goto 400
           End If
           If (iblock==222) Then
              Goto 400
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
4          Continue
           cutoff = em1 + em2 + 0.14
           If (srt>2.245) Then
              sig = ppt(srt)
              signn = sig - pp1(srt)
           Else
              sig = xpp(srt)
              If (zet(lb(i1))*zet(lb(i2))<=0) sig = xnp(srt)
              If (zet(lb(i1))*zet(lb(i2))>0) sig = xpp(srt)
              If (zet(lb(i1))==0 .And. zet(lb(i2))==0) sig = xpp(srt)
              If ((lb(i1)==-1 .And. lb(i2)==-2) .Or. (lb(i2)==-1 .And. lb(i1)==-2)) sig = xnp(srt)
              If (srt<1.897) Then
                 signn = sig
              Else
                 signn = 35.0/(1.+(srt-1.897)*100.0) + 20.0
              End If
           End If
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           ipdflag = 0
           If (idpert==1) Then
              ipert1 = 1
              ec = 2.012**2
              sigr0 = sig
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 365
              signn0 = 0.
              Call crnn(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn0, sigr0, nt, ipert1)
              ipdflag = 1
365           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           If (signn<=0) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           ec = 3.59709
           ds = sqrt(sig/pi/10.)
           dsr = ds + 0.1
           If ((e(i1)>=1.) .And. (e(i2)>=1.)) ec = 4.75
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           Goto 362
362        ekaon(1, iss) = ekaon(1, iss) + 1
           Call crnn(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
           If (iblock==0 .And. ipdflag==1) iblock = 501
           If (iblock==4 .Or. iblock==9 .Or. iblock>=44 .Or. iblock==-9 .Or. iblock==222 .Or. iblock==501) Then
              lcoll = lcoll + 1
              If (iblock==4) Then
                 ldirt = ldirt + 1
              Else If (iblock==44) Then
                 lddrho = lddrho + 1
              Else If (iblock==45) Then
                 lnnrho = lnnrho + 1
              Else If (iblock==46) Then
                 lnnom = lnnom + 1
              Else If (iblock==222) Then
              Else If (iblock==9) Then
                 lnnk = lnnk + 1
              Else If (iblock==-9) Then
              End If
              Goto 400
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
505        Continue
           ianti = 0
           If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
           Call sdmbb(srt, sdm, ianti)
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = 2.012**2
           ds = sqrt(sdm/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crdmbb(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, sdm, nt, ianti)
           lcoll = lcoll + 1
           Goto 400
506        Continue
           ianti = 0
           If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
           Call sdbelastic(srt, sdb)
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = 2.012**2
           ds = sqrt(sdb/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crdbel(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, sdb, nt, ianti)
           lcoll = lcoll + 1
           Goto 400
444        Continue
           cutoff = em1 + em2 + 0.02
           If (srt<=cutoff) Goto 400
           If (srt>2.245) Then
              signn = pp2(srt)
           Else
              signn = 35.0/(1.+(srt-cutoff)*100.0) + 20.0
           End If
           If (signn<=0) Goto 400
           Call xddin(pcx, pcy, pcz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
           sig = signn + xinel
           ec = (em1+em2+0.02)**2
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           ipdflag = 0
           If (idpert==1) Then
              ipert1 = 1
              sigr0 = sig
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 367
              signn0 = 0.
              Call crdd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn0, sigr0, nt, ipert1)
              ipdflag = 1
367           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           ds = sqrt(sig/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           Goto 364
364        ekaon(2, iss) = ekaon(2, iss) + 1
           Call crdd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
           If (iblock==0 .And. ipdflag==1) iblock = 501
           If (iabs(iblock)==10) Then
              lcoll = lcoll + 1
              If (iblock==10) Then
                 lddk = lddk + 1
              Else If (iblock==-10) Then
              End If
              Goto 400
           End If
           If (iblock==222 .Or. iblock==501) Then
              Goto 400
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
777        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec0 = em1 + em2 + 0.02
           If (srt<=ec0) Goto 400
           ec = (em1+em2+0.02)**2
           ppel = 20.
           ipp = 1
           If (lb1<3 .Or. lb1>5 .Or. lb2<3 .Or. lb2>5) Goto 778
           Call ppxs(lb1, lb2, srt, ppsig, spprho, ipp)
           ppel = ppsig
778        ppink = pipik(srt)
           ppink = 2.0*ppink
           If (lb1>=25 .And. lb2>=25) ppink = rrkk
           If (((lb1==0 .Or. (lb1>=3 .And. lb1<=5)) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==0 .Or. (lb2>=3 .And. lb2<=5)) .And. (lb1>=25 .And. lb1<=28))) Then
              ppink = 0.
              If (srt>=(aka+aks)) ppink = prkk
           End If
           Call spprr(lb1, lb2, srt)
           Call sppee(lb1, lb2, srt)
           Call spppe(lb1, lb2, srt)
           Call srpre(lb1, lb2, srt)
           Call sopoe(lb1, lb2, srt)
           Call srree(lb1, lb2, srt)
           ppinnb = 0.
           If (srt>thresh(1)) Then
              Call getnst(srt)
              If (lb1>=3 .And. lb1<=5 .And. lb2>=3 .And. lb2<=5) Then
                 ppinnb = ppbbar(srt)
              Else If ((lb1>=3 .And. lb1<=5 .And. lb2>=25 .And. lb2<=27) .Or. (lb2>=3 .And. lb2<=5 .And. lb1>=25 .And. lb1<=27)) Then
                 ppinnb = prbbar(srt)
              Else If (lb1>=25 .And. lb1<=27 .And. lb2>=25 .And. lb2<=27) Then
                 ppinnb = rrbbar(srt)
              Else If ((lb1>=3 .And. lb1<=5 .And. lb2==28) .Or. (lb2>=3 .And. lb2<=5 .And. lb1==28)) Then
                 ppinnb = pobbar(srt)
              Else If ((lb1>=25 .And. lb1<=27 .And. lb2==28) .Or. (lb2>=25 .And. lb2<=27 .And. lb1==28)) Then
                 ppinnb = robbar(srt)
              Else If (lb1==28 .And. lb2==28) Then
                 ppinnb = oobbar(srt)
              Else
                 If (lb1/=0 .And. lb2/=0) Write (6, *) 'missed MM lb1,lb2=', lb1, lb2
              End If
           End If
           ppin = ppink + ppinnb + pprr + ppee + pppe + rpre + xopoe + rree
           If ((ppel+ppin)<=0.01) Goto 400
           dspp = sqrt((ppel+ppin)/31.4)
           dsppr = dspp + 0.1
           Call distce(i1, i2, dsppr, dspp, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           If (ppel==0) Goto 400
           ekaon(5, iss) = ekaon(5, iss) + 1
           Call crpp(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ppel, ppin, spprho, ipp)
           If (iblock==666) Goto 555
           If (iblock==6) lpp = lpp + 1
           If (iblock==66) Then
              lppk = lppk + 1
           Else If (iblock==366) Then
              lppk = lppk + 1
           Else If (iblock==367) Then
              lppk = lppk + 1
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
2799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           dsppb = sqrt(xppbar(srt)/pi/10.)
           dsppbr = dsppb + 0.1
           Call distce(i1, i2, dsppbr, dsppb, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crppba(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
3555       px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           dskk = sqrt(sig/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crlaba(px1cm, py1cm, pz1cm, srt, brel, brsgm, i1, i2, nt, iblock, nchrg, icase)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
3455       px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           Call pertur(px1cm, py1cm, pz1cm, srt, irun, i1, i2, nt, kp, icontp)
           If (icontp==0) Then
              em1 = e(i1)
              em2 = e(i2)
              iblock = 727
              Goto 440
           End If
           If (e(i1)==0.) Goto 800
           If (e(i2)==0.) Goto 600
           Goto 400
7222       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call xphib(lb1, lb2, em1, em2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, sigp)
           dskk = sqrt(sigp/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crphib(px1cm, py1cm, pz1cm, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, sigp, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
7444       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call phimes(i1, i2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, sigphi)
           dskk = sqrt(sigphi/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           pzrt = p(3, i1) + p(3, i2)
           er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
           er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
           ert = er1 + er2
           yy = 0.5*log((ert+pzrt)/(ert-pzrt))
           Call crphim(px1cm, py1cm, pz1cm, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, sigphi, ikkg, ikkl, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
7799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call lambar(i1, i2, srt, siglab)
           dshn = sqrt(siglab/pi/10.)
           dshnr = dshn + 0.1
           Call distce(i1, i2, dshnr, dshn, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crhb(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
5699       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call xkhype(i1, i2, srt, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk)
           dskk = sqrt(sigk/pi)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           If (lb(i1)==23 .Or. lb(i2)==23) Then
              ikmp = 1
           Else
              ikmp = -1
           End If
           Call crkhyp(px1cm, py1cm, pz1cm, srt, i1, i2, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk, ikmp, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
5999       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           sigkp = 15.
           dskk = sqrt(sigkp/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crlan(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
8699       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call crkphi(px1cm, py1cm, pz1cm, ec, srt, iblock, emm1, emm2, lbp1, lbp2, i1, i2, ikk, icase, rrkk, prkk)
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If
           If (lbp1==29 .Or. lbp2==29) Then
              pzrt = p(3, i1) + p(3, i2)
              er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
              er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
              ert = er1 + er2
              yy = 0.5*log((ert+pzrt)/(ert-pzrt))
              iblock = 222
              ntag = 0
           End If
           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
8799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call crksph(px1cm, py1cm, pz1cm, ec, srt, emm1, emm2, lbp1, lbp2, i1, i2, ikkg, ikkl, iblock, icase, srhoks)
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If
           If (lbp1==29 .Or. lbp2==20) Then
              pzrt = p(3, i1) + p(3, i2)
              er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
              er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
              ert = er1 + er2
              yy = 0.5*log((ert+pzrt)/(ert-pzrt))
           End If
           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
888        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           sig = 10.
           If (iabs(lb1)==14 .Or. iabs(lb2)==14 .Or. iabs(lb1)==30 .Or. iabs(lb2)==30) sig = 20.
           If (lb1==29 .Or. lb2==29) sig = 5.0
           dskn = sqrt(sig/pi/10.)
           dsknr = dskn + 0.1
           Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crkn(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
440        Continue
           If (iblock==0) Goto 400
           lcoll = lcoll + 1
           ntag = 0
           e1cm = sqrt(em1**2+px1cm**2+py1cm**2+pz1cm**2)
           p1beta = px1cm*betax + py1cm*betay + pz1cm*betaz
           transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
           pt1i1 = betax*transf + px1cm
           pt2i1 = betay*transf + py1cm
           pt3i1 = betaz*transf + pz1cm
           Goto 90002
90002      Continue
           e2cm = sqrt(em2**2+px1cm**2+py1cm**2+pz1cm**2)
           transf = gamma*(-gamma*p1beta/(gamma+1.)+e2cm)
           pt1i2 = betax*transf - px1cm
           pt2i2 = betay*transf - py1cm
           pt3i2 = betaz*transf - pz1cm
           Goto 90003
90003      If (iblock==1) lcnne = lcnne + 1
           If (iblock==5) ldd = ldd + 1
           If (iblock==2) lcnnd = lcnnd + 1
           If (iblock==8) lkn = lkn + 1
           If (iblock==43) ldou = ldou + 1
           If (iblock==3) lcndn = lcndn + 1
           p(1, i1) = pt1i1
           p(2, i1) = pt2i1
           p(3, i1) = pt3i1
           p(1, i2) = pt1i2
           p(2, i2) = pt2i2
           p(3, i2) = pt3i2
           px1 = p(1, i1)
           py1 = p(2, i1)
           pz1 = p(3, i1)
           em1 = e(i1)
           em2 = e(i2)
           lb1 = lb(i1)
           lb2 = lb(i2)
           id(i1) = 2
           id(i2) = 2
           e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
           id1 = id(i1)
           Goto 90004
90004      Continue
           am1 = em1
           am2 = em2
400        Continue
555        Continue
600     End Do
798     If (nt==ntmax .And. ipi0dcy==1 .And. i1==(massr(irun)+msum)) Then
           Do ipion = 1, nnn
              If (lpion(ipion,irun)==4) Then
                 wid = 7.85E-9
                 Call resdec(i1, nt, nnn, wid, idecay, ipion)
              End If
           End Do
        End If
800  End Do
     n0 = mass + msum
     Do n = n0 + 1, massr(irun) + msum
        If (e(n)>0. .Or. lb(n)>5000) Then
           nnn = nnn + 1
           rpion(1, nnn, irun) = r(1, n)
           rpion(2, nnn, irun) = r(2, n)
           rpion(3, nnn, irun) = r(3, n)
           If (nt==ntmax) Then
              ftpisv(nnn, irun) = ftsv(n)
              tfdpi(nnn, irun) = tfdcy(n)
           End If
           ppion(1, nnn, irun) = p(1, n)
           ppion(2, nnn, irun) = p(2, n)
           ppion(3, nnn, irun) = p(3, n)
           epion(nnn, irun) = e(n)
           lpion(nnn, irun) = lb(n)
           propi(nnn, irun) = proper(n)
           dppion(nnn, irun) = dpertp(n)
        End If
     End Do
     massrn(irun) = nnn + mass
  End Do
  ia = 0
  ib = 0
  Do irun = 1, num
     ia = ia + massr(irun-1)
     ib = ib + massrn(irun-1)
     Do ic = 1, massrn(irun)
        ie = ia + ic
        ig = ib + ic
        If (ic<=mass) Then
           rt(1, ig) = r(1, ie)
           rt(2, ig) = r(2, ie)
           rt(3, ig) = r(3, ie)
           If (nt==ntmax) Then
              fttemp(ig) = ftsv(ie)
              tft(ig) = tfdcy(ie)
           End If
           pt(1, ig) = p(1, ie)
           pt(2, ig) = p(2, ie)
           pt(3, ig) = p(3, ie)
           et(ig) = e(ie)
           lt(ig) = lb(ie)
           prot(ig) = proper(ie)
           dptemp(ig) = dpertp(ie)
        Else
           i0 = ic - mass
           rt(1, ig) = rpion(1, i0, irun)
           rt(2, ig) = rpion(2, i0, irun)
           rt(3, ig) = rpion(3, i0, irun)
           If (nt==ntmax) Then
              fttemp(ig) = ftpisv(i0, irun)
              tft(ig) = tfdpi(i0, irun)
           End If
           pt(1, ig) = ppion(1, i0, irun)
           pt(2, ig) = ppion(2, i0, irun)
           pt(3, ig) = ppion(3, i0, irun)
           et(ig) = epion(i0, irun)
           lt(ig) = lpion(i0, irun)
           prot(ig) = propi(i0, irun)
           dptemp(ig) = dppion(i0, irun)
        End If
     End Do
  End Do
  il = 0
  Do irun = 1, num
     massr(irun) = massrn(irun)
     il = il + massr(irun-1)
     Do im = 1, massr(irun)
        in = il + im
        r(1, in) = rt(1, in)
        r(2, in) = rt(2, in)
        r(3, in) = rt(3, in)
        If (nt==ntmax) Then
           ftsv(in) = fttemp(in)
           tfdcy(in) = tft(in)
        End If
        p(1, in) = pt(1, in)
        p(2, in) = pt(2, in)
        p(3, in) = pt(3, in)
        e(in) = et(in)
        lb(in) = lt(in)
        proper(in) = prot(in)
        dpertp(in) = dptemp(in)
        If (lb(in)<1 .Or. lb(in)>2) id(in) = 0
     End Do
  End Do
  Return
End Subroutine relcol
