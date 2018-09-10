Subroutine artmn
  Parameter (maxstr=150001, maxr=1, amu=0.9383, aka=0.498, etam=0.5475)
  Parameter (maxx=20, maxz=24)
  Parameter (isum=1001, igam=1100)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Integer outpar, zta, zpr
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Common /hh/proper(maxstr)
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /input/nstar, ndirct, dir
  Common /pp/prho(-20:20, -24:24)
  Common /qq/phrho(-maxz:maxz, -24:24)
  Common /rr/massr(0:maxr)
  Common /ss/inout(20)
  Common /zz/zta, zpr
  Common /run/num
  Common /kkk/tkaon(7), ekaon(7, 0:2000)
  Common /kaon/ak(3, 50, 36), speck(50, 36, 7), mf
  Common /table/xarray(0:1000), earray(0:1000)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /ddpi/pirho(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Dimension temp(3, maxstr), skaon(7), sekaon(7, 0:2000)
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /input3/plab, elab, zeropt, b0, bi, bm, dencut, cycbox
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arercp/pro1(maxstr, maxr)
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Dimension npi(maxr)
  Dimension rt(3, maxstr, maxr), pt(3, maxstr, maxr), et(maxstr, maxr), lt(maxstr, maxr), prot(maxstr, maxr)
  External iarflv, invflv
  Common /lastt/itimeh, bimp
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  Common /resdcy/nsav, iksdcy
  Common /rndf77/nseed
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./
  nlast = 0
  Do i = 1, maxstr
     ftsv(i) = 0.
     Do irun = 1, maxr
        ftsvt(i, irun) = 0.
     End Do
     lblast(i) = 999
     Do j = 1, 4
        xlast(j, i) = 0.
        plast(j, i) = 0.
     End Do
  End Do
  Call tablem
  ikaon = 1
  nstar = 1
  ndirct = 0
  dir = 0.02
  asy = 0.032
  esbin = 0.04
  mf = 36
  radta = 1.124*float(massta)**(1./3.)
  radpr = 1.124*float(masspr)**(1./3.)
  zdist = radta + radpr
  bmax = radta + radpr
  mass = massta + masspr
  ntotal = num*mass
  If (ntotal>maxstr) Then
     Write (12, '(//10X,''**** FATAL ERROR: TOO MANY TEST PART. ****'//' '')')
     Stop
  End If
  eta = float(massta)*amu
  pzta = 0.0
  betata = 0.0
  gammta = 1.0
  epr = float(masspr)*(amu+0.001*elab)
  pzpr = sqrt(epr**2-(amu*float(masspr))**2)
  betapr = pzpr/epr
  gammpr = 1.0/sqrt(1.0-betapr**2)
  betac = (pzpr+pzta)/(epr+eta)
  gammc = 1.0/sqrt(1.-betac**2)
  If (insys/=0) Then
     s = (epr+eta)**2 - pzpr**2
     xx1 = 4.*alog(float(massta))
     xx2 = 4.*alog(float(masspr))
     xx1 = exp(xx1)
     xx2 = exp(xx2)
     psqare = (s**2+(xx1+xx2)*amu**4-2.0*s*amu**2*float(massta**2+masspr**2)-2.0*float(massta**2*masspr**2)*amu**4)/(4.0*s)
     eta = sqrt(psqare+(float(massta)*amu)**2)
     pzta = -sqrt(psqare)
     betata = pzta/eta
     gammta = 1.0/sqrt(1.0-betata**2)
     epr = sqrt(psqare+(float(masspr)*amu)**2)
     pzpr = sqrt(psqare)
     betapr = pzpr/epr
     gammpr = 1.0/sqrt(1.0-betapr**2)
  Else
  End If
  pzta = pzta/float(massta)
  pzpr = pzpr/float(masspr)
  ecms0 = eta + epr
  Do imany = 1, manyb
     If (manyb>1) Then
111     bx = 1.0 - 2.0*ranart(nseed)
        by = 1.0 - 2.0*ranart(nseed)
        b2 = bx*bx + by*by
        If (b2>1.0) Goto 111
        b = sqrt(b2)*(bm-bi) + bi
     Else
        b = b0
     End If
     Call coulin(masspr, massta, num)
     Call init(1, massta, num, radta, b/2., zeropt+zdist/2., pzta, gammta, iseed, mass, imomen)
     Call init(1+massta, mass, num, radpr, -b/2., zeropt-zdist/2., pzpr, gammpr, iseed, mass, imomen)
     outpar = 0
     massr(0) = 0
     Do ir = 1, num
        massr(ir) = mass
     End Do
     Call dens(ipot, mass, num, outpar)
     If (icoll/=-1) Then
        Do i = 1, ntotal
           ix = nint(r(1,i))
           iy = nint(r(2,i))
           iz = nint(r(3,i))
           If (ix>=maxx .Or. iy>=maxx .Or. iz>=maxz .Or. ix<=-maxx .Or. iy<=-maxx .Or. iz<=-maxz) Goto 700
           Call gradu(ipot, ix, iy, iz, gradx, grady, gradz)
           p(1, i) = p(1, i) - (0.5*dt)*gradx
           p(2, i) = p(2, i) - (0.5*dt)*grady
           p(3, i) = p(3, i) - (0.5*dt)*gradz
700     End Do
     End If
     rcnne = 0
     rdd = 0
     rpp = 0
     rppk = 0
     rpn = 0
     rpd = 0
     rkn = 0
     rnnk = 0
     rddk = 0
     rndk = 0
     rcnnd = 0
     rcndn = 0
     rcoll = 0
     rbloc = 0
     rdirt = 0
     rdecay = 0
     rres = 0
     Do kkk = 1, 5
        skaon(kkk) = 0
        Do is = 1, 2000
           sekaon(kkk, is) = 0
        End Do
     End Do
     pr0 = 0.
     pr1 = 0.
     ska0 = 0.
     ska1 = 0.
     If (iapar2(1)/=1) Then
        Do i = 1, maxstr
           Do j = 1, 3
              r(j, i) = 0.
              p(j, i) = 0.
           End Do
           e(i) = 0.
           lb(i) = 0
           id(i) = 0
           proper(i) = 1.
        End Do
        mass = 0
        np = 0
        Do j = 1, num
           massr(j) = 0
           npi(j) = 1
        End Do
        Do i = 1, maxr
           Do j = 1, maxstr
              rt(1, j, i) = 0.
              rt(2, j, i) = 0.
              rt(3, j, i) = 0.
              pt(1, j, i) = 0.
              pt(2, j, i) = 0.
              pt(3, j, i) = 0.
              et(j, i) = 0.
              lt(j, i) = 0
              prot(j, i) = 1.
           End Do
        End Do
     End If
     Do nt = 1, ntmax
        lp1 = 0
        lp2 = 0
        lp3 = 0
        ld1 = 0
        ld2 = 0
        ld3 = 0
        ld4 = 0
        ln1 = 0
        ln2 = 0
        ln5 = 0
        le = 0
        lkaon = 0
        lkaons = 0
        If (icoll/=1) Then
           numnt = nt
           Call relcol(lcoll, lbloc, lcnne, ldd, lpp, lppk, lpn, lpd, lrho, lomega, lkn, lnnk, lddk, lndk, lcnnd, lcndn, ldirt, ldecay, lres, ldou, lddrho, lnnrho, lnnom, numnt, ntmax, sp, akaon, sk)
           rcoll = rcoll + float(lcoll)/num
           rbloc = rbloc + float(lbloc)/num
           rcnne = rcnne + float(lcnne)/num
           rdd = rdd + float(ldd)/num
           rpp = rpp + float(lpp)/num
           rppk = rppk + float(lppk)/num
           rpn = rpn + float(lpn)/num
           rpd = rpd + float(lpd)/num
           rkn = rkn + float(lkn)/num
           rnnk = rnnk + float(lnnk)/num
           rddk = rddk + float(lddk)/num
           rndk = rndk + float(lndk)/num
           rcnnd = rcnnd + float(lcnnd)/num
           rcndn = rcndn + float(lcndn)/num
           rdirt = rdirt + float(ldirt)/num
           rdecay = rdecay + float(ldecay)/num
           rres = rres + float(lres)/num
           adirt = ldirt/dt/num
           acoll = (lcoll-lbloc)/dt/num
           acnnd = lcnnd/dt/num
           acndn = lcndn/dt/num
           adecay = ldecay/dt/num
           ares = lres/dt/num
           adou = ldou/dt/num
           addrho = lddrho/dt/num
           annrho = lnnrho/dt/num
           annom = lnnom/dt/num
           add = ldd/dt/num
           app = lpp/dt/num
           appk = lppk/dt/num
           apn = lpn/dt/num
           apd = lpd/dt/num
           arh = lrho/dt/num
           aom = lomega/dt/num
           akn = lkn/dt/num
           annk = lnnk/dt/num
           addk = lddk/dt/num
           andk = lndk/dt/num
        End If
        Call dens(ipot, mass, num, outpar)
        sumene = 0
        iso = 0
        Do mrun = 1, num
           iso = iso + massr(mrun-1)
           Do i0 = 1, massr(mrun)
              i = i0 + iso
              etotal = sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
              sumene = sumene + etotal
              If (kpoten/=0 .And. lb(i)==23) Then
                 den = 0.
                 ix = nint(r(1,i))
                 iy = nint(r(2,i))
                 iz = nint(r(3,i))
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) den = rho(ix, iy, iz)
                 akg = 0.1727
                 bkg = 0.333
                 rnsg = den
                 ecor = -akg*rnsg + (bkg*den)**2
                 etotal = sqrt(etotal**2+ecor)
              End If
              If (kpoten/=0 .And. lb(i)==21) Then
                 den = 0.
                 ix = nint(r(1,i))
                 iy = nint(r(2,i))
                 iz = nint(r(3,i))
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) den = rho(ix, iy, iz)
                 akg = 0.1727
                 bkg = 0.333
                 rnsg = den
                 ecor = -akg*rnsg + (bkg*den)**2
                 etotal = sqrt(etotal**2+ecor)
              End If
              r(1, i) = r(1, i) + dt*p(1, i)/etotal
              r(2, i) = r(2, i) + dt*p(2, i)/etotal
              r(3, i) = r(3, i) + dt*p(3, i)/etotal
              If (cycbox/=0) Then
                 If (r(1,i)>cycbox/2) r(1, i) = r(1, i) - cycbox
                 If (r(1,i)<=-cycbox/2) r(1, i) = r(1, i) + cycbox
                 If (r(2,i)>cycbox/2) r(2, i) = r(2, i) - cycbox
                 If (r(2,i)<=-cycbox/2) r(2, i) = r(2, i) + cycbox
                 If (r(3,i)>cycbox/2) r(3, i) = r(3, i) - cycbox
                 If (r(3,i)<=-cycbox/2) r(3, i) = r(3, i) + cycbox
              End If
              lb1 = lb(i)
              If (lb1==9) ld1 = ld1 + 1
              If (lb1==8) ld2 = ld2 + 1
              If (lb1==7) ld3 = ld3 + 1
              If (lb1==6) ld4 = ld4 + 1
              If (lb1==11) ln1 = ln1 + 1
              If (lb1==10) ln2 = ln2 + 1
              If ((lb1==13) .Or. (lb1==12)) ln5 = ln5 + 1
              If (lb1==0) le = le + 1
              If (lb1==23) lkaon = lkaon + 1
              If (lb1==30) lkaons = lkaons + 1
              If (lb1==5) lp1 = lp1 + 1
              If (lb1==4) lp2 = lp2 + 1
              If (lb1==3) lp3 = lp3 + 1
           End Do
        End Do
        lp = lp1 + lp2 + lp3
        ld = ld1 + ld2 + ld3 + ld4
        ln = ln1 + ln2
        alp = float(lp)/float(num)
        ald = float(ld)/float(num)
        aln = float(ln)/float(num)
        aln5 = float(ln5)/float(num)
        atotal = alp + ald + aln + 0.5*aln5
        ale = float(le)/float(num)
        alkaon = float(lkaon)/float(num)
        If (icou==1) Then
           iso = 0
           Do irun = 1, num
              iso = iso + massr(irun-1)
              Do il = 1, massr(irun)
                 temp(1, il) = 0.
                 temp(2, il) = 0.
                 temp(3, il) = 0.
              End Do
              Do il = 1, massr(irun)
                 i = iso + il
                 If (zet(lb(i))/=0) Then
                    Do jl = 1, il - 1
                       j = iso + jl
                       If (zet(lb(j))/=0) Then
                          ddx = r(1, i) - r(1, j)
                          ddy = r(2, i) - r(2, j)
                          ddz = r(3, i) - r(3, j)
                          rdiff = sqrt(ddx**2+ddy**2+ddz**2)
                          If (rdiff<=1.) rdiff = 1.
                          grp = zet(lb(i))*zet(lb(j))/rdiff**3
                          ddx = ddx*grp
                          ddy = ddy*grp
                          ddz = ddz*grp
                          temp(1, il) = temp(1, il) + ddx
                          temp(2, il) = temp(2, il) + ddy
                          temp(3, il) = temp(3, il) + ddz
                          temp(1, jl) = temp(1, jl) - ddx
                          temp(2, jl) = temp(2, jl) - ddy
                          temp(3, jl) = temp(3, jl) - ddz
                       End If
                    End Do
                 End If
              End Do
              Do il = 1, massr(irun)
                 i = iso + il
                 If (zet(lb(i))/=0) Then
                    Do idir = 1, 3
                       p(idir, i) = p(idir, i) + temp(idir, il)*dt*0.00144
                    End Do
                 End If
              End Do
           End Do
        End If
        spt = 0
        spz = 0
        ncen = 0
        ekin = 0
        nlost = 0
        mean = 0
        nquark = 0
        nbaryn = 0
        rads = 2.
        zras = 0.1
        denst = 0.
        edenst = 0.
        Do irun = 1, num
           mean = mean + massr(irun-1)
           Do j = 1, massr(irun)
              i = j + mean
              radut = sqrt(r(1,i)**2+r(2,i)**2)
              If (radut<=rads) Then
                 If (abs(r(3,i))<=zras*nt*dt) Then
                    vols = 3.14159*rads**2*zras
                    engs = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
                    gammas = 1.
                    If (e(i)/=0.) gammas = engs/e(i)
                    denst = denst + 1./gammas/vols
                    edenst = edenst + engs/gammas/gammas/vols
                 End If
              End If
              drr = sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
              If (drr<=2.0) Then
                 spt = spt + p(1, i)**2 + p(2, i)**2
                 spz = spz + p(3, i)**2
                 ncen = ncen + 1
                 ekin = ekin + sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2) - e(i)
              End If
              ix = nint(r(1,i))
              iy = nint(r(2,i))
              iz = nint(r(3,i))
              If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                 If (rho(ix,iy,iz)/0.168>dencut) Goto 5800
                 If ((rho(ix,iy,iz)/0.168>5.) .And. (e(i)>0.9)) nbaryn = nbaryn + 1
                 If (pel(ix,iy,iz)>2.0) nquark = nquark + 1
              End If
              If (kpoten/=0 .And. lb(i)==23) Then
                 den = 0.
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    den = rho(ix, iy, iz)
                    akg = 0.1727
                    bkg = 0.333
                    rnsg = den
                    ecor = -akg*rnsg + (bkg*den)**2
                    etotal = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2+ecor)
                    ecor = -akg + 2.*bkg**2*den + 2.*bkg*etotal
                    Call graduk(ix, iy, iz, gradxk, gradyk, gradzk)
                    p(1, i) = p(1, i) - dt*gradxk*ecor/(2.*etotal)
                    p(2, i) = p(2, i) - dt*gradyk*ecor/(2.*etotal)
                    p(3, i) = p(3, i) - dt*gradzk*ecor/(2.*etotal)
                 End If
              End If
              If (kpoten/=0 .And. lb(i)==21) Then
                 den = 0.
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    den = rho(ix, iy, iz)
                    Call graduk(ix, iy, iz, gradxk, gradyk, gradzk)
                    akg = 0.1727
                    bkg = 0.333
                    rnsg = den
                    ecor = -akg*rnsg + (bkg*den)**2
                    etotal = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2+ecor)
                    ecor = -akg + 2.*bkg**2*den - 2.*bkg*etotal
                    p(1, i) = p(1, i) - dt*gradxk*ecor/(2.*etotal)
                    p(2, i) = p(2, i) - dt*gradyk*ecor/(2.*etotal)
                    p(3, i) = p(3, i) - dt*gradzk*ecor/(2.*etotal)
                 End If
              End If
              If (j>mass) Goto 5800
              If (icoll/=-1) Then
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    Call gradu(ipot, ix, iy, iz, gradx, grady, gradz)
                    tz = 0.
                    gradxn = 0
                    gradyn = 0
                    gradzn = 0
                    gradxp = 0
                    gradyp = 0
                    gradzp = 0
                    If (icou==1) Then
                       Call gradup(ix, iy, iz, gradxp, gradyp, gradzp)
                       Call gradun(ix, iy, iz, gradxn, gradyn, gradzn)
                       If (zet(lb(i))/=0) tz = -1
                       If (zet(lb(i))==0) tz = 1
                    End If
                    If (iabs(lb(i))>=14 .And. iabs(lb(i))<=17) Then
                       facl = 2./3.
                    Else If (iabs(lb(i))==40 .Or. iabs(lb(i))==41) Then
                       facl = 1./3.
                    Else
                       facl = 1.
                    End If
                    p(1, i) = p(1, i) - facl*dt*(gradx+asy*(gradxn-gradxp)*tz)
                    p(2, i) = p(2, i) - facl*dt*(grady+asy*(gradyn-gradyp)*tz)
                    p(3, i) = p(3, i) - facl*dt*(gradz+asy*(gradzn-gradzp)*tz)
                 End If
              End If
5800       End Do
        End Do
        cden = rho(0, 0, 0)/0.168
        If ((nt/nfreq)*nfreq==nt) Then
           If (icflow==1) Call flow(nt)
        End If
        If (iapar2(1)/=1) Then
           ct = nt*dt
           ia = 0
           Do irun = 1, num
              Do ic = 1, massr(irun)
                 ie = ia + ic
                 rt(1, ic, irun) = r(1, ie)
                 rt(2, ic, irun) = r(2, ie)
                 rt(3, ic, irun) = r(3, ie)
                 pt(1, ic, irun) = p(1, ie)
                 pt(2, ic, irun) = p(2, ie)
                 pt(3, ic, irun) = p(3, ie)
                 et(ic, irun) = e(ie)
                 lt(ic, irun) = lb(ie)
                 prot(ic, irun) = proper(ie)
                 dpertt(ic, irun) = dpertp(ie)
              End Do
              np = massr(irun)
              np1 = npi(irun)
              ctlong = ct
              If (nt==(ntmax-1)) Then
                 ctlong = 1.E30
              Else If (nt==ntmax) Then
                 Goto 1111
              End If
              Do While (np1<=multi1(irun) .And. ft1(np1,irun)>((nt-1)*dt) .And. ft1(np1,irun)<=ctlong)
                 np = np + 1
                 udt = (ct-ft1(np1,irun))/ee1(np1, irun)
                 If (nt==(ntmax-1)) Then
                    ftsvt(np, irun) = ft1(np1, irun)
                    If (ft1(np1,irun)>ct) udt = 0.
                 End If
                 rt(1, np, irun) = gx1(np1, irun) + px1(np1, irun)*udt
                 rt(2, np, irun) = gy1(np1, irun) + py1(np1, irun)*udt
                 rt(3, np, irun) = gz1(np1, irun) + pz1(np1, irun)*udt
                 pt(1, np, irun) = px1(np1, irun)
                 pt(2, np, irun) = py1(np1, irun)
                 pt(3, np, irun) = pz1(np1, irun)
                 et(np, irun) = xm1(np1, irun)
                 lt(np, irun) = iarflv(ityp1(np1,irun))
                 dpertt(np, irun) = dpp1(np1, irun)
                 np1 = np1 + 1
                 prot(np, irun) = 1.
              End Do
1111          Continue
              npi(irun) = np1
              ia = ia + massr(irun)
              massr(irun) = np
           End Do
           ia = 0
           Do irun = 1, num
              ia = ia + massr(irun-1)
              Do ic = 1, massr(irun)
                 ie = ia + ic
                 r(1, ie) = rt(1, ic, irun)
                 r(2, ie) = rt(2, ic, irun)
                 r(3, ie) = rt(3, ic, irun)
                 p(1, ie) = pt(1, ic, irun)
                 p(2, ie) = pt(2, ic, irun)
                 p(3, ie) = pt(3, ic, irun)
                 e(ie) = et(ic, irun)
                 lb(ie) = lt(ic, irun)
                 proper(ie) = prot(ic, irun)
                 If (nt==(ntmax-1)) ftsv(ie) = ftsvt(ic, irun)
                 dpertp(ie) = dpertt(ic, irun)
              End Do
              Call hbtout(massr(irun), nt, ntmax)
           End Do
        End If
     End Do
     iss = 0
     Do lrun = 1, num
        iss = iss + massr(lrun-1)
        Do l0 = 1, massr(lrun)
           ipart = iss + l0
        End Do
     End Do
     If (iapar2(1)/=1) Then
        ia = 0
        Do irun = 1, num
           ia = ia + massr(irun-1)
           np1 = npi(irun)
           nsh = massr(irun) - np1 + 1
           multi1(irun) = multi1(irun) + nsh
           If (nsh>0) Then
              ib = multi1(irun)
              ie = massr(irun) + 1
              ii = -1
           Else If (nsh<0) Then
              ib = massr(irun) + 1
              ie = multi1(irun)
              ii = 1
           End If
           If (nsh/=0) Then
              Do i = ib, ie, ii
                 j = i - nsh
                 ityp1(i, irun) = ityp1(j, irun)
                 gx1(i, irun) = gx1(j, irun)
                 gy1(i, irun) = gy1(j, irun)
                 gz1(i, irun) = gz1(j, irun)
                 ft1(i, irun) = ft1(j, irun)
                 px1(i, irun) = px1(j, irun)
                 py1(i, irun) = py1(j, irun)
                 pz1(i, irun) = pz1(j, irun)
                 ee1(i, irun) = ee1(j, irun)
                 xm1(i, irun) = xm1(j, irun)
                 pro1(i, irun) = pro1(j, irun)
                 dpp1(i, irun) = dpp1(j, irun)
              End Do
           End If
           Do i = 1, massr(irun)
              ib = ia + i
              ityp1(i, irun) = invflv(lb(ib))
              gx1(i, irun) = r(1, ib)
              gy1(i, irun) = r(2, ib)
              gz1(i, irun) = r(3, ib)
              If (ft1(i,irun)<ct) ft1(i, irun) = ct
              px1(i, irun) = p(1, ib)
              py1(i, irun) = p(2, ib)
              pz1(i, irun) = p(3, ib)
              xm1(i, irun) = e(ib)
              ee1(i, irun) = sqrt(px1(i,irun)**2+py1(i,irun)**2+pz1(i,irun)**2+xm1(i,irun)**2)
              pro1(i, irun) = proper(ib)
           End Do
        End Do
     End If
  End Do
  Return
End Subroutine artmn
