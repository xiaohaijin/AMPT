Subroutine resdec(i1, nt, nnn, wid, idecay, ipion)
  Parameter (hbarc=0.19733)
  Parameter (ak0=0.498, apich=0.140, api0=0.135, an=0.940, addm=0.02)
  Parameter (maxstr=150001, maxr=1)
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /resdcy/nsav, iksdcy
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  External iarflv, invflv
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
  irun = idecay
  If (nt==ntmax .And. ipi0dcy==1 .And. ((lb1==4 .And. ipion==0) .Or. ipion>=1)) Then
     kf = 111
  Else If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. lb1==24 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=9) .Or. iabs(lb1)==16) Then
     kf = invflv(lb1)
  Else
     Return
  End If
  ip = 1
  n = 1
  k(ip, 1) = 1
  k(ip, 3) = 0
  k(ip, 4) = 0
  k(ip, 5) = 0
  k(ip, 2) = kf
  If (ipion==0) Then
     p(ip, 1) = px1
     p(ip, 2) = py1
     p(ip, 3) = pz1
     If ((lb1==0 .Or. lb1==28) .And. em1<(2*apich+api0+addm)) Then
        em1 = 2*apich + api0 + addm
     Else If (lb1>=25 .And. lb1<=27 .And. em1<(2*apich+addm)) Then
        em1 = 2*apich + addm
     Else If (iabs(lb1)==30 .And. em1<(apich+ak0+addm)) Then
        em1 = apich + ak0 + addm
     Else If (iabs(lb1)>=6 .And. iabs(lb1)<=9 .And. em1<(apich+an+addm)) Then
        em1 = apich + an + addm
     End If
     e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
     p(ip, 4) = e1
     p(ip, 5) = em1
     dpdecp = dpertp(i1)
  Else If (nt==ntmax .And. ipi0dcy==1 .And. ipion>=1) Then
     p(ip, 1) = ppion(1, ipion, irun)
     p(ip, 2) = ppion(2, ipion, irun)
     p(ip, 3) = ppion(3, ipion, irun)
     p(ip, 5) = epion(ipion, irun)
     p(ip, 4) = sqrt(p(ip,5)**2+p(ip,1)**2+p(ip,2)**2+p(ip,3)**2)
     dpdecp = dppion(ipion, irun)
  Else
     Print *, 'stopped in resdec() a'
     Stop
  End If
  Call ludecy(ip)
  If (nt==ntmax) Then
     tau0 = hbarc/wid
     taudcy = tau0*(-1.)*alog(1.-ranart(nseed))
     ndaut = n - nsav
     If (ndaut<=1) Then
        Write (10, *) 'note: ndaut(<1)=', ndaut
        Call lulist(2)
        Stop
     End If
     If (ipion==0) Then
        taudcy = taudcy*e1/em1
        tfnl = tfnl + taudcy
        xfnl = xfnl + px1/e1*taudcy
        yfnl = yfnl + py1/e1*taudcy
        zfnl = zfnl + pz1/e1*taudcy
     Else If (ipion>=1) Then
        taudcy = taudcy*p(ip, 4)/p(ip, 5)
        tfnl = tfdpi(ipion, irun) + taudcy
        xfnl = rpion(1, ipion, irun) + p(ip, 1)/p(ip, 4)*taudcy
        yfnl = rpion(2, ipion, irun) + p(ip, 2)/p(ip, 4)*taudcy
        zfnl = rpion(3, ipion, irun) + p(ip, 3)/p(ip, 4)*taudcy
     Else
        Print *, 'stopped in resdec() b', ipion, wid, p(ip, 4)
        Stop
     End If
     If (n>=(nsav+2) .And. ipion==0) Then
        Do idau = nsav + 2, n
           kdaut = k(idau, 2)
           If (kdaut==221 .Or. kdaut==113 .Or. kdaut==213 .Or. kdaut==-213 .Or. kdaut==310) Then
              ksave = kdaut
              pxsave = p(idau, 1)
              pysave = p(idau, 2)
              pzsave = p(idau, 3)
              esave = p(idau, 4)
              xmsave = p(idau, 5)
              k(idau, 2) = k(nsav+1, 2)
              p(idau, 1) = p(nsav+1, 1)
              p(idau, 2) = p(nsav+1, 2)
              p(idau, 3) = p(nsav+1, 3)
              p(idau, 4) = p(nsav+1, 4)
              p(idau, 5) = p(nsav+1, 5)
              k(nsav+1, 2) = ksave
              p(nsav+1, 1) = pxsave
              p(nsav+1, 2) = pysave
              p(nsav+1, 3) = pzsave
              p(nsav+1, 4) = esave
              p(nsav+1, 5) = xmsave
              Goto 111
           End If
        End Do
     End If
111  Continue
     enet = 0.
     Do idau = nsav + 1, n
        enet = enet + p(idau, 4)
     End Do
  End If
  Do idau = nsav + 1, n
     kdaut = k(idau, 2)
     lbdaut = iarflv(kdaut)
     If (nt==ntmax .And. (kdaut==130 .Or. kdaut==310 .Or. iabs(kdaut)==311)) Then
        If (kdaut==130) Then
           lbdaut = 22
        Else If (kdaut==310) Then
           lbdaut = 24
        Else If (iabs(kdaut)==311) Then
           If (ranart(nseed)<0.5) Then
              lbdaut = 22
           Else
              lbdaut = 24
           End If
        End If
     End If
     If (idau==(nsav+1)) Then
        If (ipion==0) Then
           lb(i1) = lbdaut
           e(i1) = p(idau, 5)
           px1n = p(idau, 1)
           py1n = p(idau, 2)
           pz1n = p(idau, 3)
           dp1n = dpdecp
        Else If (ipion>=1) Then
           lpion(ipion, irun) = lbdaut
           epion(ipion, irun) = p(idau, 5)
           ppion(1, ipion, irun) = p(idau, 1)
           ppion(2, ipion, irun) = p(idau, 2)
           ppion(3, ipion, irun) = p(idau, 3)
           rpion(1, ipion, irun) = xfnl
           rpion(2, ipion, irun) = yfnl
           rpion(3, ipion, irun) = zfnl
           tfdpi(ipion, irun) = tfnl
           dppion(ipion, irun) = dpdecp
        End If
     Else
        nnn = nnn + 1
        lpion(nnn, irun) = lbdaut
        epion(nnn, irun) = p(idau, 5)
        ppion(1, nnn, irun) = p(idau, 1)
        ppion(2, nnn, irun) = p(idau, 2)
        ppion(3, nnn, irun) = p(idau, 3)
        rpion(1, nnn, irun) = xfnl
        rpion(2, nnn, irun) = yfnl
        rpion(3, nnn, irun) = zfnl
        tfdpi(nnn, irun) = tfnl
        dppion(nnn, irun) = dpdecp
     End If
  End Do
  Return
End Subroutine resdec
