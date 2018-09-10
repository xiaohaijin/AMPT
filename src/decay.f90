Subroutine decay(irun, i, nnn, iseed, wid, nt)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, etam=0.5475, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  lbanti = lb(i)
  dm = e(i)
  If (iabs(lb(i))==11) Then
     x3 = ranart(nseed)
     If (x3>(1./3.)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
     Else
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
     End If
  Else If (iabs(lb(i))==10) Then
     x4 = ranart(nseed)
     If (x4>(1./3.)) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
     Else
        lb(i) = 2
        nalb = 2
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
     End If
  Else If (iabs(lb(i))==12) Then
     ctrl = 0.65
     If (dm<=1.49) ctrl = -1.
     x5 = ranart(nseed)
     If (x5>=ctrl) Then
        x6 = ranart(nseed)
        If (x6>(1./3.)) Then
           lb(i) = 1
           nlab = 1
           lpion(nnn, irun) = 3
           epion(nnn, irun) = ap2
        Else
           lb(i) = 2
           nalb = 2
           lpion(nnn, irun) = 4
           epion(nnn, irun) = ap1
        End If
     Else
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 0
        epion(nnn, irun) = etam
     End If
  Else If (iabs(lb(i))==13) Then
     ctrl = 0.65
     If (dm<=1.49) ctrl = -1.
     x5 = ranart(nseed)
     If (x5>=ctrl) Then
        x8 = ranart(nseed)
        If (x8>(1./3.)) Then
           lb(i) = 2
           nlab = 2
           lpion(nnn, irun) = 5
           epion(nnn, irun) = ap2
        Else
           lb(i) = 1
           nlab = 1
           lpion(nnn, irun) = 4
           epion(nnn, irun) = ap1
        End If
     Else
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 0
        epion(nnn, irun) = etam
     End If
  End If
  Call dkine(irun, i, nnn, nlab, iseed, wid, nt)
  If (lbanti<0) Then
     lbi = lb(i)
     If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     Else If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     End If
     lb(i) = lbi
     lbi = lpion(nnn, irun)
     If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     Else If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     End If
     lpion(nnn, irun) = lbi
  End If
  If (nt==ntmax) Then
     lbm = lpion(nnn, irun)
     If (lbm==0 .Or. lbm==25 .Or. lbm==26 .Or. lbm==27) Then
        lbsave = lbm
        xmsave = epion(nnn, irun)
        pxsave = ppion(1, nnn, irun)
        pysave = ppion(2, nnn, irun)
        pzsave = ppion(3, nnn, irun)
        dpsave = dppion(nnn, irun)
        lpion(nnn, irun) = lb(i)
        epion(nnn, irun) = e(i)
        ppion(1, nnn, irun) = p(1, i)
        ppion(2, nnn, irun) = p(2, i)
        ppion(3, nnn, irun) = p(3, i)
        dppion(nnn, irun) = dpertp(i)
        lb(i) = lbsave
        e(i) = xmsave
        p(1, i) = pxsave
        p(2, i) = pysave
        p(3, i) = pzsave
        dpertp(i) = dpsave
     End If
  End If
  Return
End Subroutine decay
