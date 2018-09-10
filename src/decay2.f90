Subroutine decay2(irun, i, nnn, iseed, wid, nt)
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
  Common /rndf77/nseed
  Save
  lbanti = lb(i)
  dm = e(i)
  If (iabs(lb(i))==11) Then
     x3 = ranart(nseed)
     If (x3<(1./3)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else If (x3<2./3 .And. x3>1./3.) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 3
        epion(nnn+1, irun) = ap2
     Else
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     End If
  Else If (iabs(lb(i))==10) Then
     x3 = ranart(nseed)
     If (x3<(1./3)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else If (x3<2./3 .And. x3>1./3.) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 3
        epion(nnn+1, irun) = ap2
     End If
  End If
  Call dkine2(irun, i, nnn, nlab, iseed, wid, nt)
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
     lbi = lpion(nnn+1, irun)
     If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     Else If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     End If
     lpion(nnn+1, irun) = lbi
  End If
  Return
End Subroutine decay2
