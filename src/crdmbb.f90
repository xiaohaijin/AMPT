Subroutine crdmbb(px, py, pz, srt, i1, i2, iblock, ntag, sig, nt, ianti)
  Parameter (maxstr=150001, maxr=1)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /arevt/iaevt, iarun, miss
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /dpifsl/lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1, lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2, lbpp1, lbpp2
  Common /dpifsm/xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2, xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1, xmpp2
  Common /dpisig/sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp, sdmpp
  Common /rndf77/nseed
  Save
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  s = srt**2
  If (sig<=0) Return
  If (iabs(lb1)==42) Then
    ideut = i1
    lbm = lb2
    idm = i2
  Else
    ideut = i2
    lbm = lb1
    idm = i1
  End If
  If ((idpert==1 .Or. idpert==2) .And. dpertp(ideut)/=1.) Then
    x1 = ranart(nseed)
    If (x1<=sdmel/sig) Then
      If (ianti==0) Then
        Write (91, *) '  d+', lbm, ' (pert d M elastic) @nt=', nt, ' @prob=', dpertp(ideut)
      Else
        Write (91, *) '  d+', lbm, ' (pert dbar M elastic) @nt=', nt, ' @prob=', dpertp(ideut)
      End If
      scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
      If (scheck<0) Then
        Write (99, *) 'scheck51: ', scheck
        scheck = 0.
      End If
      pfinal = sqrt(scheck)/2./srt
      Call dmelangle(pxn, pyn, pzn, pfinal)
      Call rotate(px, py, pz, pxn, pyn, pzn)
      edcm = sqrt(e(ideut)**2+pxn**2+pyn**2+pzn**2)
      pdbeta = pxn*betax + pyn*betay + pzn*betaz
      transf = gamma*(gamma*pdbeta/(gamma+1.)+edcm)
      pt1d = betax*transf + pxn
      pt2d = betay*transf + pyn
      pt3d = betaz*transf + pzn
      p(1, ideut) = pt1d
      p(2, ideut) = pt2d
      p(3, ideut) = pt3d
      iblock = 504
      px1 = p(1, i1)
      py1 = p(2, i1)
      pz1 = p(3, i1)
      id(i1) = 2
      id(i2) = 2
      r(1, ideut) = r(1, idm)
      r(2, ideut) = r(2, idm)
      r(3, ideut) = r(3, idm)
    Else
      If (ianti==0) Then
        Write (91, *) '  d+', lbm, ' ->BB (pert d destrn) @nt=', nt, ' @prob=', dpertp(ideut)
      Else
        Write (91, *) '  d+', lbm, ' ->BB (pert dbar destrn) @nt=', nt, ' @prob=', dpertp(ideut)
      End If
      e(ideut) = 0.
      iblock = 502
    End If
    Return
  End If
  iblock = 502
  x1 = ranart(nseed)
  If (x1<=sdmnn/sig) Then
    lbb1 = lbnn1
    lbb2 = lbnn2
    xmb1 = xmnn1
    xmb2 = xmnn2
  Else If (x1<=(sdmnn+sdmnd)/sig) Then
    lbb1 = lbnd1
    lbb2 = lbnd2
    xmb1 = xmnd1
    xmb2 = xmnd2
  Else If (x1<=(sdmnn+sdmnd+sdmns)/sig) Then
    lbb1 = lbns1
    lbb2 = lbns2
    xmb1 = xmns1
    xmb2 = xmns2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp)/sig) Then
    lbb1 = lbnp1
    lbb2 = lbnp2
    xmb1 = xmnp1
    xmb2 = xmnp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd)/sig) Then
    lbb1 = lbdd1
    lbb2 = lbdd2
    xmb1 = xmdd1
    xmb2 = xmdd2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds)/sig) Then
    lbb1 = lbds1
    lbb2 = lbds2
    xmb1 = xmds1
    xmb2 = xmds2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp)/sig) Then
    lbb1 = lbdp1
    lbb2 = lbdp2
    xmb1 = xmdp1
    xmb2 = xmdp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss)/sig) Then
    lbb1 = lbss1
    lbb2 = lbss2
    xmb1 = xmss1
    xmb2 = xmss2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss+sdmsp)/sig) Then
    lbb1 = lbsp1
    lbb2 = lbsp2
    xmb1 = xmsp1
    xmb2 = xmsp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss+sdmsp+sdmpp)/sig) Then
    lbb1 = lbpp1
    lbb2 = lbpp2
    xmb1 = xmpp1
    xmb2 = xmpp2
  Else
    lbb1 = lb1
    lbb2 = lb2
    xmb1 = em1
    xmb2 = em2
    iblock = 504
  End If
  lb(i1) = lbb1
  e(i1) = xmb1
  lb(i2) = lbb2
  e(i2) = xmb2
  lb1 = lb(i1)
  lb2 = lb(i2)
  scheck = (s-(xmb1+xmb2)**2)*(s-(xmb1-xmb2)**2)
  If (scheck<0) Then
    Write (99, *) 'scheck52: ', scheck
    scheck = 0.
  End If
  pfinal = sqrt(scheck)/2./srt
  If (iblock==502) Then
    Call dmangle(pxn, pyn, pzn, nt, ianti, pfinal, lbm)
  Else If (iblock==504) Then
    If (ianti==0) Then
      Write (91, *) ' d+', lbm, ' (regular d M elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
    Else
      Write (91, *) ' d+', lbm, ' (regular dbar M elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
    End If
    Call dmelangle(pxn, pyn, pzn, pfinal)
  Else
    Print *, 'Wrong iblock number in crdmbb()'
    Stop
  End If
  Call rotate(px, py, pz, pxn, pyn, pzn)
  e1cm = sqrt(e(i1)**2+pxn**2+pyn**2+pzn**2)
  p1beta = pxn*betax + pyn*betay + pzn*betaz
  transf = gamma*(gamma*p1beta/(gamma+1.)+e1cm)
  pt1i1 = betax*transf + pxn
  pt2i1 = betay*transf + pyn
  pt3i1 = betaz*transf + pzn
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e2cm = sqrt(e(i2)**2+pxn**2+pyn**2+pzn**2)
  p2beta = -pxn*betax - pyn*betay - pzn*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf - pxn
  pt2i2 = betay*transf - pyn
  pt3i2 = betaz*transf - pzn
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  em2 = e(i2)
  id(i1) = 2
  id(i2) = 2
  Return
End Subroutine crdmbb
