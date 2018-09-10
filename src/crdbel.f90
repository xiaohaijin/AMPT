Subroutine crdbel(px, py, pz, srt, i1, i2, iblock, ntag, sig, nt, ianti)
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
  Save
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  s = srt**2
  If (sig<=0) Return
  iblock = 503
  If (iabs(lb1)==42) Then
    ideut = i1
    lbb = lb2
    idb = i2
  Else
    ideut = i2
    lbb = lb1
    idb = i1
  End If
  If ((idpert==1 .Or. idpert==2) .And. dpertp(ideut)/=1.) Then
    If (ianti==0) Then
      Write (91, *) '  d+', lbb, ' (pert d B elastic) @nt=', nt, ' @prob=', dpertp(ideut), p(1, idb), p(2, idb), p(1, ideut), p(2, ideut)
    Else
      Write (91, *) '  d+', lbb, ' (pert dbar Bbar elastic) @nt=', nt, ' @prob=', dpertp(ideut), p(1, idb), p(2, idb), p(1, ideut), p(2, ideut)
    End If
    scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
    If (scheck<0) Then
      Write (99, *) 'scheck53: ', scheck
      scheck = 0.
    End If
    pfinal = sqrt(scheck)/2./srt
    Call dbelangle(pxn, pyn, pzn, pfinal)
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
    px1 = p(1, i1)
    py1 = p(2, i1)
    pz1 = p(3, i1)
    id(i1) = 2
    id(i2) = 2
    r(1, ideut) = r(1, idb)
    r(2, ideut) = r(2, idb)
    r(3, ideut) = r(3, idb)
    Return
  End If
  If (ianti==0) Then
    Write (91, *) ' d+', lbb, ' (regular d B elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  Else
    Write (91, *) ' d+', lbb, ' (regular dbar Bbar elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  End If
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<0) Then
    Write (99, *) 'scheck54: ', scheck
    scheck = 0.
  End If
  pfinal = sqrt(scheck)/2./srt
  Call dbelangle(pxn, pyn, pzn, pfinal)
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
End Subroutine crdbel
