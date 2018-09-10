Subroutine dmangle(pxn, pyn, pzn, nt, ianti, pfinal, lbm)
  Parameter (pi=3.1415926)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /arevt/iaevt, iarun, miss
  Common /rndf77/nseed
  Save
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pzn = pfinal*c1
  pxn = pfinal*s1*ct1
  pyn = pfinal*s1*st1
  If (ianti==0) Then
    Write (91, *) ' d+', lbm, ' ->BB (regular d destrn) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  Else
    Write (91, *) ' d+', lbm, ' ->BB (regular dbar destrn) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  End If
  Return
End Subroutine dmangle
