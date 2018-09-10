Subroutine bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
  Parameter (pi=3.1415926)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /rndf77/nseed
  Common /para8/idpert, npertd, idxsec
  Common /arevt/iaevt, iarun, miss
  Save
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pzd = pfinal*c1
  pxd = pfinal*s1*ct1
  pyd = pfinal*s1*st1
  If (idpert==1 .And. npertd>=1) Then
    dprob = dprob1
  Else If (idpert==2 .And. npertd>=1) Then
    dprob = 1./float(npertd)
  End If
  If (ianti==0) Then
    If (idpert==0 .Or. (idpert==1 .And. ipert1==0) .Or. (idpert==2 .And. idloop==(npertd+1))) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (regular d prodn)    @evt#', iaevt, ' @nt=', nt
    Else If ((idpert==1 .Or. idpert==2) .And. idloop==npertd) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (pert d prodn)       @evt#', iaevt, ' @nt=', nt, ' @prob=', dprob
    End If
  Else
    If (idpert==0 .Or. (idpert==1 .And. ipert1==0) .Or. (idpert==2 .And. idloop==(npertd+1))) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (regular dbar prodn) @evt#', iaevt, ' @nt=', nt
    Else If ((idpert==1 .Or. idpert==2) .And. idloop==npertd) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (pert dbar prodn)    @evt#', iaevt, ' @nt=', nt, ' @prob=', dprob
    End If
  End If
  Return
End Subroutine bbdangle
