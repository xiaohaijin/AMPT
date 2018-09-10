Subroutine srree(lb1, lb2, srt)
  Parameter (etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  rree = 0.
  If (lb1>=25 .And. lb1<=27 .And. lb2>=25 .And. lb2<=27) Then
    If (srt>(2*etam)) rree = rrtoee(srt)
  Else If (lb1==0 .And. lb2==0) Then
    If (srt>(2*arho)) rree = eetorr(srt)
  End If
  Return
End Subroutine srree
