Subroutine sppee(lb1, lb2, srt)
  Parameter (etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  ppee = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
    If (srt>(2*etam)) ppee = ptoe(srt)
  Else If (lb1==0 .And. lb2==0) Then
    ppee = etop(srt)
  End If
  Return
End Subroutine sppee
