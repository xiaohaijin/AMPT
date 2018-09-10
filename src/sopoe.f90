Subroutine sopoe(lb1, lb2, srt)
  Parameter (etam=0.5475, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  xopoe = 0.
  If ((lb1==28 .And. lb2>=3 .And. lb2<=5) .Or. (lb2==28 .And. lb1>=3 .And. lb1<=5)) Then
    If (srt>(aomega+etam)) xopoe = xop2oe(srt)
  Else If ((lb1==28 .And. lb2==0) .Or. (lb1==0 .And. lb2==28)) Then
    If (srt>(aomega+etam)) xopoe = xoe2op(srt)
  End If
  Return
End Subroutine sopoe
