Subroutine spppe(lb1, lb2, srt)
  Parameter (pimass=0.140, etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  pppe = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
    If (srt>(etam+pimass)) pppe = pptope(srt)
  Else If ((lb1>=3 .And. lb1<=5) .And. lb2==0) Then
    pppe = petopp(srt)
  Else If ((lb2>=3 .And. lb2<=5) .And. lb1==0) Then
    pppe = petopp(srt)
  End If
  Return
End Subroutine spppe
