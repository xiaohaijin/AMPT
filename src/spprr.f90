Subroutine spprr(lb1, lb2, srt)
  Parameter (arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  pprr = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
    If (srt>(2*arho)) pprr = ptor(srt)
  Else If ((lb1>=25 .And. lb1<=27) .And. (lb2>=25 .And. lb2<=27)) Then
    pprr = rtop(srt)
  End If
  Return
End Subroutine spprr
