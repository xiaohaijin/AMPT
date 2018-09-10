Subroutine srpre(lb1, lb2, srt)
  Parameter (pimass=0.140, etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  rpre = 0.
  If (lb1>=25 .And. lb1<=27 .And. lb2>=3 .And. lb2<=5) Then
    If (srt>(etam+arho)) rpre = rptore(srt)
  Else If (lb2>=25 .And. lb2<=27 .And. lb1>=3 .And. lb1<=5) Then
    If (srt>(etam+arho)) rpre = rptore(srt)
  Else If (lb1>=25 .And. lb1<=27 .And. lb2==0) Then
    If (srt>(pimass+arho)) rpre = retorp(srt)
  Else If (lb2>=25 .And. lb2<=27 .And. lb1==0) Then
    If (srt>(pimass+arho)) rpre = retorp(srt)
  End If
  Return
End Subroutine srpre
