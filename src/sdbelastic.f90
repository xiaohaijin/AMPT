Subroutine sdbelastic(srt, sdb)
  Parameter (srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Save
  sdb = 0.
  sdbel = 0.
  If (srt<=(em1+em2)) Return
  s = srt**2
  If (idxsec==1 .Or. idxsec==3) Then
    sdbel = fdbel(s)
  Else If (idxsec==2 .Or. idxsec==4) Then
    threshold = em1 + em2
    snew = (srt-threshold+srt0)**2
    sdbel = fdbel(snew)
  End If
  sdb = sdbel
  Return
End Subroutine sdbelastic
