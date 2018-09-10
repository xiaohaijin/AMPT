Function fdbel(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fdbel = 0.
  Else
    fdbel = 2500.*exp(-(s-7.93)**2/0.003) + 300.*exp(-(s-7.93)**2/0.1) + 10.
  End If
  Return
End Function fdbel
