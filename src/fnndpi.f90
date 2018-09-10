Function fnndpi(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fnndpi = 0.
  Else
    fnndpi = 26.*exp(-(s-4.65)**2/0.1) + 4.*exp(-(s-4.65)**2/2.) + 0.28*exp(-(s-6.)**2/10.)
  End If
  Return
End Function fnndpi
