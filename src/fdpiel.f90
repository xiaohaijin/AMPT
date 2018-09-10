Function fdpiel(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fdpiel = 0.
  Else
    fdpiel = 63.*exp(-(s-4.67)**2/0.15) + 15.*exp(-(s-6.25)**2/0.3)
  End If
  Return
End Function fdpiel
