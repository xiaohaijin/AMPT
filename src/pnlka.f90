Real Function pnlka(srt)
  Save
  ala = 1.116
  aka = 0.498
  ana = 0.939
  t1 = ala + aka
  If (srt<=t1) Then
     pnlka = 0
  Else
     If (srt<1.7) sbbk = (0.9/0.091)*(srt-t1)
     If (srt>=1.7) sbbk = 0.09/(srt-1.6)
     pnlka = 0.25*sbbk
     pnlka = pnlka/10.
  End If
  Return
End Function pnlka
