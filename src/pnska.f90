Real Function pnska(srt)
  Save
  If (srt>3.0) Then
     pnska = 0
     Return
  End If
  ala = 1.116
  aka = 0.498
  ana = 0.939
  asa = 1.197
  t1 = asa + aka
  If (srt<=t1) Then
     pnska = 0
     Return
  End If
  If (srt<1.9) sbb1 = (0.7/0.218)*(srt-t1)
  If (srt>=1.9) sbb1 = 0.14/(srt-1.7)
  sbb2 = 0.
  If (srt>1.682) sbb2 = 0.5*(1.-0.75*(srt-1.682))
  pnska = 0.25*(sbb1+sbb2)
  pnska = pnska/10.
  Return
End Function pnska
