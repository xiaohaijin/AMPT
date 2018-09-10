Real Function xrho(srt)
  Save
  pmass = 0.9383
  rmass = 0.77
  trho = 0.151
  xrho = 0.000000001
  If (srt<=2.67) Return
  esmin = 2.*0.9383 + rmass - trho/2.
  es = srt
  xrho0 = 0.24*(es-esmin)/(1.4+(es-esmin)**2)
  xrho = 3.*xrho0
  Return
End Function xrho
