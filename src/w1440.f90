Real Function w1440(dmass)
  Save
  avmass = 0.938868
  pimass = 0.137265
  aux = 0.25*(dmass**2-avmass**2-pimass**2)**2 - (avmass*pimass)**2
  If (aux>0.) Then
     qavail = sqrt(aux)/dmass
  Else
     qavail = 1.E-06
  End If
  w1440 = 0.2*(qavail/0.397)**3
  Return
End Function w1440
