Real Function xppbar(srt)
  Parameter (pmass=0.9383, xmax=400.)
  Save
  xppbar = 1.E-06
  plab2 = (srt**2/(2.*pmass)-pmass)**2 - pmass**2
  If (plab2>0) Then
    plab = sqrt(plab2)
    xppbar = 67./(plab**0.7)
    If (xppbar>xmax) xppbar = xmax
  End If
  Return
End Function xppbar
