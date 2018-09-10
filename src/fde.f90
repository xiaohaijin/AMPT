Real Function fde(dmass, srt, con)
  Save
  amn = 0.938869
  avpi = 0.13803333
  am0 = 1.232
  fd = 4.*(am0**2)*width(dmass)/((dmass**2-1.232**2)**2+am0**2*width(dmass)**2)
  If (con==1.) Then
     p11 = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (p11<=0) p11 = 1.E-06
     p1 = sqrt(p11)
  Else
     dmass = amn + avpi
     p11 = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (p11<=0) p11 = 1.E-06
     p1 = sqrt(p11)
  End If
  fde = fd*p1*dmass
  Return
End Function fde
