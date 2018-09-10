Real Function fd5(dmass, srt, con)
  Save
  amn = 0.938869
  avpi = 0.13803333
  am0 = 1.535
  fd = 4.*(am0**2)*w1535(dmass)/((dmass**2-1.535**2)**2+am0**2*w1535(dmass)**2)
  If (con==1.) Then
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck11: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
  Else
     dmass = amn + avpi
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck12: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
  End If
  fd5 = fd*p1*dmass
  Return
End Function fd5
