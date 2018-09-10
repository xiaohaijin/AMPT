Real Function fns(dmass, srt, con)
  Save
  width = 0.2
  amn = 0.938869
  avpi = 0.13803333
  an0 = 1.43
  fn = 4.*(an0**2)*width/((dmass**2-1.44**2)**2+an0**2*width**2)
  If (con==1.) Then
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck13: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
  Else
     dmass = amn + avpi
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck14: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
  End If
  fns = fn*p1*dmass
  Return
End Function fns
