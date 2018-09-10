Real Function ang(srt, iseed)
  Common /rndf77/nseed
  Save
  If ((srt>2.14) .And. (srt<=2.4)) Then
     b1s = 29.03 - 23.75*srt + 4.865*srt**2
     b2s = -30.33 + 25.53*srt - 5.301*srt**2
  End If
  If (srt>2.4) Then
     b1s = 0.06
     b2s = 0.4
  End If
  x = ranart(nseed)
  p = b1s/b2s
  q = (2.*x-1.)*(b1s+b2s)/b2s
  If ((-q/2.+sqrt((q/2.)**2+(p/3.)**3))>=0.) Then
     ang1 = (-q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  Else
     ang1 = -(q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  End If
  If ((-q/2.-sqrt((q/2.)**2+(p/3.)**3)>=0.)) Then
     ang2 = (-q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  Else
     ang2 = -(q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  End If
  ang = ang1 + ang2
  Return
End Function ang
