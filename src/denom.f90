Real Function denom(srt, con)
  Parameter (ap1=0.13496, ap2=0.13957, pi=3.1415926, avmass=0.9383)
  Save
  avpi = (ap1+2.*ap2)/3.
  am0 = 1.232
  amn = avmass
  amp = avpi
  amax = srt - avmass
  amin = avmass + avpi
  nmax = 200
  dmass = (amax-amin)/float(nmax)
  sum = 0.
  Do i = 1, nmax + 1
     dm = amin + float(i-1)*dmass
     If (con==1.) Then
        q2 = ((dm**2-amn**2+amp**2)/(2.*dm))**2 - amp**2
        If (q2>0.) Then
           q = sqrt(q2)
        Else
           q = 1.E-06
        End If
        tq = 0.47*(q**3)/(amp**2*(1.+0.6*(q/amp)**2))
     Else If (con==2) Then
        tq = 0.2
        am0 = 1.44
     Else If (con==-1.) Then
        tq = 0.1
        am0 = 1.535
     End If
     a1 = 4.*tq*am0**2/(am0**2*tq**2+(dm**2-am0**2)**2)
     s = srt**2
     p0 = (s+dm**2-amn**2)**2/(4.*s) - dm**2
     If (p0<=0.) Then
        p1 = 1.E-06
     Else
        p1 = sqrt(p0)
     End If
     f = dm*a1*p1
     If ((i==1) .Or. (i==(nmax+1))) Then
        sum = sum + f*0.5
     Else
        sum = sum + f
     End If
  End Do
  denom = sum*dmass/(2.*pi)
  Return
End Function denom
