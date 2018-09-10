Subroutine wida1(dmass, rhomp, wa1, iseed)
  Save
  pimass = 0.137265
  coupa = 14.8
  rhomax = dmass - pimass - 0.02
  If (rhomax<=0) Then
     rhomp = 0.
     wa1 = -10.
  End If
  icount = 0
711 rhomp = rhomas(rhomax, iseed)
  icount = icount + 1
  If (dmass<=(pimass+rhomp)) Then
     If (icount<=100) Then
        Goto 711
     Else
        rhomp = 0.
        wa1 = -10.
        Return
     End If
  End If
  qqp2 = (dmass**2-(rhomp+pimass)**2)*(dmass**2-(rhomp-pimass)**2)
  qqp = sqrt(qqp2)/(2.0*dmass)
  epi = sqrt(pimass**2+qqp**2)
  erho = sqrt(rhomp**2+qqp**2)
  epirho = 2.0*(epi*erho+qqp**2)**2 + rhomp**2*epi**2
  wa1 = coupa**2*qqp*epirho/(24.0*3.1416*dmass**2)
  Return
End Subroutine wida1
