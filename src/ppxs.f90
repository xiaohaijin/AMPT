Subroutine ppxs(lb1, lb2, srt, ppsig, spprho, ipp)
  Parameter (amp=0.14, pi=3.1415926)
  Save
  ppsig = 0.0
  spprho = 0.0
  ipp = 0
  If (srt<=0.3) Return
  q = sqrt((srt/2)**2-amp**2)
  esigma = 5.8*amp
  tsigma = 2.06*q
  erho = 0.77
  trho = 0.095*q*(q/amp/(1.+(q/erho)**2))**2
  esi = esigma - srt
  If (esi==0) Then
     d00 = pi/2.
     Goto 10
  End If
  d00 = atan(tsigma/2./esi)
10 erh = erho - srt
  If (erh==0.) Then
     d11 = pi/2.
     Goto 20
  End If
  d11 = atan(trho/2./erh)
20 d20 = -0.12*q/amp
  s0 = 8.*pi*sin(d00)**2/q**2
  s1 = 8*pi*3*sin(d11)**2/q**2
  s2 = 8*pi*5*sin(d20)**2/q**2
  s0 = s0*0.197**2*10.
  s1 = s1*0.197**2*10.
  s2 = s2*0.197**2*10.
  spprho = s1/2.
  If (lb1==5 .And. lb2==5) Then
     ipp = 1
     ppsig = s2
     Return
  End If
  If ((lb1==5 .And. lb2==4) .Or. (lb1==4 .And. lb2==5)) Then
     ipp = 2
     ppsig = s2/2. + s1/2.
     Return
  End If
  If ((lb1==5 .And. lb2==3) .Or. (lb1==3 .And. lb2==5)) Then
     ipp = 3
     ppsig = s2/6. + s1/2. + s0/3.
     Return
  End If
  If (lb1==4 .And. lb2==4) Then
     ipp = 4
     ppsig = 2*s2/3. + s0/3.
     Return
  End If
  If ((lb1==4 .And. lb2==3) .Or. (lb1==3 .And. lb2==4)) Then
     ipp = 5
     ppsig = s2/2. + s1/2.
     Return
  End If
  If (lb1==3 .And. lb2==3) Then
     ipp = 6
     ppsig = s2
  End If
  Return
End Subroutine ppxs
