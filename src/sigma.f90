Real Function sigma(srt, id, ioi, iof)
  Parameter (amu=0.9383, amp=0.1384, pi=3.1415926, hc=0.19733)
  Save
  If (id==1) Then
     amass0 = 1.22
     t0 = 0.12
  Else
     amass0 = 1.43
     t0 = 0.2
  End If
  If ((ioi==1) .And. (iof==1)) Then
     alfa = 3.772
     beta = 1.262
     am0 = 1.188
     t = 0.09902
  End If
  If ((ioi==1) .And. (iof==0)) Then
     alfa = 15.28
     beta = 0.
     am0 = 1.245
     t = 0.1374
  End If
  If ((ioi==0) .And. (iof==1)) Then
     alfa = 146.3
     beta = 0.
     am0 = 1.472
     t = 0.02649
  End If
  zplus = (srt-amu-amass0)*2./t0
  zminus = (amu+amp-amass0)*2./t0
  deln = atan(zplus) - atan(zminus)
  If (deln==0) deln = 1.E-06
  amass = amass0 + (t0/4.)*alog((1.+zplus**2)/(1.+zminus**2))/deln
  s = srt**2
  p2 = s/4. - amu**2
  s0 = (amu+am0)**2
  p02 = s0/4. - amu**2
  p0 = sqrt(p02)
  pr2 = (s-(amu-amass)**2)*(s-(amu+amass)**2)/(4.*s)
  If (pr2>1.E-06) Then
     pr = sqrt(pr2)
  Else
     pr = 0.
     sigma = 1.E-06
     Return
  End If
  ss = amass**2
  q2 = (ss-(amu-amp)**2)*(ss-(amu+amp)**2)/(4.*ss)
  If (q2>1.E-06) Then
     q = sqrt(q2)
  Else
     q = 0.
     sigma = 1.E-06
     Return
  End If
  ss0 = am0**2
  q02 = (ss0-(amu-amp)**2)*(ss0-(amu+amp)**2)/(4.*ss0)
  scheck = q02
  If (scheck<0) Then
     Write (99, *) 'scheck20: ', scheck
     scheck = 0.
  End If
  q0 = sqrt(scheck)
  sigma = pi*(hc)**2/(2.*p2)*alfa*(pr/p0)**beta*am0**2*t**2*(q/q0)**3/((ss-am0**2)**2+am0**2*t**2)
  sigma = sigma*10.
  If (sigma==0) sigma = 1.E-06
  Return
End Function sigma
