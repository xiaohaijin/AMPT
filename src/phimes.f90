Subroutine phimes(i1, i2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, sigphi)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895, aomega=0.7819, arho=0.77, aphi=1.02)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Save
  s = srt**2
  sigphi = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  xsk7 = 0.0
  em1 = e(i1)
  em2 = e(i2)
  lb1 = lb(i1)
  lb2 = lb(i2)
  akap = aka
  xsk1 = 5.0
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<=0) Then
    Write (99, *) 'scheck48: ', scheck
    Stop
  End If
  pii = sqrt(scheck)
  If (lb1==23 .Or. lb2==23 .Or. lb1==21 .Or. lb2==21) Then
    If (srt>(ap1+akap)) Then
      pff = sqrt((s-(ap1+akap)**2)*(s-(ap1-akap)**2))
      xsk2 = 195.639*pff/pii/32./pi/s
    End If
    If (srt>(arho+akap)) Then
      pff = sqrt((s-(arho+akap)**2)*(s-(arho-akap)**2))
      xsk3 = 526.702*pff/pii/32./pi/s
    End If
    If (srt>(aomega+akap)) Then
      pff = sqrt((s-(aomega+akap)**2)*(s-(aomega-akap)**2))
      xsk4 = 355.429*pff/pii/32./pi/s
    End If
    If (srt>(ap1+aks)) Then
      pff = sqrt((s-(ap1+aks)**2)*(s-(ap1-aks)**2))
      xsk5 = 2047.042*pff/pii/32./pi/s
    End If
    If (srt>(arho+aks)) Then
      pff = sqrt((s-(arho+aks)**2)*(s-(arho-aks)**2))
      xsk6 = 1371.257*pff/pii/32./pi/s
    End If
    If (srt>(aomega+aks)) Then
      pff = sqrt((s-(aomega+aks)**2)*(s-(aomega-aks)**2))
      xsk7 = 482.292*pff/pii/32./pi/s
    End If
  Else If (iabs(lb1)==30 .Or. iabs(lb2)==30) Then
    If (srt>(ap1+akap)) Then
      pff = sqrt((s-(ap1+akap)**2)*(s-(ap1-akap)**2))
      xsk2 = 372.378*pff/pii/32./pi/s
    End If
    If (srt>(arho+akap)) Then
      pff = sqrt((s-(arho+akap)**2)*(s-(arho-akap)**2))
      xsk3 = 1313.960*pff/pii/32./pi/s
    End If
    If (srt>(aomega+akap)) Then
      pff = sqrt((s-(aomega+akap)**2)*(s-(aomega-akap)**2))
      xsk4 = 440.558*pff/pii/32./pi/s
    End If
    If (srt>(ap1+aks)) Then
      pff = sqrt((s-(ap1+aks)**2)*(s-(ap1-aks)**2))
      xsk5 = 1496.692*pff/pii/32./pi/s
    End If
    If (srt>(arho+aks)) Then
      pff = sqrt((s-(arho+aks)**2)*(s-(arho-aks)**2))
      xsk6 = 6999.840*pff/pii/32./pi/s
    End If
    If (srt>(aomega+aks)) Then
      pff = sqrt((s-(aomega+aks)**2)*(s-(aomega-aks)**2))
      xsk7 = 1698.903*pff/pii/32./pi/s
    End If
  Else
    srr1 = em1 + em2
    If (srt>(akap+akap)) Then
      srrt = srt - srr1
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk2 = 1.69/(srrt**0.141-0.407)
      Else
        xsk2 = 3.74 + 0.008*srrt**1.9
      End If
    End If
    If (srt>(akap+aks)) Then
      srr2 = akap + aks
      srr = amax1(srr1, srr2)
      srrt = srt - srr
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk3 = 1.69/(srrt**0.141-0.407)
      Else
        xsk3 = 3.74 + 0.008*srrt**1.9
      End If
    End If
    If (srt>(aks+aks)) Then
      srr2 = aks + aks
      srr = amax1(srr1, srr2)
      srrt = srt - srr
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk4 = 1.69/(srrt**0.141-0.407)
      Else
        xsk4 = 3.74 + 0.008*srrt**1.9
      End If
    End If
  End If
  sigphi = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6 + xsk7
  Return
End Subroutine phimes
