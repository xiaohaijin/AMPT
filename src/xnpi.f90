Real Function xnpi(i1, i2, la, xmax)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Save
  avmass = 0.5*(amn+amp)
  avpi = (2.*ap2+ap1)/3.
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck19: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  If (dm<=1.1) Then
     xnpi = 1.E-09
     Return
  End If
  If (la==1) Then
     gam = width(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.232)**2)
     pdelt2 = 0.051622
     Goto 10
  End If
  If (la==0) Then
     gam = w1440(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.440)**2)
     pdelt2 = 0.157897
     Goto 10
  End If
  If (la==2) Then
     gam = w1535(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.535)**2)
     pdelt2 = 0.2181
  End If
10 pstar2 = ((dm**2-avmass**2+avpi**2)/(2.*dm))**2 - avpi**2
  If (pstar2<=0.) Then
     xnpi = 1.E-09
  Else
     xnpi = f1*(pdelt2/pstar2)*xmax/10.
  End If
  Return
End Function xnpi
