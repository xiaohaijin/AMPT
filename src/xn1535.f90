Real Function xn1535(i1, i2, la)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, etam=0.5475, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
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
     Write (99, *) 'scheck21: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  If (dm<=1.1) Then
     xn1535 = 1.E-06
     Return
  End If
  gam = w1535(dm)
  gam0 = 0.15
  f1 = 0.25*gam0**2/(0.25*gam**2+(dm-1.535)**2)
  If (la==1) Then
     xmax = 11.3
  Else
     xmax = 74.
  End If
  xn1535 = f1*xmax/10.
  Return
End Function xn1535
