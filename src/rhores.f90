Subroutine rhores(i1, i2)
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
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))
  p(1, i1) = p(1, i1) + p(1, i2)
  p(2, i1) = p(2, i1) + p(2, i2)
  p(3, i1) = p(3, i1) + p(3, i2)
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck18: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  e(i1) = dm
  e(i2) = 0
  Return
End Subroutine rhores
