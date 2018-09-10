Subroutine ksreso(i1, i2)
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
  If (lb(i2)==21 .Or. lb(i2)==23) Then
    e(i1) = 0.
    i = i2
  Else
    e(i2) = 0.
    i = i1
  End If
  If (lb(i)==23) Then
    lb(i) = 30
  Else If (lb(i)==21) Then
    lb(i) = -30
  End If
  p(1, i) = p(1, i1) + p(1, i2)
  p(2, i) = p(2, i1) + p(2, i2)
  p(3, i) = p(3, i1) + p(3, i2)
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
    Write (99, *) 'scheck49: ', scheck
    Write (99, *) 'scheck49', scheck, e10, e20, p(1, i), p(2, i), p(3, i)
    Write (99, *) 'scheck49-1', e(i1), p(1, i1), p(2, i1), p(3, i1)
    Write (99, *) 'scheck49-2', e(i2), p(1, i2), p(2, i2), p(3, i2)
  End If
  dm = sqrt(sngl(scheck))
  e(i) = dm
  Return
End Subroutine ksreso
