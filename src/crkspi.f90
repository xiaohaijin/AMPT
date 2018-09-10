Subroutine crkspi(i1, i2, xsk1, xsk2, xsk3, xsk4, sigk, iblock, lbp1, lbp2, emm1, emm2)
  Parameter (maxstr=150001, maxr=1)
  Parameter (ap1=0.13496, ap2=0.13957, rhom=0.770, pi=3.1415926)
  Parameter (aeta=0.548, amomga=0.782)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  iblock = 466
  x1 = ranart(nseed)*sigk
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  If (x1<=xsk1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = rhom
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = ap2
    e(i2) = amomga
  Else If (x1<=xsk3) Then
    lb(i1) = 0
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = aeta
    e(i2) = rhom
  Else
    lb(i1) = 0
    lb(i2) = 28
    e(i1) = aeta
    e(i2) = amomga
  End If
  If (lb(i1)==4) e(i1) = ap1
  lbp1 = lb(i1)
  lbp2 = lb(i2)
  emm1 = e(i1)
  emm2 = e(i2)
  Return
End Subroutine crkspi
