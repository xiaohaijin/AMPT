Subroutine crkkpi(i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigk, iblock, lbp1, lbp2, emm1, emm2)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ameta=0.5473, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  iblock = 1907
  x1 = ranart(nseed)*sigk
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
  xsk6 = xsk5 + xsk6
  xsk7 = xsk6 + xsk7
  xsk8 = xsk7 + xsk8
  xsk9 = xsk8 + xsk9
  xsk10 = xsk9 + xsk10
  If (x1<=xsk1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 3 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = ap2
    Goto 100
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = amrho
    Goto 100
  Else If (x1<=xsk3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = ap2
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk4) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 0
    e(i1) = ap2
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = amrho
    e(i2) = amrho
    Goto 100
  Else If (x1<=xsk6) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = amrho
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk7) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 0
    e(i1) = amrho
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk8) Then
    lb(i1) = 28
    lb(i2) = 28
    e(i1) = amomga
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk9) Then
    lb(i1) = 28
    lb(i2) = 0
    e(i1) = amomga
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk10) Then
    lb(i1) = 0
    lb(i2) = 0
    e(i1) = ameta
    e(i2) = ameta
  Else
    iblock = 222
    Call rhores(i1, i2)
    lb(i1) = 29
    e(i2) = 0.
  End If
  100 Continue
  lbp1 = lb(i1)
  lbp2 = lb(i2)
  emm1 = e(i1)
  emm2 = e(i2)
  Return
End Subroutine crkkpi
