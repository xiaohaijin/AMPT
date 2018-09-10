Subroutine crphib(px, py, pz, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, sigp, iblock)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, arho=0.77)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  px0 = px
  py0 = py
  pz0 = pz
  iblock = 223
  x1 = ranart(nseed)*sigp
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
  If (x1<=xsk1) Then
    iblock = 20
    Goto 100
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = ap1
    e(i2) = amn
    Goto 100
  Else If (x1<=xsk3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = ap1
    e(i2) = am0
    Goto 100
  Else If (x1<=xsk4) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = arho
    e(i2) = amn
    Goto 100
  Else If (x1<=xsk5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = arho
    e(i2) = am0
    Goto 100
  Else
    lb(i1) = 23
    lb(i2) = 14
    e(i1) = aka
    e(i2) = ala
    iblock = 221
  End If
  100 Continue
  em1 = e(i1)
  em2 = e(i2)
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-08
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crphib
