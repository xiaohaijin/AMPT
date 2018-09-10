Subroutine crkhyp(px, py, pz, srt, i1, i2, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk, ikmp, iblock)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, aphi=1.02, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (pimass=0.140, ameta=0.5473, aka=0.498, aml=1.116, ams=1.193, am1440=1.44, am1535=1.535)
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
  iblock = 1908
  x1 = ranart(nseed)*sigk
  xky2 = xky1 + xky2
  xky3 = xky2 + xky3
  xky4 = xky3 + xky4
  xky5 = xky4 + xky5
  xky6 = xky5 + xky6
  xky7 = xky6 + xky7
  xky8 = xky7 + xky8
  xky9 = xky8 + xky9
  xky10 = xky9 + xky10
  xky11 = xky10 + xky11
  xky12 = xky11 + xky12
  xky13 = xky12 + xky13
  xky14 = xky13 + xky14
  xky15 = xky14 + xky15
  xky16 = xky15 + xky16
  If (x1<=xky1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = amp
    Goto 100
  Else If (x1<=xky2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = pimass
    e(i2) = am0
    Goto 100
  Else If (x1<=xky3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky4) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = amp
    Goto 100
  Else If (x1<=xky6) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = amrho
    e(i2) = am0
    Goto 100
  Else If (x1<=xky7) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky8) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky9) Then
    lb(i1) = 28
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = amp
    Goto 100
  Else If (x1<=xky10) Then
    lb(i1) = 28
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = amomga
    e(i2) = am0
    Goto 100
  Else If (x1<=xky11) Then
    lb(i1) = 28
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky12) Then
    lb(i1) = 28
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky13) Then
    lb(i1) = 0
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = amp
    Goto 100
  Else If (x1<=xky14) Then
    lb(i1) = 0
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = ameta
    e(i2) = am0
    Goto 100
  Else If (x1<=xky15) Then
    lb(i1) = 0
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky16) Then
    lb(i1) = 0
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = am1535
    Goto 100
  Else
    lb(i1) = 29
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = aphi
    e(i2) = amn
    iblock = 222
    Goto 100
  End If
  100 Continue
  If (ikmp==-1) lb(i2) = -lb(i2)
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
End Subroutine crkhyp
