Subroutine crppba(px, py, pz, srt, i1, i2, iblock)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
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
  Call pbarfs(srt, npion, iseed)
  nchrg1 = 3 + int(3*ranart(nseed))
  nchrg2 = 3 + int(3*ranart(nseed))
  pmass1 = ap1
  pmass2 = ap1
  If (nchrg1==3 .Or. nchrg1==5) pmass1 = ap2
  If (nchrg2==3 .Or. nchrg2==5) pmass2 = ap2
  If (npion==2) Then
    iblock = 1902
    lb(i1) = nchrg1
    e(i1) = pmass1
    lb(i2) = nchrg2
    e(i2) = pmass2
    Goto 50
  End If
  If (npion==3) Then
    iblock = 1903
    lb(i1) = nchrg1
    e(i1) = pmass1
    lb(i2) = 22 + nchrg2
    e(i2) = amrho
    Goto 50
  End If
  If (npion==4) Then
    iblock = 1904
    If (ranart(nseed)>=0.5) Then
      lb(i1) = 22 + nchrg1
      e(i1) = amrho
      lb(i2) = 22 + nchrg2
      e(i2) = amrho
    Else
      lb(i1) = nchrg1
      e(i1) = pmass1
      lb(i2) = 28
      e(i2) = amomga
    End If
    Goto 50
  End If
  If (npion==5) Then
    iblock = 1905
    lb(i1) = 22 + nchrg1
    e(i1) = amrho
    lb(i2) = 28
    e(i2) = amomga
    Goto 50
  End If
  If (npion==6) Then
    iblock = 1906
    lb(i1) = 28
    e(i1) = amomga
    lb(i2) = 28
    e(i2) = amomga
  End If
  50 em1 = e(i1)
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
End Subroutine crppba
