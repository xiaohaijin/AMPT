Subroutine cren(px, py, pz, srt, i1, i2, iblock)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
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
  ntag = 0
  iblock = 7
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) Then
     ianti = 1
     iblock = -7
  End If
  kaonc = 0
  If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
  If (e(i1)<=0.6) Then
     lb(i1) = 23
     e(i1) = aka
     If (kaonc==1) Then
        lb(i2) = 14
        e(i2) = ala
     Else
        lb(i2) = 15 + int(3*ranart(nseed))
        e(i2) = asa
     End If
     If (ianti==1) Then
        lb(i1) = 21
        lb(i2) = -lb(i2)
     End If
  Else
     lb(i2) = 23
     e(i2) = aka
     If (kaonc==1) Then
        lb(i1) = 14
        e(i1) = ala
     Else
        lb(i1) = 15 + int(3*ranart(nseed))
        e(i1) = asa
     End If
     If (ianti==1) Then
        lb(i2) = 21
        lb(i1) = -lb(i1)
     End If
  End If
  em1 = e(i1)
  em2 = e(i2)
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Return
End Subroutine cren
