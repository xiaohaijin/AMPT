Subroutine crpp(px, py, pz, srt, i1, i2, iblock, ppel, ppin, spprho, ipp)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /rndf77/nseed
  Save
  lb1i = lb(i1)
  lb2i = lb(i2)
  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1
  If (srt>(2*aka) .And. (ppin/(ppin+ppel))>ranart(nseed)) Then
     ranpi = ranart(nseed)
     If ((pprr/ppin)>=ranpi) Then
        Call pi2ro2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If ((pprr+ppee)/ppin>=ranpi) Then
        Call pi2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe)/ppin)>=ranpi) Then
        Call pi3eta(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre)/ppin)>=ranpi) Then
        Call rpiret(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre+xopoe)/ppin)>=ranpi) Then
        Call opioet(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre+xopoe+rree)/ppin)>=ranpi) Then
        Call ro2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre+xopoe+rree+ppinnb)/ppin)>=ranpi) Then
        Call bbarfs(lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else
        iblock = 66
        ei1 = aka
        ei2 = aka
        lbb1 = 21
        lbb2 = 23
        lb1 = lb(i1)
        lb2 = lb(i2)
        If (((lb1==0 .Or. (lb1>=3 .And. lb1<=5)) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==0 .Or. (lb2>=3 .And. lb2<=5)) .And. (lb1>=25 .And. lb1<=28))) Then
           ei1 = aks
           ei2 = aka
           If (ranart(nseed)>=0.5) Then
              iblock = 366
              lbb1 = 30
              lbb2 = 21
           Else
              iblock = 367
              lbb1 = -30
              lbb2 = 23
           End If
        End If
     End If
     e(i1) = ei1
     e(i2) = ei2
     lb(i1) = lbb1
     lb(i2) = lbb2
  Else
     If ((lb(i1)<3 .Or. lb(i1)>5) .And. (lb(i2)<3 .Or. lb(i2)>5)) Return
     iblock = 6
     If (ipp==1 .Or. ipp==4 .Or. ipp==6) Goto 10
     If (spprho/ppel>ranart(nseed)) Goto 20
  End If
10 ntag = 0
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
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
20 Continue
  iblock = 666
  Call rhores(i1, i2)
  If (ipp==2) lb(i1) = 27
  If (ipp==3) lb(i1) = 26
  If (ipp==5) lb(i1) = 25
  Return
End Subroutine crpp
