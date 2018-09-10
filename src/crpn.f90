Subroutine crpn(px, py, pz, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
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
  iblock = 1
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
  If (xkaon0/(xkaon+xphi)>=x1) Then
     iblock = 7
     If (ianti==1) iblock = -7
     ntag = 0
     kaonc = 0
     If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
     If (e(i1)<=0.2) Then
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
     Goto 50
  Else If (xphi/(xkaon+xphi)>=x1) Then
     iblock = 222
     If (xphin/xphi>=ranart(nseed)) Then
        lb(i1) = 1 + int(2*ranart(nseed))
        e(i1) = amn
     Else
        lb(i1) = 6 + int(4*ranart(nseed))
        e(i1) = am0
     End If
     If (ianti==1) lb(i1) = -lb(i1)
     lb(i2) = 29
     e(i2) = aphi
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
  Else
     If (ranart(nseed)<=twopi(srt)/(twopi(srt)+threpi(srt)+fourpi(srt))) Then
        iblock = 77
     Else
        If (threpi(srt)/(threpi(srt)+fourpi(srt))>ranart(nseed)) Then
           iblock = 78
        Else
           iblock = 79
        End If
     End If
     ntag = 0
     x2 = ranart(nseed)
     If (iblock==77) Then
        dmax = srt - ap1 - 0.02
        dm = rmass(dmax, iseed)
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              Else
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 4
                 ipi = 4
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              Else
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              Else
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              Else
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
     End If
     If (iblock==78) Then
        Call rmasdd(srt, 1.232, 0.77, 1.08, 0.28, iseed, 4, dm, ameson)
        arho = ameson
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              Else
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              Else
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              Else
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              Else
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
              End If
           End If
        End If
     End If
     If (iblock==79) Then
        aomega = 0.782
        dmax = srt - 0.782 - 0.02
        dm = rmass(dmax, iseed)
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 9
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 9
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 7
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 7
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 8
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 8
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 8
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 8
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 6
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 6
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
           End If
        End If
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 7
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 7
              e(i2) = dm
              lb(i1) = 26
              e(i1) = arho
              Goto 40
           End If
        End If
     End If
40   em1 = e(i1)
     em2 = e(i2)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(ii) = -lb(ii)
        jj = i2
        If (ii==i2) jj = i1
        If (iblock==77) Then
           If (lb(jj)==3) Then
              lb(jj) = 5
           Else If (lb(jj)==5) Then
              lb(jj) = 3
           End If
        Else If (iblock==78) Then
           If (lb(jj)==25) Then
              lb(jj) = 27
           Else If (lb(jj)==27) Then
              lb(jj) = 25
           End If
        End If
     End If
  End If
50 pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 0.00000001
  pr = sqrt(pr2)/(2.*srt)
  xptr = 0.33*pr
  cc1 = ptr(xptr, iseed)
  scheck = pr**2 - cc1**2
  If (scheck<0) Then
     Write (99, *) 'scheck36: ', scheck
     scheck = 0.
  End If
  c1 = sqrt(scheck)/pr
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crpn
