Subroutine crlaba(px, py, pz, srt, brel, brsgm, i1, i2, nt, iblock, nchrg, icase)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (etam=0.5475, aomega=0.782, arho=0.77)
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
  If (icase==3) Then
     rrr = ranart(nseed)
     If (rrr<brel) Then
        iblock = 8
     Else
        iblock = 100
        If (rrr<(brel+brsgm)) Then
           lb(i1) = -15 - int(3*ranart(nseed))
           e(i1) = asa
        Else
           lb(i1) = -14
           e(i1) = ala
        End If
        lb(i2) = 3 + int(3*ranart(nseed))
        e(i2) = 0.138
     End If
  End If
  If (icase==4) Then
     rrr = ranart(nseed)
     If (rrr<brel) Then
        iblock = 8
     Else
        iblock = 102
        lb(i1) = 23
        lb(i2) = -1 - int(2*ranart(nseed))
        If (nchrg==-2) lb(i2) = -6
        If (nchrg==1) lb(i2) = -9
        e(i1) = aka
        e(i2) = 0.938
        If (nchrg==-2 .Or. nchrg==1) e(i2) = 1.232
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
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crlaba
