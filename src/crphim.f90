Subroutine crphim(px, py, pz, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, sigphi, ikkg, ikkl, iblock)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, arho=0.77, aomega=0.7819, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895)
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
  lb1 = lb(i1)
  lb2 = lb(i2)
  x1 = ranart(nseed)*sigphi
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
  xsk6 = xsk5 + xsk6
  If (x1<=xsk1) Then
    iblock = 20
    Goto 100
  Else
    If (lb1==23 .Or. lb1==21 .Or. iabs(lb1)==30 .Or. lb2==23 .Or. lb2==21 .Or. iabs(lb2)==30) Then
      If (lb1==23 .Or. lb2==23) Then
        ikkl = 1
        iblock = 224
        iad1 = 23
        iad2 = 30
      Else If (lb1==30 .Or. lb2==30) Then
        ikkl = 0
        iblock = 226
        iad1 = 23
        iad2 = 30
      Else If (lb1==21 .Or. lb2==21) Then
        ikkl = 1
        iblock = 124
        iad1 = 21
        iad2 = -30
      Else
        ikkl = 0
        iblock = 126
        iad1 = 21
        iad2 = -30
      End If
      If (x1<=xsk2) Then
        lb(i1) = 3 + int(3*ranart(nseed))
        lb(i2) = iad1
        e(i1) = ap1
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk3) Then
        lb(i1) = 25 + int(3*ranart(nseed))
        lb(i2) = iad1
        e(i1) = arho
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk4) Then
        lb(i1) = 28
        lb(i2) = iad1
        e(i1) = aomega
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk5) Then
        lb(i1) = 3 + int(3*ranart(nseed))
        lb(i2) = iad2
        e(i1) = ap1
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      Else If (x1<=xsk6) Then
        lb(i1) = 25 + int(3*ranart(nseed))
        lb(i2) = iad2
        e(i1) = arho
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      Else
        lb(i1) = 28
        lb(i2) = iad2
        e(i1) = aomega
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      End If
    Else
      iblock = 223
      If (x1<=xsk2) Then
        lb(i1) = 23
        lb(i2) = 21
        e(i1) = aka
        e(i2) = aka
        ikkg = 2
        ikkl = 0
        Goto 100
      Else If (x1<=xsk3) Then
        lb(i1) = 23
        lb(i2) = -30
        If (ranart(nseed)<=0.5) Then
          lb(i1) = 21
          lb(i2) = 30
        End If
        e(i1) = aka
        e(i2) = aks
        ikkg = 1
        ikkl = 0
        Goto 100
      Else If (x1<=xsk4) Then
        lb(i1) = 30
        lb(i2) = -30
        e(i1) = aks
        e(i2) = aks
        ikkg = 0
        ikkl = 0
        Goto 100
      End If
    End If
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
End Subroutine crphim
