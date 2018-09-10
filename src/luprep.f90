Subroutine luprep(ip)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension dps(5), dpc(5), ue(3)
  i1 = n
  Do mqgst = 1, 2
    Do i = max(1, ip), n
      If (k(i,1)/=3) Goto 120
      kc = lucomp(k(i,2))
      If (kc==0) Goto 120
      kq = kchg(kc, 2)
      If (kq==0 .Or. (mqgst==1 .And. kq==2)) Goto 120
      kcs = 4
      If (kq*isign(1,k(i,2))<0) kcs = 5
      ia = i
      nstp = 0
      100 nstp = nstp + 1
      If (nstp>4*n) Then
        Call luerrm(14, '(LUPREP:) caught in infinite loop')
        Return
      End If
      If (k(ia,1)==3) Then
        If (i1>=mstu(4)-mstu(32)-5) Then
          Call luerrm(11, '(LUPREP:) no more memory left in LUJETS')
          Return
        End If
        i1 = i1 + 1
        k(i1, 1) = 2
        If (nstp>=2 .And. iabs(k(ia,2))/=21) k(i1, 1) = 1
        k(i1, 2) = k(ia, 2)
        k(i1, 3) = ia
        k(i1, 4) = 0
        k(i1, 5) = 0
        Do j = 1, 5
          p(i1, j) = p(ia, j)
          v(i1, j) = v(ia, j)
        End Do
        k(ia, 1) = k(ia, 1) + 10
        If (k(i1,1)==1) Goto 120
      End If
      ib = ia
      If (mod(k(ib,kcs)/mstu(5)**2,2)==0 .And. mod(k(ib,kcs),mstu(5))/=0) Then
        ia = mod(k(ib,kcs), mstu(5))
        k(ib, kcs) = k(ib, kcs) + mstu(5)**2
        mrev = 0
      Else
        If (k(ib,kcs)>=2*mstu(5)**2 .Or. mod(k(ib,kcs)/mstu(5),mstu(5))==0) kcs = 9 - kcs
        ia = mod(k(ib,kcs)/mstu(5), mstu(5))
        k(ib, kcs) = k(ib, kcs) + 2*mstu(5)**2
        mrev = 1
      End If
      If (ia<=0 .Or. ia>n) Then
        Call luerrm(12, '(LUPREP:) colour rearrangement failed')
        Return
      End If
      If (mod(k(ia,4)/mstu(5),mstu(5))==ib .Or. mod(k(ia,5)/mstu(5),mstu(5))==ib) Then
        If (mrev==1) kcs = 9 - kcs
        If (mod(k(ia,kcs)/mstu(5),mstu(5))/=ib) kcs = 9 - kcs
        k(ia, kcs) = k(ia, kcs) + 2*mstu(5)**2
      Else
        If (mrev==0) kcs = 9 - kcs
        If (mod(k(ia,kcs),mstu(5))/=ib) kcs = 9 - kcs
        k(ia, kcs) = k(ia, kcs) + mstu(5)**2
      End If
      If (ia/=i) Goto 100
      k(i1, 1) = 1
    120 End Do
  End Do
  n = i1
  If (mstj(14)<=0) Goto 320
  ns = n
  140 nsin = n - ns
  pdm = 1. + parj(32)
  ic = 0
  Do i = max(1, ip), ns
    If (k(i,1)/=1 .And. k(i,1)/=2) Then
    Else If (k(i,1)==2 .And. ic==0) Then
      nsin = nsin + 1
      ic = i
      Do j = 1, 4
        dps(j) = dble(p(i,j))
      End Do
      mstj(93) = 1
      dps(5) = dble(ulmass(k(i,2)))
    Else If (k(i,1)==2) Then
      Do j = 1, 4
        dps(j) = dps(j) + dble(p(i,j))
      End Do
    Else If (ic/=0 .And. kchg(lucomp(k(i,2)),2)/=0) Then
      Do j = 1, 4
        dps(j) = dps(j) + dble(p(i,j))
      End Do
      mstj(93) = 1
      dps(5) = dps(5) + dble(ulmass(k(i,2)))
      pd = sngl(sqrt(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2))-dps(5))
      If (pd<pdm) Then
        pdm = pd
        Do j = 1, 5
          dpc(j) = dps(j)
        End Do
        ic1 = ic
        ic2 = i
      End If
      ic = 0
    Else
      nsin = nsin + 1
    End If
  End Do
  If (pdm>=parj(32)) Goto 320
  nsav = n
  pecm = sngl(sqrt(max(0D0,dpc(4)**2-dpc(1)**2-dpc(2)**2-dpc(3)**2)))
  k(n+1, 1) = 11
  k(n+1, 2) = 91
  k(n+1, 3) = ic1
  k(n+1, 4) = n + 2
  k(n+1, 5) = n + 3
  p(n+1, 1) = sngl(dpc(1))
  p(n+1, 2) = sngl(dpc(2))
  p(n+1, 3) = sngl(dpc(3))
  p(n+1, 4) = sngl(dpc(4))
  p(n+1, 5) = pecm
  k(n+2, 1) = 1
  k(n+3, 1) = 1
  If (mstu(16)/=2) Then
    k(n+2, 3) = n + 1
    k(n+3, 3) = n + 1
  Else
    k(n+2, 3) = ic1
    k(n+3, 3) = ic2
  End If
  k(n+2, 4) = 0
  k(n+3, 4) = 0
  k(n+2, 5) = 0
  k(n+3, 5) = 0
  If (iabs(k(ic1,2))/=21) Then
    kc1 = lucomp(k(ic1,2))
    kc2 = lucomp(k(ic2,2))
    If (kc1==0 .Or. kc2==0) Goto 320
    kq1 = kchg(kc1, 2)*isign(1, k(ic1,2))
    kq2 = kchg(kc2, 2)*isign(1, k(ic2,2))
    If (kq1+kq2/=0) Goto 320
    200 Call lukfdi(k(ic1,2), 0, kfln, k(n+2,2))
    Call lukfdi(k(ic2,2), -kfln, kfldmp, k(n+3,2))
    If (k(n+2,2)==0 .Or. k(n+3,2)==0) Goto 200
  Else
    If (iabs(k(ic2,2))/=21) Goto 320
    210 Call lukfdi(1+int((2.+parj(2))*rlu(0)), 0, kfln, kfdmp)
    Call lukfdi(kfln, 0, kflm, k(n+2,2))
    Call lukfdi(-kfln, -kflm, kfldmp, k(n+3,2))
    If (k(n+2,2)==0 .Or. k(n+3,2)==0) Goto 210
  End If
  p(n+2, 5) = ulmass(k(n+2,2))
  p(n+3, 5) = ulmass(k(n+3,2))
  If (p(n+2,5)+p(n+3,5)+parj(64)>=pecm .And. nsin==1) Goto 320
  If (p(n+2,5)+p(n+3,5)+parj(64)>=pecm) Goto 260
  If (dble(pecm)>=0.02D0*dpc(4)) Then
    pa = sqrt((pecm**2-(p(n+2,5)+p(n+3,5))**2)*(pecm**2-(p(n+2,5)-p(n+3,5))**2))/(2.*pecm)
    ue(3) = 2.*rlu(0) - 1.
    phi = paru(2)*rlu(0)
    ue(1) = sqrt(1.-ue(3)**2)*cos(phi)
    ue(2) = sqrt(1.-ue(3)**2)*sin(phi)
    Do j = 1, 3
      p(n+2, j) = pa*ue(j)
      p(n+3, j) = -pa*ue(j)
    End Do
    p(n+2, 4) = sqrt(pa**2+p(n+2,5)**2)
    p(n+3, 4) = sqrt(pa**2+p(n+3,5)**2)
    Call ludbrb(n+2, n+3, 0., 0., dpc(1)/dpc(4), dpc(2)/dpc(4), dpc(3)/dpc(4))
  Else
    np = 0
    Do i = ic1, ic2
      If (k(i,1)==1 .Or. k(i,1)==2) np = np + 1
    End Do
    ha = p(ic1, 4)*p(ic2, 4) - p(ic1, 1)*p(ic2, 1) - p(ic1, 2)*p(ic2, 2) - p(ic1, 3)*p(ic2, 3)
    If (np>=3 .Or. ha<=1.25*p(ic1,5)*p(ic2,5)) Goto 260
    hd1 = 0.5*(p(n+2,5)**2-p(ic1,5)**2)
    hd2 = 0.5*(p(n+3,5)**2-p(ic2,5)**2)
    hr = sqrt(max(0.,((ha-hd1-hd2)**2-(p(n+2,5)*p(n+3,5))**2)/(ha**2-(p(ic1,5)*p(ic2,5))**2))) - 1.
    hc = p(ic1, 5)**2 + 2.*ha + p(ic2, 5)**2
    hk1 = ((p(ic2,5)**2+ha)*hr+hd1-hd2)/hc
    hk2 = ((p(ic1,5)**2+ha)*hr+hd2-hd1)/hc
    Do j = 1, 4
      p(n+2, j) = (1.+hk1)*p(ic1, j) - hk2*p(ic2, j)
      p(n+3, j) = (1.+hk2)*p(ic2, j) - hk1*p(ic1, j)
    End Do
  End If
  Do j = 1, 4
    v(n+1, j) = v(ic1, j)
    v(n+2, j) = v(ic1, j)
    v(n+3, j) = v(ic2, j)
  End Do
  v(n+1, 5) = 0.
  v(n+2, 5) = 0.
  v(n+3, 5) = 0.
  n = n + 3
  Goto 300
  260 k(n+1, 5) = n + 2
  If (iabs(k(ic1,2))>100 .And. iabs(k(ic2,2))>100) Then
    Goto 320
  Else If (iabs(k(ic1,2))/=21) Then
    Call lukfdi(k(ic1,2), k(ic2,2), kfldmp, k(n+2,2))
  Else
    kfln = 1 + int((2.+parj(2))*rlu(0))
    Call lukfdi(kfln, -kfln, kfldmp, k(n+2,2))
  End If
  If (k(n+2,2)==0) Goto 260
  p(n+2, 5) = ulmass(k(n+2,2))
  ir = 0
  ha = 0.
  Do mcomb = 1, 3
    If (ir/=0) Goto 280
    Do i = max(1, ip), n
      If (k(i,1)<=0 .Or. k(i,1)>10 .Or. (i>=ic1 .And. i<=ic2 .And. k(i,1)>=1 .And. k(i,1)<=2)) Goto 270
      If (mcomb==1) kci = lucomp(k(i,2))
      If (mcomb==1 .And. kci==0) Goto 270
      If (mcomb==1 .And. kchg(kci,2)==0 .And. i<=ns) Goto 270
      If (mcomb==2 .And. iabs(k(i,2))>10 .And. iabs(k(i,2))<=100) Goto 270
      hcr = sngl(dpc(4))*p(i, 4) - sngl(dpc(1))*p(i, 1) - sngl(dpc(2))*p(i, 2) - sngl(dpc(3))*p(i, 3)
      If (hcr>ha) Then
        ir = i
        ha = hcr
      End If
    270 End Do
  280 End Do
  hb = pecm**2 + ha
  hc = p(n+2, 5)**2 + ha
  hd = p(ir, 5)**2 + ha
  hk2 = 0.0
  If (ha**2-(pecm*p(ir,5))**2==0.0 .Or. hb+hd==0.0) Goto 285
  hk2 = 0.5*(hb*sqrt(((hb+hc)**2-4.*(hb+hd)*p(n+2,5)**2)/(ha**2-(pecm*p(ir,5))**2))-(hb+hc))/(hb+hd)
  285 hk1 = (0.5*(p(n+2,5)**2-pecm**2)+hd*hk2)/hb
  Do j = 1, 4
    p(n+2, j) = (1.+hk1)*sngl(dpc(j)) - hk2*p(ir, j)
    p(ir, j) = (1.+hk2)*p(ir, j) - hk1*sngl(dpc(j))
    v(n+1, j) = v(ic1, j)
    v(n+2, j) = v(ic1, j)
  End Do
  v(n+1, 5) = 0.
  v(n+2, 5) = 0.
  n = n + 2
  300 Do i = ic1, ic2
    If ((k(i,1)==1 .Or. k(i,1)==2) .And. kchg(lucomp(k(i,2)),2)/=0) Then
      k(i, 1) = k(i, 1) + 10
      If (mstu(16)/=2) Then
        k(i, 4) = nsav + 1
        k(i, 5) = nsav + 1
      Else
        k(i, 4) = nsav + 2
        k(i, 5) = n
      End If
    End If
  End Do
  If (n<mstu(4)-mstu(32)-5) Goto 140
  320 np = 0
  kfn = 0
  kqs = 0
  Do j = 1, 5
    dps(j) = 0D0
  End Do
  Do i = max(1, ip), n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 360
    kc = lucomp(k(i,2))
    If (kc==0) Goto 360
    kq = kchg(kc, 2)*isign(1, k(i,2))
    If (kq==0) Goto 360
    np = np + 1
    If (kq/=2) Then
      kfn = kfn + 1
      kqs = kqs + kq
      mstj(93) = 1
      dps(5) = dps(5) + dble(ulmass(k(i,2)))
    End If
    Do j = 1, 4
      dps(j) = dps(j) + dble(p(i,j))
    End Do
    If (k(i,1)==1) Then
      If (np/=1 .And. (kfn==1 .Or. kfn>=3 .Or. kqs/=0)) Call luerrm(2, '(LUPREP:) unphysical flavour combination')
      If (np/=2 .And. dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2<(0.9D0*dble(parj(32))+dps(5))**2) Then
        Call luerrm(3, '(LUPREP:) too small mass in jet system')
        Write (6, *) 'DPS(1-5),KI1-5=', dps(1), dps(2), dps(3), dps(4), dps(5), '*', k(i, 1), k(i, 2), k(i, 3), k(i, 4), k(i, 5)
      End If
      np = 0
      kfn = 0
      kqs = 0
      Do j = 1, 5
        dps(j) = 0D0
      End Do
    End If
  360 End Do
  Return
End Subroutine luprep
