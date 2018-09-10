Subroutine luindf(ip)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension dps(5), psi(4), nfi(3), nfl(3), ifet(3), kflf(3), kflo(2), pxo(2), pyo(2), wo(2)
  nsav = n
  njet = 0
  kqsum = 0
  Do j = 1, 5
    dps(j) = 0.D0
  End Do
  i = ip - 1
  110 i = i + 1
  If (i>min(n,mstu(4)-mstu(32))) Then
    Call luerrm(12, '(LUINDF:) failed to reconstruct jet system')
    If (mstu(21)>=1) Return
  End If
  If (k(i,1)/=1 .And. k(i,1)/=2) Goto 110
  kc = lucomp(k(i,2))
  If (kc==0) Goto 110
  kq = kchg(kc, 2)*isign(1, k(i,2))
  If (kq==0) Goto 110
  njet = njet + 1
  If (kq/=2) kqsum = kqsum + kq
  Do j = 1, 5
    k(nsav+njet, j) = k(i, j)
    p(nsav+njet, j) = p(i, j)
    dps(j) = dps(j) + dble(p(i,j))
  End Do
  k(nsav+njet, 3) = i
  If (k(i,1)==2 .Or. (mstj(3)<=5 .And. n>i .And. k(i+1,1)==2)) Goto 110
  If (njet/=1 .And. kqsum/=0) Then
    Call luerrm(12, '(LUINDF:) unphysical flavour combination')
    If (mstu(21)>=1) Return
  End If
  If (njet/=1) Call ludbrb(nsav+1, nsav+njet, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))
  pecm = 0.
  Do j = 1, 3
    nfi(j) = 0
  End Do
  Do i = nsav + 1, nsav + njet
    pecm = pecm + p(i, 4)
    kfa = iabs(k(i,2))
    If (kfa<=3) Then
      nfi(kfa) = nfi(kfa) + isign(1, k(i,2))
    Else If (kfa>1000) Then
      kfla = mod(kfa/1000, 10)
      kflb = mod(kfa/100, 10)
      If (kfla<=3) nfi(kfla) = nfi(kfla) + isign(1, k(i,2))
      If (kflb<=3) nfi(kflb) = nfi(kflb) + isign(1, k(i,2))
    End If
  End Do
  ntry = 0
  150 ntry = ntry + 1
  n = nsav + njet
  If (ntry>200) Then
    Call luerrm(14, '(LUINDF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  Do j = 1, 3
    nfl(j) = nfi(j)
    ifet(j) = 0
    kflf(j) = 0
  End Do
  Do ip1 = nsav + 1, nsav + njet
    mstj(91) = 0
    nsav1 = n
    kflh = iabs(k(ip1,2))
    If (kflh>10) kflh = mod(kflh/1000, 10)
    kflo(2) = 0
    wf = p(ip1, 4) + sqrt(p(ip1,1)**2+p(ip1,2)**2+p(ip1,3)**2)
    170 If (iabs(k(ip1,2))/=21) Then
      nstr = 1
      kflo(1) = k(ip1, 2)
      Call luptdi(0, pxo(1), pyo(1))
      wo(1) = wf
    Else If (mstj(2)<=2) Then
      nstr = 1
      If (mstj(2)==2) mstj(91) = 1
      kflo(1) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
      Call luptdi(0, pxo(1), pyo(1))
      wo(1) = wf
    Else
      nstr = 2
      If (mstj(2)==4) mstj(91) = 1
      kflo(1) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
      kflo(2) = -kflo(1)
      Call luptdi(0, pxo(1), pyo(1))
      pxo(2) = -pxo(1)
      pyo(2) = -pyo(1)
      wo(1) = wf*rlu(0)**(1./3.)
      wo(2) = wf - wo(1)
    End If
    Do istr = 1, nstr
      180 i = n
      irank = 0
      kfl1 = kflo(istr)
      px1 = pxo(istr)
      py1 = pyo(istr)
      w = wo(istr)
      190 i = i + 1
      If (i>=mstu(4)-mstu(32)-njet-5) Then
        Call luerrm(11, '(LUINDF:) no more memory left in LUJETS')
        If (mstu(21)>=1) Return
      End If
      irank = irank + 1
      k(i, 1) = 1
      k(i, 3) = ip1
      k(i, 4) = 0
      k(i, 5) = 0
      200 Call lukfdi(kfl1, 0, kfl2, k(i,2))
      If (k(i,2)==0) Goto 180
      If (mstj(12)>=3 .And. irank==1 .And. iabs(kfl1)<=10 .And. iabs(kfl2)>10) Then
        If (rlu(0)>parj(19)) Goto 200
      End If
      p(i, 5) = ulmass(k(i,2))
      Call luptdi(kfl1, px2, py2)
      p(i, 1) = px1 + px2
      p(i, 2) = py1 + py2
      pr = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
      Call luzdis(kfl1, kfl2, pr, z)
      p(i, 3) = 0.5*(z*w-pr/(z*w))
      p(i, 4) = 0.5*(z*w+pr/(z*w))
      If (mstj(3)>=1 .And. irank==1 .And. kflh>=4 .And. p(i,3)<=0.001) Then
        If (w>=p(i,5)+0.5*parj(32)) Goto 180
        p(i, 3) = 0.0001
        p(i, 4) = sqrt(pr)
        z = p(i, 4)/w
      End If
      kfl1 = -kfl2
      px1 = -px2
      py1 = -py2
      w = (1.-z)*w
      Do j = 1, 5
        v(i, j) = 0.
      End Do
      If (mstj(3)>=0 .And. p(i,3)<0.) i = i - 1
      If (w>parj(31)) Goto 190
      n = i
    End Do
    If (mod(mstj(3),5)==4 .And. n==nsav1) wf = wf + 0.1*parj(32)
    If (mod(mstj(3),5)==4 .And. n==nsav1) Goto 170
    the = ulangl(p(ip1,3), sqrt(p(ip1,1)**2+p(ip1,2)**2))
    phi = ulangl(p(ip1,1), p(ip1,2))
    Call ludbrb(nsav1+1, n, the, phi, 0D0, 0D0, 0D0)
    k(k(ip1,3), 4) = nsav1 + 1
    k(k(ip1,3), 5) = n
  End Do
  If (njet==1 .Or. mstj(3)<=0) Goto 470
  If (mod(mstj(3),5)/=0 .And. n-nsav-njet<2) Goto 150
  Do i = nsav + njet + 1, n
    kfa = iabs(k(i,2))
    kfla = mod(kfa/1000, 10)
    kflb = mod(kfa/100, 10)
    kflc = mod(kfa/10, 10)
    If (kfla==0) Then
      If (kflb<=3) nfl(kflb) = nfl(kflb) - isign(1, k(i,2))*(-1)**kflb
      If (kflc<=3) nfl(kflc) = nfl(kflc) + isign(1, k(i,2))*(-1)**kflb
    Else
      If (kfla<=3) nfl(kfla) = nfl(kfla) - isign(1, k(i,2))
      If (kflb<=3) nfl(kflb) = nfl(kflb) - isign(1, k(i,2))
      If (kflc<=3) nfl(kflc) = nfl(kflc) - isign(1, k(i,2))
    End If
  End Do
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nreq==0) Goto 320
  nrem = 0
  250 irem = 0
  p2min = pecm**2
  Do i = nsav + njet + 1, n
    p2 = p(i, 1)**2 + p(i, 2)**2 + p(i, 3)**2
    If (k(i,1)==1 .And. p2<p2min) irem = i
    If (k(i,1)==1 .And. p2<p2min) p2min = p2
  End Do
  If (irem==0) Goto 150
  k(irem, 1) = 7
  kfa = iabs(k(irem,2))
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  If (kfla>=4 .Or. kflb>=4) k(irem, 1) = 8
  If (k(irem,1)==8) Goto 250
  If (kfla==0) Then
    isgn = isign(1, k(irem,2))*(-1)**kflb
    If (kflb<=3) nfl(kflb) = nfl(kflb) + isgn
    If (kflc<=3) nfl(kflc) = nfl(kflc) - isgn
  Else
    If (kfla<=3) nfl(kfla) = nfl(kfla) + isign(1, k(irem,2))
    If (kflb<=3) nfl(kflb) = nfl(kflb) + isign(1, k(irem,2))
    If (kflc<=3) nfl(kflc) = nfl(kflc) + isign(1, k(irem,2))
  End If
  nrem = nrem + 1
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nreq>nrem) Goto 250
  Do i = nsav + njet + 1, n
    If (k(i,1)==8) k(i, 1) = 1
  End Do
  280 nfet = 2
  If (nfl(1)+nfl(2)+nfl(3)/=0) nfet = 3
  If (nreq<nrem) nfet = 1
  If (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))==0) nfet = 0
  Do j = 1, nfet
    ifet(j) = 1 + int((iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)))*rlu(0))
    kflf(j) = isign(1, nfl(1))
    If (ifet(j)>iabs(nfl(1))) kflf(j) = isign(2, nfl(2))
    If (ifet(j)>iabs(nfl(1))+iabs(nfl(2))) kflf(j) = isign(3, nfl(3))
  End Do
  If (nfet==2 .And. (ifet(1)==ifet(2) .Or. kflf(1)*kflf(2)>0)) Goto 280
  If (nfet==3 .And. (ifet(1)==ifet(2) .Or. ifet(1)==ifet(3) .Or. ifet(2)==ifet(3) .Or. kflf(1)*kflf(2)<0 .Or. kflf(1)*kflf(3)<0 .Or. kflf(1)*(nfl(1)+nfl(2)+nfl(3))<0)) Goto 280
  If (nfet==0) kflf(1) = 1 + int((2.+parj(2))*rlu(0))
  If (nfet==0) kflf(2) = -kflf(1)
  If (nfet==1) kflf(2) = isign(1+int((2.+parj(2))*rlu(0)), -kflf(1))
  If (nfet<=2) kflf(3) = 0
  If (kflf(3)/=0) Then
    kflfc = isign(1000*max(iabs(kflf(1)),iabs(kflf(3)))+100*min(iabs(kflf(1)),iabs(kflf(3)))+1, kflf(1))
    If (kflf(1)==kflf(3) .Or. (1.+3.*parj(4))*rlu(0)>1.) kflfc = kflfc + isign(2, kflfc)
  Else
    kflfc = kflf(1)
  End If
  Call lukfdi(kflfc, kflf(2), kfldmp, kf)
  If (kf==0) Goto 280
  Do j = 1, max(2, nfet)
    nfl(iabs(kflf(j))) = nfl(iabs(kflf(j))) - isign(1, kflf(j))
  End Do
  npos = min(1+int(rlu(0)*nrem), nrem)
  Do i = nsav + njet + 1, n
    If (k(i,1)==7) npos = npos - 1
    If (k(i,1)==1 .Or. npos/=0) Goto 310
    k(i, 1) = 1
    k(i, 2) = kf
    p(i, 5) = ulmass(k(i,2))
    p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  310 End Do
  nrem = nrem - 1
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nrem>0) Goto 280
  320 If (mod(mstj(3),5)/=0 .And. mod(mstj(3),5)/=4) Then
    Do j = 1, 3
      psi(j) = 0.
      Do i = nsav + njet + 1, n
        psi(j) = psi(j) + p(i, j)
      End Do
    End Do
    psi(4) = psi(1)**2 + psi(2)**2 + psi(3)**2
    pws = 0.
    Do i = nsav + njet + 1, n
      If (mod(mstj(3),5)==1) pws = pws + p(i, 4)
      If (mod(mstj(3),5)==2) pws = pws + sqrt(p(i,5)**2+(psi(1)*p(i,1)+psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
      If (mod(mstj(3),5)==3) pws = pws + 1.
    End Do
    Do i = nsav + njet + 1, n
      If (mod(mstj(3),5)==1) pw = p(i, 4)
      If (mod(mstj(3),5)==2) pw = sqrt(p(i,5)**2+(psi(1)*p(i,1)+psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
      If (mod(mstj(3),5)==3) pw = 1.
      Do j = 1, 3
        p(i, j) = p(i, j) - psi(j)*pw/pws
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
    End Do
  Else If (mod(mstj(3),5)==4) Then
    Do i = n + 1, n + njet
      k(i, 1) = 0
      Do j = 1, 5
        p(i, j) = 0.
      End Do
    End Do
    Do i = nsav + njet + 1, n
      ir1 = k(i, 3)
      ir2 = n + ir1 - nsav
      k(ir2, 1) = k(ir2, 1) + 1
      pls = (p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/(p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
      Do j = 1, 3
        p(ir2, j) = p(ir2, j) + p(i, j) - pls*p(ir1, j)
      End Do
      p(ir2, 4) = p(ir2, 4) + p(i, 4)
      p(ir2, 5) = p(ir2, 5) + pls
    End Do
    pss = 0.
    Do i = n + 1, n + njet
      If (k(i,1)/=0) pss = pss + p(i, 4)/(pecm*(0.8*p(i,5)+0.2))
    End Do
    Do i = nsav + njet + 1, n
      ir1 = k(i, 3)
      ir2 = n + ir1 - nsav
      pls = (p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/(p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
      Do j = 1, 3
        p(i, j) = p(i, j) - p(ir2, j)/k(ir2, 1) + (1./(p(ir2,5)*pss)-1.)*pls*p(ir1, j)
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
    End Do
  End If
  If (mod(mstj(3),5)/=0) Then
    pms = 0.
    pes = 0.
    pqs = 0.
    Do i = nsav + njet + 1, n
      pms = pms + p(i, 5)
      pes = pes + p(i, 4)
      pqs = pqs + p(i, 5)**2/p(i, 4)
    End Do
    If (pms>=pecm) Goto 150
    neco = 0
    440 neco = neco + 1
    pfac = (pecm-pqs)/(pes-pqs)
    pes = 0.
    pqs = 0.
    Do i = nsav + njet + 1, n
      Do j = 1, 3
        p(i, j) = pfac*p(i, j)
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      pes = pes + p(i, 4)
      pqs = pqs + p(i, 5)**2/p(i, 4)
    End Do
    If (neco<10 .And. abs(pecm-pes)>2E-6*pecm) Goto 440
  End If
  470 Do i = nsav + njet + 1, n
    If (mstu(16)/=2) k(i, 3) = nsav + 1
    If (mstu(16)==2) k(i, 3) = k(k(i,3), 3)
  End Do
  Do i = nsav + 1, nsav + njet
    i1 = k(i, 3)
    k(i1, 1) = k(i1, 1) + 10
    If (mstu(16)/=2) Then
      k(i1, 4) = nsav + 1
      k(i1, 5) = nsav + 1
    Else
      k(i1, 4) = k(i1, 4) - njet + 1
      k(i1, 5) = k(i1, 5) - njet + 1
      If (k(i1,5)<k(i1,4)) Then
        k(i1, 4) = 0
        k(i1, 5) = 0
      End If
    End If
  End Do
  nsav = nsav + 1
  k(nsav, 1) = 11
  k(nsav, 2) = 93
  k(nsav, 3) = ip
  k(nsav, 4) = nsav + 1
  k(nsav, 5) = n - njet + 1
  Do j = 1, 4
    p(nsav, j) = sngl(dps(j))
    v(nsav, j) = v(ip, j)
  End Do
  p(nsav, 5) = sqrt(sngl(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2)))
  v(nsav, 5) = 0.
  Do i = nsav + njet, n
    Do j = 1, 5
      k(i-njet+1, j) = k(i, j)
      p(i-njet+1, j) = p(i, j)
      v(i-njet+1, j) = v(i, j)
    End Do
  End Do
  n = n - njet + 1
  If (njet/=1) Call ludbrb(nsav+1, n, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))
  Do i = nsav + 1, n
    Do j = 1, 4
      v(i, j) = v(ip, j)
    End Do
  End Do
  Return
End Subroutine luindf
