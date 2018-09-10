Subroutine hijfrg(jtp, ntp, ierror)
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /rndf77/nseed
  Common /anim/nevent, isoft, isflag, izpc
  Save
  ierror = 0
  Call luedit(0)
  n = 0
  If (ntp==3) Then
    isg = jtp
    n = njsg(isg)
    Do i = 1, njsg(isg)
      k(i, 1) = k1sg(isg, i)
      k(i, 2) = k2sg(isg, i)
      p(i, 1) = pxsg(isg, i)
      p(i, 2) = pysg(isg, i)
      p(i, 3) = pzsg(isg, i)
      p(i, 4) = pesg(isg, i)
      p(i, 5) = pmsg(isg, i)
    End Do
    If (isoft/=2 .Or. isflag/=0) Call luexec
    Return
  End If
  If (ntp==2) Goto 200
  If (jtp>ihnt2(1)) Return
  If (nfp(jtp,5)/=3 .And. nfp(jtp,3)/=0 .And. npj(jtp)==0 .And. nfp(jtp,10)==0) Goto 1000
  If (nfp(jtp,15)==-1) Then
    kf1 = nfp(jtp, 2)
    kf2 = nfp(jtp, 1)
    pq21 = pp(jtp, 6)
    pq22 = pp(jtp, 7)
    pq11 = pp(jtp, 8)
    pq12 = pp(jtp, 9)
    am1 = pp(jtp, 15)
    am2 = pp(jtp, 14)
  Else
    kf1 = nfp(jtp, 1)
    kf2 = nfp(jtp, 2)
    pq21 = pp(jtp, 8)
    pq22 = pp(jtp, 9)
    pq11 = pp(jtp, 6)
    pq12 = pp(jtp, 7)
    am1 = pp(jtp, 14)
    am2 = pp(jtp, 15)
  End If
  pb1 = pq11 + pq21
  pb2 = pq12 + pq22
  pb3 = pp(jtp, 3)
  pecm = pp(jtp, 5)
  btz = pb3/pp(jtp, 4)
  If ((abs(pb1-pp(jtp,1))>0.01 .Or. abs(pb2-pp(jtp,2))>0.01) .And. ihpr2(10)/=0) Write (6, *) '  Pt of Q and QQ do not sum to the total', jtp, ntp, pq11, pq21, pb1, '*', pq12, pq22, pb2, '*', pp(jtp, 1), pp(jtp, 2)
  Goto 300
  200 If (jtp>ihnt2(3)) Return
  If (nft(jtp,5)/=3 .And. nft(jtp,3)/=0 .And. ntj(jtp)==0 .And. nft(jtp,10)==0) Goto 1200
  If (nft(jtp,15)==1) Then
    kf1 = nft(jtp, 1)
    kf2 = nft(jtp, 2)
    pq11 = pt(jtp, 6)
    pq12 = pt(jtp, 7)
    pq21 = pt(jtp, 8)
    pq22 = pt(jtp, 9)
    am1 = pt(jtp, 14)
    am2 = pt(jtp, 15)
  Else
    kf1 = nft(jtp, 2)
    kf2 = nft(jtp, 1)
    pq11 = pt(jtp, 8)
    pq12 = pt(jtp, 9)
    pq21 = pt(jtp, 6)
    pq22 = pt(jtp, 7)
    am1 = pt(jtp, 15)
    am2 = pt(jtp, 14)
  End If
  pb1 = pq11 + pq21
  pb2 = pq12 + pq22
  pb3 = pt(jtp, 3)
  pecm = pt(jtp, 5)
  btz = pb3/pt(jtp, 4)
  If ((abs(pb1-pt(jtp,1))>0.01 .Or. abs(pb2-pt(jtp,2))>0.01) .And. ihpr2(10)/=0) Write (6, *) '  Pt of Q and QQ do not sum to the total', jtp, ntp, pq11, pq21, pb1, '*', pq12, pq22, pb2, '*', pt(jtp, 1), pt(jtp, 2)
  300 If (pecm<hipr1(1)) Then
    ierror = 1
    If (ihpr2(10)==0) Return
    Write (6, *) ' ECM=', pecm, ' energy of the string is too small'
    Write (6, *) 'JTP,NTP,pq=', jtp, ntp, pq11, pq12, pq21, pq22
    Return
  End If
  amt = pecm**2 + pb1**2 + pb2**2
  amt1 = am1**2 + pq11**2 + pq12**2
  amt2 = am2**2 + pq21**2 + pq22**2
  pzcm = sqrt(abs(amt**2+amt1**2+amt2**2-2.0*amt*amt1-2.0*amt*amt2-2.0*amt1*amt2))/2.0/sqrt(amt)
  k(1, 1) = 2
  k(1, 2) = kf1
  p(1, 1) = pq11
  p(1, 2) = pq12
  p(1, 3) = pzcm
  p(1, 4) = sqrt(amt1+pzcm**2)
  p(1, 5) = am1
  k(2, 1) = 1
  k(2, 2) = kf2
  p(2, 1) = pq21
  p(2, 2) = pq22
  p(2, 3) = -pzcm
  p(2, 4) = sqrt(amt2+pzcm**2)
  p(2, 5) = am2
  n = 2
  Call hirobo(0.0, 0.0, 0.0, 0.0, btz)
  jetot = 0
  If ((pq21**2+pq22**2)>(pq11**2+pq12**2)) Then
    pmax1 = p(2, 1)
    pmax2 = p(2, 2)
    pmax3 = p(2, 3)
  Else
    pmax1 = p(1, 1)
    pmax2 = p(1, 2)
    pmax3 = p(1, 3)
  End If
  If (ntp==1) Then
    pp(jtp, 10) = pmax1
    pp(jtp, 11) = pmax2
    pp(jtp, 12) = pmax3
  Else If (ntp==2) Then
    pt(jtp, 10) = pmax1
    pt(jtp, 11) = pmax2
    pt(jtp, 12) = pmax3
  End If
  If (ntp==1 .And. npj(jtp)/=0) Then
    jetot = npj(jtp)
    iex = 0
    If ((abs(kf1)>1000 .And. kf1<0) .Or. (abs(kf1)<1000 .And. kf1>0)) iex = 1
    Do i = n, 2, -1
      Do j = 1, 5
        ii = npj(jtp) + i
        k(ii, j) = k(i, j)
        p(ii, j) = p(i, j)
        v(ii, j) = v(i, j)
      End Do
    End Do
    Do i = 1, npj(jtp)
      Do j = 1, 5
        k(i+1, j) = 0
        v(i+1, j) = 0
      End Do
      i0 = i
      If (iex==1 .And. (isoft/=2 .Or. isflag/=0)) i0 = npj(jtp) - i + 1
      kk1 = kfpj(jtp, i0)
      k(i+1, 1) = 2
      k(i+1, 2) = kk1
      If (kk1/=21 .And. kk1/=0) k(i+1, 1) = 1 + (abs(kk1)+(2*iex-1)*kk1)/2/abs(kk1)
      p(i+1, 1) = pjpx(jtp, i0)
      p(i+1, 2) = pjpy(jtp, i0)
      p(i+1, 3) = pjpz(jtp, i0)
      p(i+1, 4) = pjpe(jtp, i0)
      p(i+1, 5) = pjpm(jtp, i0)
    End Do
    n = n + npj(jtp)
  Else If (ntp==2 .And. ntj(jtp)/=0) Then
    jetot = ntj(jtp)
    iex = 1
    If ((abs(kf2)>1000 .And. kf2<0) .Or. (abs(kf2)<1000 .And. kf2>0)) iex = 0
    Do i = n, 2, -1
      Do j = 1, 5
        ii = ntj(jtp) + i
        k(ii, j) = k(i, j)
        p(ii, j) = p(i, j)
        v(ii, j) = v(i, j)
      End Do
    End Do
    Do i = 1, ntj(jtp)
      Do j = 1, 5
        k(i+1, j) = 0
        v(i+1, j) = 0
      End Do
      i0 = i
      If (iex==1 .And. (isoft/=2 .Or. isflag/=0)) i0 = ntj(jtp) - i + 1
      kk1 = kftj(jtp, i0)
      k(i+1, 1) = 2
      k(i+1, 2) = kk1
      If (kk1/=21 .And. kk1/=0) k(i+1, 1) = 1 + (abs(kk1)+(2*iex-1)*kk1)/2/abs(kk1)
      p(i+1, 1) = pjtx(jtp, i0)
      p(i+1, 2) = pjty(jtp, i0)
      p(i+1, 3) = pjtz(jtp, i0)
      p(i+1, 4) = pjte(jtp, i0)
      p(i+1, 5) = pjtm(jtp, i0)
    End Do
    n = n + ntj(jtp)
  End If
  If (ihpr2(1)>0 .And. ranart(nseed)<=hidat(3)) Then
    hdat20 = hidat(2)
    hpr150 = hipr1(5)
    If (ihpr2(8)==0 .And. ihpr2(3)==0 .And. ihpr2(9)==0) hidat(2) = 2.0
    If (hint1(1)>=1000.0 .And. jetot==0) Then
      hidat(2) = 3.0
      hipr1(5) = 5.0
    End If
    Call attrad(ierror)
    hidat(2) = hdat20
    hipr1(5) = hpr150
  Else If (jetot==0 .And. ihpr2(1)>0 .And. hint1(1)>=1000.0 .And. ranart(nseed)<=0.8) Then
    hdat20 = hidat(2)
    hpr150 = hipr1(5)
    hidat(2) = 3.0
    hipr1(5) = 5.0
    If (ihpr2(8)==0 .And. ihpr2(3)==0 .And. ihpr2(9)==0) hidat(2) = 2.0
    Call attrad(ierror)
    hidat(2) = hdat20
    hipr1(5) = hpr150
  End If
  If (ierror/=0) Return
  If (isoft/=2 .Or. isflag/=0) Call luexec
  Return
  1000 n = 1
  k(1, 1) = 1
  k(1, 2) = nfp(jtp, 3)
  Do jj = 1, 5
    p(1, jj) = pp(jtp, jj)
  End Do
  If (isoft/=2 .Or. isflag/=0) Call luexec
  Return
  1200 n = 1
  k(1, 1) = 1
  k(1, 2) = nft(jtp, 3)
  Do jj = 1, 5
    p(1, jj) = pt(jtp, jj)
  End Do
  If (isoft/=2 .Or. isflag/=0) Call luexec
  Return
End Subroutine hijfrg
