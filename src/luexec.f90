Subroutine luexec
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension ps(2, 6)
  mstu(24) = 0
  If (mstu(12)>=1) Call lulist(0)
  mstu(31) = mstu(31) + 1
  mstu(1) = 0
  mstu(2) = 0
  mstu(3) = 0
  mcons = 1
  nsav = n
  Do i = 1, 2
    Do j = 1, 6
      ps(i, j) = 0.
    End Do
  End Do
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 120
    Do j = 1, 4
      ps(1, j) = ps(1, j) + p(i, j)
    End Do
    ps(1, 6) = ps(1, 6) + luchge(k(i,2))
  120 End Do
  paru(21) = ps(1, 4)
  Call luprep(0)
  mbe = 0
  130 mbe = mbe + 1
  ip = 0
  140 ip = ip + 1
  kc = 0
  If (k(ip,1)>0 .And. k(ip,1)<=10) kc = lucomp(k(ip,2))
  If (kc==0) Then
  Else If (kchg(kc,2)==0) Then
    If (mstj(21)>=1 .And. mdcy(kc,1)>=1) Then
      If (mstj(51)<=0 .Or. mbe==2 .Or. pmas(kc,2)>=parj(91) .Or. iabs(k(ip,2))==311) Call ludecy(ip)
    End If
    If (mstj(92)>0) Then
      ip1 = mstj(92)
      qmax = sqrt(max(0.,(p(ip1,4)+p(ip1+1,4))**2-(p(ip1,1)+p(ip1+1,1))**2-(p(ip1,2)+p(ip1+1,2))**2-(p(ip1,3)+p(ip1+1,3))**2))
      Call lushow(ip1, ip1+1, qmax)
      Call luprep(ip1)
      mstj(92) = 0
    Else If (mstj(92)<0) Then
      ip1 = -mstj(92)
      pip5 = p(ip, 5)
      Call lushow(ip1, -3, pip5)
      Call luprep(ip1)
      mstj(92) = 0
    End If
  Else If (k(ip,1)==1 .Or. k(ip,1)==2) Then
    mfrag = mstj(1)
    If (mfrag>=1 .And. k(ip,1)==1) mfrag = 2
    If (mstj(21)>=2 .And. k(ip,1)==2 .And. n>ip) Then
      If (k(ip+1,1)==1 .And. k(ip+1,3)==k(ip,3) .And. k(ip,3)>0 .And. k(ip,3)<ip) Then
        If (kchg(lucomp(k(k(ip,3),2)),2)==0) mfrag = min(1, mfrag)
      End If
    End If
    If (mfrag==1) Then
      Call lustrf(ip)
    End If
    If (mfrag==2) Call luindf(ip)
    If (mfrag==2 .And. k(ip,1)==1) mcons = 0
    If (mfrag==2 .And. (mstj(3)<=0 .Or. mod(mstj(3),5)==0)) mcons = 0
  End If
  If (mstu(24)/=0 .And. mstu(21)>=2) Then
  Else If (ip<n .And. n<mstu(4)-20-mstu(32)) Then
    Goto 140
  Else If (ip<n) Then
    Call luerrm(11, '(LUEXEC:) no more memory left in LUJETS')
  End If
  If (mbe==1 .And. mstj(51)>=1) Then
    Call luboei(nsav)
    Goto 130
  End If
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 160
    Do j = 1, 4
      ps(2, j) = ps(2, j) + p(i, j)
    End Do
    ps(2, 6) = ps(2, 6) + luchge(k(i,2))
  160 End Do
  pdev = (abs(ps(2,1)-ps(1,1))+abs(ps(2,2)-ps(1,2))+abs(ps(2,3)-ps(1,3))+abs(ps(2,4)-ps(1,4)))/(1.+abs(ps(2,4))+abs(ps(1,4)))
  If (mcons==1 .And. pdev>paru(11)) Call luerrm(15, '(LUEXEC:) four-momentum was not conserved')
  If (mcons==1 .And. abs(ps(2,6)-ps(1,6))>0.1) Call luerrm(15, '(LUEXEC:) charge was not conserved')
  Return
End Subroutine luexec
