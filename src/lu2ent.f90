Subroutine lu2ent(ip, kf1, kf2, pecm)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  mstu(28) = 0
  If (mstu(12)>=1) Call lulist(0)
  ipa = max(1, iabs(ip))
  If (ipa>mstu(4)-1) Call luerrm(21, '(LU2ENT:) writing outside LUJETS memory')
  kc1 = lucomp(kf1)
  kc2 = lucomp(kf2)
  If (kc1==0 .Or. kc2==0) Call luerrm(12, '(LU2ENT:) unknown flavour code')
  pm1 = 0.
  If (mstu(10)==1) pm1 = p(ipa, 5)
  If (mstu(10)>=2) pm1 = ulmass(kf1)
  pm2 = 0.
  If (mstu(10)==1) pm2 = p(ipa+1, 5)
  If (mstu(10)>=2) pm2 = ulmass(kf2)
  Do i = ipa, ipa + 1
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do
  kq1 = kchg(kc1, 2)*isign(1, kf1)
  kq2 = kchg(kc2, 2)*isign(1, kf2)
  If (kq1+kq2/=0 .And. kq1+kq2/=4) Call luerrm(2, '(LU2ENT:) unphysical flavour combination')
  k(ipa, 2) = kf1
  k(ipa+1, 2) = kf2
  If (ip>=0) Then
    k(ipa, 1) = 1
    If (kq1/=0 .And. kq2/=0) k(ipa, 1) = 2
    k(ipa+1, 1) = 1
  Else
    If (kq1==0 .Or. kq2==0) Call luerrm(2, '(LU2ENT:) requested flavours can not develop parton shower')
    k(ipa, 1) = 3
    k(ipa+1, 1) = 3
    k(ipa, 4) = mstu(5)*(ipa+1)
    k(ipa, 5) = k(ipa, 4)
    k(ipa+1, 4) = mstu(5)*ipa
    k(ipa+1, 5) = k(ipa+1, 4)
  End If
  If (pecm<=pm1+pm2) Call luerrm(13, '(LU2ENT:) energy smaller than sum of masses')
  pa = sqrt(max(0.,(pecm**2-pm1**2-pm2**2)**2-(2.*pm1*pm2)**2))/(2.*pecm)
  p(ipa, 3) = pa
  p(ipa, 4) = sqrt(pm1**2+pa**2)
  p(ipa, 5) = pm1
  p(ipa+1, 3) = -pa
  p(ipa+1, 4) = sqrt(pm2**2+pa**2)
  p(ipa+1, 5) = pm2
  n = ipa + 1
  If (ip==0) Call luexec
  Return
End Subroutine lu2ent
