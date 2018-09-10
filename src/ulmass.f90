Function ulmass(kf)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  ulmass = 0.
  kfa = iabs(kf)
  kc = lucomp(kf)
  If (kc==0) Return
  parf(106) = pmas(6, 1)
  parf(107) = pmas(7, 1)
  parf(108) = pmas(8, 1)
  If ((mstj(93)==1 .Or. mstj(93)==2) .And. kfa<=10) Then
    ulmass = parf(100+kfa)
    If (mstj(93)==2) ulmass = max(0., ulmass-parf(121))
  Else If (kfa<=100 .Or. kc<=80 .Or. kc>100) Then
    ulmass = pmas(kc, 1)
  Else
    kfla = mod(kfa/1000, 10)
    kflb = mod(kfa/100, 10)
    kflc = mod(kfa/10, 10)
    kfls = mod(kfa, 10)
    kflr = mod(kfa/10000, 10)
    pma = parf(100+kfla)
    pmb = parf(100+kflb)
    pmc = parf(100+kflc)
    If (kfla==0 .And. kflr==0 .And. kfls<=3) Then
      If (kfls==1) pmspl = -3./(pmb*pmc)
      If (kfls>=3) pmspl = 1./(pmb*pmc)
      ulmass = parf(111) + pmb + pmc + parf(113)*parf(101)**2*pmspl
    Else If (kfla==0) Then
      kmul = 2
      If (kfls==1) kmul = 3
      If (kflr==2) kmul = 4
      If (kfls==5) kmul = 5
      ulmass = parf(113+kmul) + pmb + pmc
    Else If (kflc==0) Then
      If (kfls==1) pmspl = -3./(pma*pmb)
      If (kfls==3) pmspl = 1./(pma*pmb)
      ulmass = 2.*parf(112)/3. + pma + pmb + parf(114)*parf(101)**2*pmspl
      If (mstj(93)==1) ulmass = pma + pmb
      If (mstj(93)==2) ulmass = max(0., ulmass-parf(122)-2.*parf(112)/3.)
    Else
      If (kfls==2 .And. kfla==kflb) Then
        pmspl = 1./(pma*pmb) - 2./(pma*pmc) - 2./(pmb*pmc)
      Else If (kfls==2 .And. kflb>=kflc) Then
        pmspl = -2./(pma*pmb) - 2./(pma*pmc) + 1./(pmb*pmc)
      Else If (kfls==2) Then
        pmspl = -3./(pmb*pmc)
      Else
        pmspl = 1./(pma*pmb) + 1./(pma*pmc) + 1./(pmb*pmc)
      End If
      ulmass = parf(112) + pma + pmb + pmc + parf(114)*parf(101)**2*pmspl
    End If
  End If
  If (mstj(24)>=1 .And. pmas(kc,2)>1E-4) Then
    If (mstj(24)==1 .Or. (mstj(24)==2 .And. kfa>100)) Then
      ulmass = ulmass + 0.5*pmas(kc, 2)*tan((2.*rlu(0)-1.)*atan(2.*pmas(kc,3)/pmas(kc,2)))
    Else
      pm0 = ulmass
      pmlow = atan((max(0.,pm0-pmas(kc,3))**2-pm0**2)/(pm0*pmas(kc,2)))
      pmupp = atan((pm0+pmas(kc,3))**2-pm0**2)/(pm0*pmas(kc,2))
      ulmass = sqrt(max(0.,pm0**2+pm0*pmas(kc,2)*tan(pmlow+(pmupp-pmlow)*rlu(0))))
    End If
  End If
  mstj(93) = 0
  Return
End Function ulmass
