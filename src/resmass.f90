Function resmass(kf)
  Parameter (arho=0.775, aomega=0.783, aeta=0.548, aks=0.894, aphi=1.019, adelta=1.232)
  Parameter (wrho=0.149, womega=0.00849, weta=1.30E-6, wks=0.0498, wphi=0.00426, wdelta=0.118)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  If (kf==113 .Or. abs(kf)==213) Then
     amass = arho
     wid = wrho
  Else If (kf==221) Then
     amass = aeta
     wid = weta
  Else If (kf==223) Then
     amass = aomega
     wid = womega
  Else If (abs(kf)==313 .Or. abs(kf)==323) Then
     amass = aks
     wid = wks
  Else If (kf==333) Then
     amass = aphi
     wid = wphi
  Else If (abs(kf)==1114 .Or. abs(kf)==2114 .Or. abs(kf)==2214 .Or. abs(kf)==2224) Then
     amass = adelta
     wid = wdelta
  End If
  dmin = amass - 2*wid
  dmax = amass + 2*wid
  If (amass==adelta) dmin = 1.078
  fm = 1.
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  fmass = (amass*wid)**2/((dm**2-amass**2)**2+(amass*wid)**2)
  If ((ranart(nseed)>fmass/fm) .And. (ntry1<=10)) Goto 10
  resmass = dm
  Return
End Function resmass
