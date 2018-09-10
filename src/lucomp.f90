Function lucomp(kf)
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  lucomp = 0
  kfa = iabs(kf)
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  kfls = mod(kfa, 10)
  kflr = mod(kfa/10000, 10)
  If (kfa==0 .Or. kfa>=100000) Then
  Else If (kfa<=100) Then
    lucomp = kfa
    If (kf<0 .And. kchg(kfa,3)==0) lucomp = 0
  Else If (kfls==0) Then
    If (kf==130) lucomp = 221
    If (kf==310) lucomp = 222
    If (kfa==210) lucomp = 281
    If (kfa==2110) lucomp = 282
    If (kfa==2210) lucomp = 283
  Else If (kfa-10000*kflr<1000) Then
    If (kflb==0 .Or. kflb==9 .Or. kflc==0 .Or. kflc==9) Then
    Else If (kflb<kflc) Then
    Else If (kf<0 .And. kflb==kflc) Then
    Else If (kflb==kflc) Then
      If (kflr==0 .And. kfls==1) Then
        lucomp = 110 + kflb
      Else If (kflr==0 .And. kfls==3) Then
        lucomp = 130 + kflb
      Else If (kflr==1 .And. kfls==3) Then
        lucomp = 150 + kflb
      Else If (kflr==1 .And. kfls==1) Then
        lucomp = 170 + kflb
      Else If (kflr==2 .And. kfls==3) Then
        lucomp = 190 + kflb
      Else If (kflr==0 .And. kfls==5) Then
        lucomp = 210 + kflb
      End If
    Else If (kflb<=5 .And. kflc<=3) Then
      If (kflr==0 .And. kfls==1) Then
        lucomp = 100 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==0 .And. kfls==3) Then
        lucomp = 120 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==1 .And. kfls==3) Then
        lucomp = 140 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==1 .And. kfls==1) Then
        lucomp = 160 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==2 .And. kfls==3) Then
        lucomp = 180 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==0 .And. kfls==5) Then
        lucomp = 200 + ((kflb-1)*(kflb-2))/2 + kflc
      End If
    Else If ((kfls==1 .And. kflr<=1) .Or. (kfls==3 .And. kflr<=2) .Or. (kfls==5 .And. kflr==0)) Then
      lucomp = 80 + kflb
    End If
  Else If ((kflr==0 .Or. kflr==1) .And. kflc==0) Then
    If (kfls/=1 .And. kfls/=3) Then
    Else If (kfla==9 .Or. kflb==0 .Or. kflb==9) Then
    Else If (kfla<kflb) Then
    Else If (kfls==1 .And. kfla==kflb) Then
    Else
      lucomp = 90
    End If
  Else If (kflr==0 .And. kfls==2) Then
    If (kfla==9 .Or. kflb==0 .Or. kflb==9 .Or. kflc==9) Then
    Else If (kfla<=kflc .Or. kfla<kflb) Then
    Else If (kfla>=6 .Or. kflb>=4 .Or. kflc>=4) Then
      lucomp = 80 + kfla
    Else If (kflb<kflc) Then
      lucomp = 300 + ((kfla+1)*kfla*(kfla-1))/6 + (kflc*(kflc-1))/2 + kflb
    Else
      lucomp = 330 + ((kfla+1)*kfla*(kfla-1))/6 + (kflb*(kflb-1))/2 + kflc
    End If
  Else If (kflr==0 .And. kfls==4) Then
    If (kfla==9 .Or. kflb==0 .Or. kflb==9 .Or. kflc==9) Then
    Else If (kfla<kflb .Or. kflb<kflc) Then
    Else If (kfla>=6 .Or. kflb>=4) Then
      lucomp = 80 + kfla
    Else
      lucomp = 360 + ((kfla+1)*kfla*(kfla-1))/6 + (kflb*(kflb-1))/2 + kflc
    End If
  End If
  Return
End Function lucomp
