Subroutine luname(kf, chau)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Character chau*16
  chau = ' '
  kfa = iabs(kf)
  kc = lucomp(kf)
  If (kc==0) Return
  kq = luchge(kf)
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  kfls = mod(kfa, 10)
  kflr = mod(kfa/10000, 10)
  If (kfa<=100 .Or. (kfa>100 .And. kc>100)) Then
    chau = chaf(kc)
    len = 0
    Do lem = 1, 8
      If (chau(lem:lem)/=' ') len = lem
    End Do
  Else If (kflc==0) Then
    chau(1:2) = chaf(kfla)(1:1) // chaf(kflb)(1:1)
    If (kfls==1) chau(3:4) = '_0'
    If (kfls==3) chau(3:4) = '_1'
    len = 4
  Else If (kfla==0) Then
    If (kflb==5) chau(1:1) = 'B'
    If (kflb==6) chau(1:1) = 'T'
    If (kflb==7) chau(1:1) = 'L'
    If (kflb==8) chau(1:1) = 'H'
    len = 1
    If (kflr==0 .And. kfls==1) Then
    Else If (kflr==0 .And. kfls==3) Then
      chau(2:2) = '*'
      len = 2
    Else If (kflr==1 .And. kfls==3) Then
      chau(2:3) = '_1'
      len = 3
    Else If (kflr==1 .And. kfls==1) Then
      chau(2:4) = '*_0'
      len = 4
    Else If (kflr==2) Then
      chau(2:4) = '*_1'
      len = 4
    Else If (kfls==5) Then
      chau(2:4) = '*_2'
      len = 4
    End If
    If (kflc>=3 .And. kflr==0 .And. kfls<=3) Then
      chau(len+1:len+2) = '_' // chaf(kflc)(1:1)
      len = len + 2
    Else If (kflc>=3) Then
      chau(len+1:len+1) = chaf(kflc)(1:1)
      len = len + 1
    End If
  Else
    If (kflb<=2 .And. kflc<=2) Then
      chau = 'Sigma '
      If (kflc>kflb) chau = 'Lambda'
      If (kfls==4) chau = 'Sigma*'
      len = 5
      If (chau(6:6)/=' ') len = 6
    Else If (kflb<=2 .Or. kflc<=2) Then
      chau = 'Xi '
      If (kfla>kflb .And. kflb>kflc) chau = 'Xi'''
      If (kfls==4) chau = 'Xi*'
      len = 2
      If (chau(3:3)/=' ') len = 3
    Else
      chau = 'Omega '
      If (kfla>kflb .And. kflb>kflc) chau = 'Omega'''
      If (kfls==4) chau = 'Omega*'
      len = 5
      If (chau(6:6)/=' ') len = 6
    End If
    chau(len+1:len+2) = '_' // chaf(kfla)(1:1)
    len = len + 2
    If (kflb>=kflc .And. kflc>=4) Then
      chau(len+1:len+2) = chaf(kflb)(1:1) // chaf(kflc)(1:1)
      len = len + 2
    Else If (kflb>=kflc .And. kflb>=4) Then
      chau(len+1:len+1) = chaf(kflb)(1:1)
      len = len + 1
    Else If (kflc>kflb .And. kflb>=4) Then
      chau(len+1:len+2) = chaf(kflc)(1:1) // chaf(kflb)(1:1)
      len = len + 2
    Else If (kflc>kflb .And. kflc>=4) Then
      chau(len+1:len+1) = chaf(kflc)(1:1)
      len = len + 1
    End If
  End If
  If (kf>0 .Or. len==0) Then
  Else If (kfa>10 .And. kfa<=40 .And. kq/=0) Then
  Else If (kfa==89 .Or. (kfa>=91 .And. kfa<=99)) Then
  Else If (kfa>100 .And. kfla==0 .And. kq/=0) Then
  Else If (mstu(15)<=1) Then
    chau(len+1:len+1) = '~'
    len = len + 1
  Else
    chau(len+1:len+3) = 'bar'
    len = len + 3
  End If
  If (kq==6) chau(len+1:len+2) = '++'
  If (kq==-6) chau(len+1:len+2) = '--'
  If (kq==3) chau(len+1:len+1) = '+'
  If (kq==-3) chau(len+1:len+1) = '-'
  If (kq==0 .And. (kfa<=22 .Or. len==0)) Then
  Else If (kq==0 .And. (kfa>=81 .And. kfa<=100)) Then
  Else If (kfa>100 .And. kfla==0 .And. kflb==kflc .And. kflb/=1) Then
  Else If (kq==0) Then
    chau(len+1:len+1) = '0'
  End If
  Return
End Subroutine luname
