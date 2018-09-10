!> @brief
!!       初始化整个事件的设定(一些读进去的事件设定参数)
Subroutine hijset(efrm, frame, proj, targ, iap, izp, iat, izt)
  Character frame*4, proj*4, targ*4, eframe*4
  Double Precision dd1, dd2, dd3, dd4
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  External fnkick, fnkc2, fnstru, fnstrm, fnstrs
  Save
  Call title
  ihnt2(1) = iap
  ihnt2(2) = izp
  ihnt2(3) = iat
  ihnt2(4) = izt
  ihnt2(5) = 0
  ihnt2(6) = 0
  hint1(8) = max(ulmass(2112), ulmass(2212))
  hint1(9) = hint1(8)
  If (proj/='A') Then
    If (proj=='P') Then
      ihnt2(5) = 2212
    Else If (proj=='PBAR') Then
      ihnt2(5) = -2212
    Else If (proj=='PI+') Then
      ihnt2(5) = 211
    Else If (proj=='PI-') Then
      ihnt2(5) = -211
    Else If (proj=='K+') Then
      ihnt2(5) = 321
    Else If (proj=='K-') Then
      ihnt2(5) = -321
    Else If (proj=='N') Then
      ihnt2(5) = 2112
    Else If (proj=='NBAR') Then
      ihnt2(5) = -2112
    Else
      Write (6, *) proj, 'wrong or unavailable proj name'
      Stop
    End If
    hint1(8) = ulmass(ihnt2(5))
  End If
  If (targ/='A') Then
    If (targ=='P') Then
      ihnt2(6) = 2212
    Else If (targ=='PBAR') Then
      ihnt2(6) = -2212
    Else If (targ=='PI+') Then
      ihnt2(6) = 211
    Else If (targ=='PI-') Then
      ihnt2(6) = -211
    Else If (targ=='K+') Then
      ihnt2(6) = 321
    Else If (targ=='K-') Then
      ihnt2(6) = -321
    Else If (targ=='N') Then
      ihnt2(6) = 2112
    Else If (targ=='NBAR') Then
      ihnt2(6) = -2112
    Else
      Write (6, *) targ, 'wrong or unavailable targ name'
      Stop
    End If
    hint1(9) = ulmass(ihnt2(6))
  End If
  If (ihpr2(12)>0) Then
    Call lugive('MDCY(C221,1)=0')
    Call lugive('MDCY(C313,1)=0')
    Call lugive('MDCY(C-313,1)=0')
    Call lugive('MDCY(C323,1)=0')
    Call lugive('MDCY(C-323,1)=0')
    Call lugive('MDCY(C311,1)=0')
    Call lugive('MDCY(C-311,1)=0')
    Call lugive('MDCY(C1114,1)=0')
    Call lugive('MDCY(C2114,1)=0')
    Call lugive('MDCY(C2214,1)=0')
    Call lugive('MDCY(C2224,1)=0')
    Call lugive('MDCY(C-1114,1)=0')
    Call lugive('MDCY(C-2114,1)=0')
    Call lugive('MDCY(C-2214,1)=0')
    Call lugive('MDCY(C-2224,1)=0')
    Call lugive('MDCY(C213,1)=0')
    Call lugive('MDCY(C-213,1)=0')
    Call lugive('MDCY(C113,1)=0')
    Call lugive('MDCY(C223,1)=0')
    Call lugive('MDCY(C333,1)=0')
    Call lugive('MDCY(C111,1)=0')
    Call lugive('MDCY(C310,1)=0')
    Call lugive('MDCY(C411,1)=0;MDCY(C-411,1)=0')
    Call lugive('MDCY(C421,1)=0;MDCY(C-421,1)=0')
    Call lugive('MDCY(C431,1)=0;MDCY(C-431,1)=0')
    Call lugive('MDCY(C511,1)=0;MDCY(C-511,1)=0')
    Call lugive('MDCY(C521,1)=0;MDCY(C-521,1)=0')
    Call lugive('MDCY(C531,1)=0;MDCY(C-531,1)=0')
    Call lugive('MDCY(C3122,1)=0;MDCY(C-3122,1)=0')
    Call lugive('MDCY(C3112,1)=0;MDCY(C-3112,1)=0')
    Call lugive('MDCY(C3212,1)=0;MDCY(C-3212,1)=0')
    Call lugive('MDCY(C3222,1)=0;MDCY(C-3222,1)=0')
    Call lugive('MDCY(C3312,1)=0;MDCY(C-3312,1)=0')
    Call lugive('MDCY(C3322,1)=0;MDCY(C-3322,1)=0')
    Call lugive('MDCY(C3334,1)=0;MDCY(C-3334,1)=0')
    Call lugive('MDCY(C441,1)=0')
    Call lugive('MDCY(C443,1)=0')
    Call lugive('MDCY(C413,1)=0;MDCY(C-413,1)=0')
    Call lugive('MDCY(C423,1)=0;MDCY(C-423,1)=0')
    Call lugive('MDCY(C433,1)=0;MDCY(C-433,1)=0')
    Call lugive('MDCY(C4112,1)=0;MDCY(C-4112,1)=0')
    Call lugive('MDCY(C4114,1)=0;MDCY(C-4114,1)=0')
    Call lugive('MDCY(C4122,1)=0;MDCY(C-4122,1)=0')
    Call lugive('MDCY(C4212,1)=0;MDCY(C-4212,1)=0')
    Call lugive('MDCY(C4214,1)=0;MDCY(C-4214,1)=0')
    Call lugive('MDCY(C4222,1)=0;MDCY(C-4222,1)=0')
    Call lugive('MDCY(C4224,1)=0;MDCY(C-4224,1)=0')
    Call lugive('MDCY(C4132,1)=0;MDCY(C-4132,1)=0')
    Call lugive('MDCY(C4312,1)=0;MDCY(C-4312,1)=0')
    Call lugive('MDCY(C4314,1)=0;MDCY(C-4314,1)=0')
    Call lugive('MDCY(C4232,1)=0;MDCY(C-4232,1)=0')
    Call lugive('MDCY(C4322,1)=0;MDCY(C-4322,1)=0')
    Call lugive('MDCY(C4324,1)=0;MDCY(C-4324,1)=0')
    Call lugive('MDCY(C4332,1)=0;MDCY(C-4332,1)=0')
    Call lugive('MDCY(C4334,1)=0;MDCY(C-4334,1)=0')
    Call lugive('MDCY(C551,1)=0')
    Call lugive('MDCY(C553,1)=0')
    Call lugive('MDCY(C513,1)=0;MDCY(C-513,1)=0')
    Call lugive('MDCY(C523,1)=0;MDCY(C-523,1)=0')
    Call lugive('MDCY(C533,1)=0;MDCY(C-533,1)=0')
    Call lugive('MDCY(C5112,1)=0;MDCY(C-5112,1)=0')
    Call lugive('MDCY(C5114,1)=0;MDCY(C-5114,1)=0')
    Call lugive('MDCY(C5122,1)=0;MDCY(C-5122,1)=0')
    Call lugive('MDCY(C5212,1)=0;MDCY(C-5212,1)=0')
    Call lugive('MDCY(C5214,1)=0;MDCY(C-5214,1)=0')
    Call lugive('MDCY(C5222,1)=0;MDCY(C-5222,1)=0')
    Call lugive('MDCY(C5224,1)=0;MDCY(C-5224,1)=0')
  End If
  mstu(12) = 0
  mstu(21) = 1
  If (ihpr2(10)==0) Then
    mstu(22) = 0
    mstu(25) = 0
    mstu(26) = 0
  End If
  mstj(12) = ihpr2(11)
  parj(21) = hipr1(2)
  rkp = hipr1(4)*(2+hipr1(3))/parj(42)/(2+parj(41))
  parj(2) = parj(2)**(1./rkp)
  parj(21) = parj(21)*sqrt(rkp)
  hipr1(2) = parj(21)
  parj(2) = min(parj(2), 0.4)
  If (frame=='LAB') Then
    dd1 = dble(efrm)
    dd2 = dble(hint1(8))
    dd3 = dble(hint1(9))
    hint1(1) = sqrt(hint1(8)**2+2.0*hint1(9)*efrm+hint1(9)**2)
    dd4 = dsqrt(dd1**2-dd2**2)/(dd1+dd3)
    hint1(2) = sngl(dd4)
    hint1(3) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    dd4 = dsqrt(dd1**2-dd2**2)/dd1
    hint1(4) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    hint1(5) = 0.0
    hint1(6) = efrm
    hint1(7) = hint1(9)
  Else If (frame=='CMS') Then
    hint1(1) = efrm
    hint1(2) = 0.0
    hint1(3) = 0.0
    dd1 = dble(hint1(1))
    dd2 = dble(hint1(8))
    dd3 = dble(hint1(9))
    dd4 = dsqrt(1.D0-4.D0*dd2**2/dd1**2)
    hint1(4) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    dd4 = dsqrt(1.D0-4.D0*dd3**2/dd1**2)
    hint1(5) = -0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    hint1(6) = hint1(1)/2.0
    hint1(7) = hint1(1)/2.0
  End If
  If (ihnt2(1)>1) Then
    Call hijwds(ihnt2(1), 1, rmax)
    hipr1(34) = rmax
  End If
  If (ihnt2(3)>1) Then
    Call hijwds(ihnt2(3), 2, rmax)
    hipr1(35) = rmax
  End If
  i = 0
  20 i = i + 1
  If (i==10) Goto 30
  If (hidat0(10,i)<=hint1(1)) Goto 20
  30 If (i==1) i = 2
  Do j = 1, 9
    hidat(j) = hidat0(j, i-1) + (hidat0(j,i)-hidat0(j,i-1))*(hint1(1)-hidat0(10,i-1))/(hidat0(10,i)-hidat0(10,i-1))
  End Do
  hipr1(31) = hidat(5)
  hipr1(30) = 2.0*hidat(5)
  Call hijcrs
  If (ihpr2(5)/=0) Then
    Call hifun(3, 0.0, 36.0, fnkick)
  End If
  Call hifun(7, 0.0, 6.0, fnkc2)
  Call hifun(4, 0.0, 1.0, fnstru)
  Call hifun(5, 0.0, 1.0, fnstrm)
  Call hifun(6, 0.0, 1.0, fnstrs)
  eframe = 'Ecm'
  If (frame=='LAB') eframe = 'Elab'
  Write (6, 100) eframe, efrm, proj, ihnt2(1), ihnt2(2), targ, ihnt2(3), ihnt2(4)
  Return
  100 Format (10X, '**************************************************'/10X, '*', 48X, '*'/10X, '*         HIJING has been initialized at         *'/10X, '*', 13X, A4, '= ', F10.2, ' GeV/n', 13X, '*'/10X, '*', 48X, '*'/10X, '*', 8X, 'for ', A4, '(', I3, ',', I3, ')', ' + ', A4, '(', I3, ',', I3, ')', 7X, '*'/10X, '**************************************************')
End Subroutine hijset
