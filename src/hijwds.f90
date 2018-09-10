Subroutine hijwds(ia, idh, xhigh)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /wood/r, d, fnorm, w
  Dimension iaa(20), rr(20), dd(20), ww(20)
  External rwdsax, wdsax
  Save
  Data iaa/2, 4, 12, 16, 27, 32, 40, 56, 63, 93, 184, 197, 208, 7*0./
  Data rr/0.01, .964, 2.355, 2.608, 2.84, 3.458, 3.766, 3.971, 4.214, 4.87, 6.51, 6.38, 6.624, 7*0./
  Data dd/0.5882, .322, .522, .513, .569, .61, .586, .5935, .586, .573, .535, .535, .549, 7*0./
  Data ww/0.0, .517, -0.149, -0.051, 0., -0.208, -0.161, 13*0./
  a = ia
  d = 0.54
  r = 1.19*a**(1./3.) - 1.61*a**(-1./3.)
  w = 0.
  Do i = 1, 13
    If (ia==iaa(i)) Then
      r = rr(i)
      d = dd(i)
      w = ww(i)
    End If
  End Do
  fnorm = 1.0
  xlow = 0.
  xhigh = r + 12.*d
  If (w<-0.01) Then
    If (xhigh>r/sqrt(abs(w))) xhigh = r/sqrt(abs(w))
  End If
  fgaus = gauss1(rwdsax, xlow, xhigh, 0.001)
  fnorm = 1./fgaus
  If (idh==1) Then
    hint1(72) = r
    hint1(73) = d
    hint1(74) = w
    hint1(75) = fnorm/4.0/hipr1(40)
  Else If (idh==2) Then
    hint1(76) = r
    hint1(77) = d
    hint1(78) = w
    hint1(79) = fnorm/4.0/hipr1(40)
  End If
  Call hifun(idh, xlow, xhigh, rwdsax)
  Return
End Subroutine hijwds
