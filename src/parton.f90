Subroutine parton(f, x1, x2, qq)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Common /njet/n, ipcrs
  Common /cmsflag/dshadow, ishadow
  Dimension f(2, 7)
  Save
  dlam = dble(hipr1(15))
  q0 = dble(hipr1(16))
  s = dlog(dlog(qq/dlam**2)/dlog(q0**2/dlam**2))
  If (ihpr2(7)==2) Goto 200
  at1 = 0.419D0 + 0.004D0*s - 0.007D0*s**2
  at2 = 3.460D0 + 0.724D0*s - 0.066D0*s**2
  gmud = 4.40D0 - 4.86D0*s + 1.33D0*s**2
  at3 = 0.763D0 - 0.237D0*s + 0.026D0*s**2
  at4 = 4.00D0 + 0.627D0*s - 0.019D0*s**2
  gmd = -0.421D0*s + 0.033D0*s**2
  cas = 1.265D0 - 1.132D0*s + 0.293D0*s**2
  as = -0.372D0*s - 0.029D0*s**2
  bs = 8.05D0 + 1.59D0*s - 0.153D0*s**2
  aphs = 6.31D0*s - 0.273D0*s**2
  btas = -10.5D0*s - 3.17D0*s**2
  gms = 14.7D0*s + 9.80D0*s**2
  cag = 1.56D0 - 1.71D0*s + 0.638D0*s**2
  ag = -0.949D0*s + 0.325D0*s**2
  bg = 6.0D0 + 1.44D0*s - 1.05D0*s**2
  aphg = 9.0D0 - 7.19D0*s + 0.255D0*s**2
  btag = -16.5D0*s + 10.9D0*s**2
  gmg = 15.3D0*s - 10.1D0*s**2
  Goto 300
  200 at1 = 0.374D0 + 0.014D0*s
  at2 = 3.33D0 + 0.753D0*s - 0.076D0*s**2
  gmud = 6.03D0 - 6.22D0*s + 1.56D0*s**2
  at3 = 0.761D0 - 0.232D0*s + 0.023D0*s**2
  at4 = 3.83D0 + 0.627D0*s - 0.019D0*s**2
  gmd = -0.418D0*s + 0.036D0*s**2
  cas = 1.67D0 - 1.92D0*s + 0.582D0*s**2
  as = -0.273D0*s - 0.164D0*s**2
  bs = 9.15D0 + 0.530D0*s - 0.763D0*s**2
  aphs = 15.7D0*s - 2.83D0*s**2
  btas = -101.0D0*s + 44.7D0*s**2
  gms = 223.0D0*s - 117.0D0*s**2
  cag = 0.879D0 - 0.971D0*s + 0.434D0*s**2
  ag = -1.16D0*s + 0.476D0*s**2
  bg = 4.0D0 + 1.23D0*s - 0.254D0*s**2
  aphg = 9.0D0 - 5.64D0*s - 0.817D0*s**2
  btag = -7.54D0*s + 5.50D0*s**2
  gmg = -0.596D0*s + 1.26D0*s**2
  300 b12 = dexp(gmre(at1)+gmre(at2+1.D0)-gmre(at1+at2+1.D0))
  b34 = dexp(gmre(at3)+gmre(at4+1.D0)-gmre(at3+at4+1.D0))
  cnud = 3.D0/b12/(1.D0+gmud*at1/(at1+at2+1.D0))
  cnd = 1.D0/b34/(1.D0+gmd*at3/(at3+at4+1.D0))
  fud1 = cnud*x1**at1*(1.D0-x1)**at2*(1.D0+gmud*x1)
  fs1 = cas*x1**as*(1.D0-x1)**bs*(1.D0+aphs*x1+btas*x1**2+gms*x1**3)
  f(1, 3) = cnd*x1**at3*(1.D0-x1)**at4*(1.D0+gmd*x1) + fs1/6.D0
  f(1, 1) = fud1 - f(1, 3) + fs1/3.D0
  f(1, 2) = fs1/6.D0
  f(1, 4) = fs1/6.D0
  f(1, 5) = fs1/6.D0
  f(1, 6) = fs1/6.D0
  f(1, 7) = cag*x1**ag*(1.D0-x1)**bg*(1.D0+aphg*x1+btag*x1**2+gmg*x1**3)
  fud2 = cnud*x2**at1*(1.D0-x2)**at2*(1.D0+gmud*x2)
  fs2 = cas*x2**as*(1.D0-x2)**bs*(1.D0+aphs*x2+btas*x2**2+gms*x2**3)
  f(2, 3) = cnd*x2**at3*(1.D0-x2)**at4*(1.D0+gmd*x2) + fs2/6.D0
  f(2, 1) = fud2 - f(2, 3) + fs2/3.D0
  f(2, 2) = fs2/6.D0
  f(2, 4) = fs2/6.D0
  f(2, 5) = fs2/6.D0
  f(2, 6) = fs2/6.D0
  f(2, 7) = cag*x2**ag*(1.D0-x2)**bg*(1.D0+aphg*x2+btag*x2**2+gmg*x2**3)
  If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
    aax = 1.193D0*dble(alog(float(ihnt2(1)))**0.16666666)
    rrx = aax*(x1**3-1.2D0*x1**2+0.21D0*x1) + 1.D0 + dble(1.079*(float(ihnt2(1))**0.33333333-1.0))/dble(alog(float(ihnt2(1))+1.0))*dsqrt(x1)*dexp(-x1**2/0.01D0)
    If (ishadow==1) rrx = 1.D0 + dshadow*(rrx-1.D0)
    If (ipcrs==1 .Or. ipcrs==3) rrx = dexp(-x1**2/0.01D0)
    If ((ipcrs==1 .Or. ipcrs==3) .And. ishadow==1) rrx = dexp(-x1**2/0.01D0)*dshadow
    Do i = 1, 7
      f(1, i) = rrx*f(1, i)
    End Do
  End If
  If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
    aax = 1.193D0*dble(alog(float(ihnt2(3)))**0.16666666)
    rrx = aax*(x2**3-1.2D0*x2**2+0.21D0*x2) + 1.D0 + dble(1.079*(float(ihnt2(3))**0.33333-1.0))/dble(alog(float(ihnt2(3))+1.0))*dsqrt(x2)*dexp(-x2**2/0.01D0)
    If (ishadow==1) rrx = 1.D0 + dshadow*(rrx-1.D0)
    If (ipcrs==2 .Or. ipcrs==3) rrx = dexp(-x2**2/0.01D0)
    If ((ipcrs==2 .Or. ipcrs==3) .And. ishadow==1) rrx = dexp(-x2**2/0.01D0)*dshadow
    Do i = 1, 7
      f(2, i) = rrx*f(2, i)
    End Do
  End If
  Return
End Subroutine parton
