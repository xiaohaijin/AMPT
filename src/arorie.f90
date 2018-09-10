Subroutine arorie(s, x1, x3, jl)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /rndf77/nseed
  Save
  w = sqrt(s)
  x2 = 2. - x1 - x3
  e1 = .5*x1*w
  e3 = .5*x3*w
  p1 = sqrt(e1**2-p(jl,5)**2)
  p3 = sqrt(e3**2-p(jl+1,5)**2)
  cbet = 1.
  If (p1>0. .And. p3>0.) cbet = (p(jl,5)**2+p(jl+1,5)**2+2.*e1*e3-s*(1.-x2))/(2.*p1*p3)
  If (abs(cbet)>1.0) cbet = max(-1., min(1.,cbet))
  bet = acos(cbet)
  If (p1>=p3) Then
    psi = .5*ulangl(p1**2+p3**2*cos(2.*bet), -p3**2*sin(2.*bet))
    pt1 = p1*sin(psi)
    pz1 = p1*cos(psi)
    pt3 = p3*sin(psi+bet)
    pz3 = p3*cos(psi+bet)
  Else If (p3>p1) Then
    psi = .5*ulangl(p3**2+p1**2*cos(2.*bet), -p1**2*sin(2.*bet))
    pt1 = p1*sin(bet+psi)
    pz1 = -p1*cos(bet+psi)
    pt3 = p3*sin(psi)
    pz3 = -p3*cos(psi)
  End If
  del = 2.0*hipr1(40)*ranart(nseed)
  p(jl, 4) = e1
  p(jl, 1) = pt1*sin(del)
  p(jl, 2) = -pt1*cos(del)
  p(jl, 3) = pz1
  p(jl+2, 4) = e3
  p(jl+2, 1) = pt3*sin(del)
  p(jl+2, 2) = -pt3*cos(del)
  p(jl+2, 3) = pz3
  p(jl+1, 4) = w - e1 - e3
  p(jl+1, 1) = -p(jl, 1) - p(jl+2, 1)
  p(jl+1, 2) = -p(jl, 2) - p(jl+2, 2)
  p(jl+1, 3) = -p(jl, 3) - p(jl+2, 3)
  Return
End Subroutine arorie
