Subroutine attrad(ierror)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Common /rndf77/nseed
  Save
  ierror = 0
  40 sm = 0.
  jl = 1
  Do i = 1, n - 1
    s = 2.*(p(i,4)*p(i+1,4)-p(i,1)*p(i+1,1)-p(i,2)*p(i+1,2)-p(i,3)*p(i+1,3)) + p(i, 5)**2 + p(i+1, 5)**2
    If (s<0.) s = 0.
    wp = sqrt(s) - 1.5*(p(i,5)+p(i+1,5))
    If (wp>sm) Then
      pbt1 = p(i, 1) + p(i+1, 1)
      pbt2 = p(i, 2) + p(i+1, 2)
      pbt3 = p(i, 3) + p(i+1, 3)
      pbt4 = p(i, 4) + p(i+1, 4)
      btt = (pbt1**2+pbt2**2+pbt3**2)/pbt4**2
      If (btt>=1.0-1.0E-10) Goto 30
      If ((i/=1 .Or. i/=n-1) .And. (k(i,2)/=21 .And. k(i+1,2)/=21)) Goto 30
      jl = i
      sm = wp
    End If
  30 End Do
  s = (sm+1.5*(p(jl,5)+p(jl+1,5)))**2
  If (sm<hipr1(5)) Goto 2
  If (jl+1==n) Goto 190
  Do j = n, jl + 2, -1
    k(j+1, 1) = k(j, 1)
    k(j+1, 2) = k(j, 2)
    Do m = 1, 5
      p(j+1, m) = p(j, m)
    End Do
  End Do
  190 n = n + 1
  p1 = p(jl, 1) + p(jl+1, 1)
  p2 = p(jl, 2) + p(jl+1, 2)
  p3 = p(jl, 3) + p(jl+1, 3)
  p4 = p(jl, 4) + p(jl+1, 4)
  bex = -p1/p4
  bey = -p2/p4
  bez = -p3/p4
  imin = jl
  imax = jl + 1
  Call atrobo(0., 0., bex, bey, bez, imin, imax, ierror)
  If (ierror/=0) Return
  cth = p(jl, 3)/sqrt(p(jl,4)**2-p(jl,5)**2)
  If (abs(cth)>1.0) cth = max(-1., min(1.,cth))
  theta = acos(cth)
  phi = ulangl(p(jl,1), p(jl,2))
  Call atrobo(0., -phi, 0., 0., 0., imin, imax, ierror)
  Call atrobo(-theta, 0., 0., 0., 0., imin, imax, ierror)
  1 Call ar3jet(s, x1, x3, jl)
  Call arorie(s, x1, x3, jl)
  If (hidat(2)>0.0) Then
    ptg1 = sqrt(p(jl,1)**2+p(jl,2)**2)
    ptg2 = sqrt(p(jl+1,1)**2+p(jl+1,2)**2)
    ptg3 = sqrt(p(jl+2,1)**2+p(jl+2,2)**2)
    ptg = max(ptg1, ptg2, ptg3)
    If (ptg>hidat(2)) Then
      fmfact = exp(-(ptg**2-hidat(2)**2)/hipr1(2)**2)
      If (ranart(nseed)>fmfact) Goto 1
    End If
  End If
  imin = jl
  imax = jl + 2
  Call atrobo(theta, phi, -bex, -bey, -bez, imin, imax, ierror)
  If (ierror/=0) Return
  k(jl+2, 1) = k(jl+1, 1)
  k(jl+2, 2) = k(jl+1, 2)
  k(jl+2, 3) = k(jl+1, 3)
  k(jl+2, 4) = k(jl+1, 4)
  k(jl+2, 5) = k(jl+1, 5)
  p(jl+2, 5) = p(jl+1, 5)
  k(jl+1, 1) = 2
  k(jl+1, 2) = 21
  k(jl+1, 3) = 0
  k(jl+1, 4) = 0
  k(jl+1, 5) = 0
  p(jl+1, 5) = 0.
  If (sm>=hipr1(5)) Goto 40
  2 k(1, 1) = 2
  k(1, 3) = 0
  k(1, 4) = 0
  k(1, 5) = 0
  k(n, 1) = 1
  k(n, 3) = 0
  k(n, 4) = 0
  k(n, 5) = 0
  Return
End Subroutine attrad
