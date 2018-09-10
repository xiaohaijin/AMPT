Subroutine pydiff
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Do jt = 1, mstp(126) + 10
    i = mint(83) + jt
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do
  n = mint(84)
  mint(3) = 0
  mint(21) = 0
  mint(22) = 0
  mint(23) = 0
  mint(24) = 0
  mint(4) = 4
  Do jt = 1, 2
    i = mint(83) + jt
    k(i, 1) = 21
    k(i, 2) = mint(10+jt)
    p(i, 5) = vint(2+jt)
    p(i, 3) = vint(5)*(-1)**(jt+1)
    p(i, 4) = sqrt(p(i,3)**2+p(i,5)**2)
  End Do
  mint(6) = 2
  isub = mint(1)
  sqlam = (vint(2)-vint(63)-vint(64))**2 - 4.*vint(63)*vint(64)
  pz = sqrt(sqlam)/(2.*vint(1))
  Do jt = 1, 2
    i = mint(83) + jt
    pe = (vint(2)+vint(62+jt)-vint(65-jt))/(2.*vint(1))
    If (mint(16+jt)<=0) Then
      n = n + 1
      k(n, 1) = 1
      k(n, 2) = k(i, 2)
      k(n, 3) = i + 2
      p(n, 3) = pz*(-1)**(jt+1)
      p(n, 4) = pe
      p(n, 5) = p(i, 5)
    Else If (mstp(101)==1) Then
      n = n + 2
      k(n-1, 1) = 2
      k(n, 1) = 1
      k(n-1, 3) = i + 2
      k(n, 3) = i + 2
      Call pyspli(k(i,2), 21, k(n,2), k(n-1,2))
      p(n-1, 5) = ulmass(k(n-1,2))
      p(n, 5) = ulmass(k(n,2))
      sqlam = (vint(62+jt)-p(n-1,5)**2-p(n,5)**2)**2 - 4.*p(n-1, 5)**2*p(n, 5)**2
      p(n-1, 3) = (pe*sqrt(sqlam)+pz*(vint(62+jt)+p(n-1,5)**2-p(n,5)**2))/(2.*vint(62+jt))*(-1)**(jt+1)
      p(n-1, 4) = sqrt(p(n-1,3)**2+p(n-1,5)**2)
      p(n, 3) = pz*(-1)**(jt+1) - p(n-1, 3)
      p(n, 4) = sqrt(p(n,3)**2+p(n,5)**2)
    Else
      n = n + 3
      k(n-2, 1) = 2
      k(n-1, 1) = 2
      k(n, 1) = 1
      k(n-2, 3) = i + 2
      k(n-1, 3) = i + 2
      k(n, 3) = i + 2
      Call pyspli(k(i,2), 21, k(n,2), k(n-2,2))
      k(n-1, 2) = 21
      p(n-2, 5) = ulmass(k(n-2,2))
      p(n-1, 5) = 0.
      p(n, 5) = ulmass(k(n,2))
      120 imb = 1
      If (mod(k(i,2)/1000,10)/=0) imb = 2
      chik = parp(92+2*imb)
      If (mstp(92)<=1) Then
        If (imb==1) chi = rlu(0)
        If (imb==2) chi = 1. - sqrt(rlu(0))
      Else If (mstp(92)==2) Then
        chi = 1. - rlu(0)**(1./(1.+chik))
      Else If (mstp(92)==3) Then
        cut = 2.*0.3/vint(1)
        130 chi = rlu(0)**2
        If ((chi**2/(chi**2+cut**2))**0.25*(1.-chi)**chik<rlu(0)) Goto 130
      Else
        cut = 2.*0.3/vint(1)
        cutr = (1.+sqrt(1.+cut**2))/cut
        140 chir = cut*cutr**rlu(0)
        chi = (chir**2-cut**2)/(2.*chir)
        If ((1.-chi)**chik<rlu(0)) Goto 140
      End If
      If (chi<p(n,5)**2/vint(62+jt) .Or. chi>1.-p(n-2,5)**2/vint(62+jt)) Goto 120
      sqm = p(n-2, 5)**2/(1.-chi) + p(n, 5)**2/chi
      If ((sqrt(sqm)+parj(32))**2>=vint(62+jt)) Goto 120
      pzi = (pe*(vint(62+jt)-sqm)+pz*(vint(62+jt)+sqm))/(2.*vint(62+jt))
      pei = sqrt(pzi**2+sqm)
      pqqp = (1.-chi)*(pei+pzi)
      p(n-2, 3) = 0.5*(pqqp-p(n-2,5)**2/pqqp)*(-1)**(jt+1)
      p(n-2, 4) = sqrt(p(n-2,3)**2+p(n-2,5)**2)
      p(n-1, 3) = (pz-pzi)*(-1)**(jt+1)
      p(n-1, 4) = abs(p(n-1,3))
      p(n, 3) = pzi*(-1)**(jt+1) - p(n-2, 3)
      p(n, 4) = sqrt(p(n,3)**2+p(n,5)**2)
    End If
    k(i+2, 1) = 21
    If (mint(16+jt)==0) k(i+2, 2) = mint(10+jt)
    If (mint(16+jt)/=0) k(i+2, 2) = 10*(mint(10+jt)/10)
    k(i+2, 3) = i
    p(i+2, 3) = pz*(-1)**(jt+1)
    p(i+2, 4) = pe
    p(i+2, 5) = sqrt(vint(62+jt))
  End Do
  Call ludbrb(mint(83)+3, n, acos(vint(23)), vint(24), 0D0, 0D0, 0D0)
  Return
End Subroutine pydiff
