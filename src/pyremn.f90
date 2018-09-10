Subroutine pyremn(ipu1, ipu2)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save /hparnt/
  Common /hstrng/nfp(300, 15), pphi(300, 15), nft(300, 15), pthi(300, 15)
  Save /hstrng/
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Dimension kflch(2), kflsp(2), chi(2), pms(6), is(2), robo(5)
  If (mint(43)==1) Then
    Do jt = 1, 2
      i = mint(83) + jt + 2
      k(i, 1) = 21
      k(i, 2) = k(i-2, 2)
      k(i, 3) = i - 2
      Do j = 1, 5
        p(i, j) = p(i-2, j)
      End Do
    End Do
  End If
  If (ipu1==0 .And. ipu2==0) Return
  isub = mint(1)
  ilep = 0
  If (ipu1==0) ilep = 1
  If (ipu2==0) ilep = 2
  If (isub==95) ilep = -1
  If (ilep==1) iq = mint(84) + 1
  If (ilep==2) iq = mint(84) + 2
  ip = max(ipu1, ipu2)
  ilepr = mint(83) + 5 - ilep
  ns = n
  110 Do jt = 1, 2
    i = mint(83) + jt + 2
    If (jt==1) ipu = ipu1
    If (jt==2) ipu = ipu2
    k(i, 1) = 21
    k(i, 3) = i - 2
    If (isub==95) Then
      k(i, 2) = 21
      shs = 0.
    Else If (mint(40+jt)==1 .And. ipu/=0) Then
      k(i, 2) = k(ipu, 2)
      p(i, 5) = p(ipu, 5)
      p(i, 1) = 0.
      p(i, 2) = 0.
      pms(jt) = p(i, 5)**2
    Else If (ipu/=0) Then
      k(i, 2) = k(ipu, 2)
      p(i, 5) = p(ipu, 5)
      rpt1 = 0.0
      rpt2 = 0.0
      ssw2 = (pphi(ihnt2(11),4)+pthi(ihnt2(12),4))**2 - (pphi(ihnt2(11),1)+pthi(ihnt2(12),1))**2 - (pphi(ihnt2(11),2)+pthi(ihnt2(12),2))**2 - (pphi(ihnt2(11),3)+pthi(ihnt2(12),3))**2
      If (ssw2<=4.0*parp(93)**2) Goto 1211
      If (ihpr2(5)<=0) Then
        120 If (mstp(91)<=0) Then
          pt = 0.
        Else If (mstp(91)==1) Then
          pt = parp(91)*sqrt(-log(rlu(0)))
        Else
          rpt1 = rlu(0)
          rpt2 = rlu(0)
          pt = -parp(92)*log(rpt1*rpt2)
        End If
        If (pt>parp(93)) Goto 120
        phi = paru(2)*rlu(0)
        rpt1 = pt*cos(phi)
        rpt2 = pt*sin(phi)
      Else If (ihpr2(5)==1) Then
        If (jt==1) jpt = nfp(ihnt2(11), 11)
        If (jt==2) jpt = nft(ihnt2(12), 11)
        1205 ptgs = parp(91)*sqrt(-log(rlu(0)))
        If (ptgs>parp(93)) Goto 1205
        phi = 2.0*hipr1(40)*rlu(0)
        rpt1 = ptgs*cos(phi)
        rpt2 = ptgs*sin(phi)
        Do iint = 1, jpt - 1
          pkcsq = parp(91)*sqrt(-log(rlu(0)))
          phi = 2.0*hipr1(40)*rlu(0)
          rpt1 = rpt1 + pkcsq*cos(phi)
          rpt2 = rpt2 + pkcsq*sin(phi)
        End Do
        If (rpt1**2+rpt2**2>=ssw2/4.0) Goto 1205
      End If
      1211 p(i, 1) = rpt1
      p(i, 2) = rpt2
      pms(jt) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
    Else
      k(i, 2) = k(iq, 2)
      q2 = vint(52)
      p(i, 5) = -sqrt(q2)
      pms(jt) = -q2
      shs = (1.-vint(43-jt))*q2/vint(43-jt) + vint(5-jt)**2
    End If
  End Do
  i1 = mint(83) + 3
  i2 = mint(83) + 4
  If (ilep==0) shs = vint(141)*vint(142)*vint(2) + (p(i1,1)+p(i2,1))**2 + (p(i1,2)+p(i2,2))**2
  shr = sqrt(max(0.,shs))
  If (ilep==0) Then
    If ((shs-pms(1)-pms(2))**2-4.*pms(1)*pms(2)<=0.) Goto 110
    p(i1, 4) = 0.5*(shr+(pms(1)-pms(2))/shr)
    p(i1, 3) = sqrt(max(0.,p(i1,4)**2-pms(1)))
    p(i2, 4) = shr - p(i1, 4)
    p(i2, 3) = -p(i1, 3)
  Else If (ilep==1) Then
    p(i1, 4) = p(iq, 4)
    p(i1, 3) = p(iq, 3)
    p(i2, 4) = p(ip, 4)
    p(i2, 3) = p(ip, 3)
  Else If (ilep==2) Then
    p(i1, 4) = p(ip, 4)
    p(i1, 3) = p(ip, 3)
    p(i2, 4) = p(iq, 4)
    p(i2, 3) = p(iq, 3)
  End If
  If (mint(43)==1) Return
  If (ilep==0) Then
    robo(3) = (p(i1,1)+p(i2,1))/shr
    robo(4) = (p(i1,2)+p(i2,2))/shr
    Call ludbrb(i1, i2, 0., 0., -dble(robo(3)), -dble(robo(4)), 0D0)
    robo(2) = ulangl(p(i1,1), p(i1,2))
    Call ludbrb(i1, i2, 0., -robo(2), 0D0, 0D0, 0D0)
    robo(1) = ulangl(p(i1,3), p(i1,1))
    Call ludbrb(i1, i2, -robo(1), 0., 0D0, 0D0, 0D0)
    nmax = max(mint(52), ipu1, ipu2)
    Call ludbrb(i1, nmax, robo(1), robo(2), dble(robo(3)), dble(robo(4)), 0D0)
    robo(5) = max(-0.999999, min(0.999999,(vint(141)-vint(142))/(vint(141)+vint(142))))
    Call ludbrb(i1, nmax, 0., 0., 0D0, 0D0, dble(robo(5)))
  End If
  If (ilep<=0) Then
    If (mstp(81)<=0 .Or. mstp(82)<=0 .Or. isub==95) Then
      vint(151) = 0.
      vint(152) = 0.
    End If
    peh = p(i1, 4) + p(i2, 4) + 0.5*vint(1)*(vint(151)+vint(152))
    pzh = p(i1, 3) + p(i2, 3) + 0.5*vint(1)*(vint(151)-vint(152))
    shh = (vint(1)-peh)**2 - (p(i1,1)+p(i2,1))**2 - (p(i1,2)+p(i2,2))**2 - pzh**2
    pmmin = p(mint(83)+1, 5) + p(mint(83)+2, 5) + ulmass(k(i1,2)) + ulmass(k(i2,2))
    If (shr>=vint(1) .Or. shh<=(pmmin+parp(111))**2) Then
      mint(51) = 1
      Return
    End If
    shr = sqrt(shh+(p(i1,1)+p(i2,1))**2+(p(i1,2)+p(i2,2))**2)
  Else
    pei = p(iq, 4) + p(ip, 4)
    pzi = p(iq, 3) + p(ip, 3)
    pms(ilep) = max(0., pei**2-pzi**2)
    pmmin = p(ilepr-2, 5) + ulmass(k(ilepr,2)) + sqrt(pms(ilep))
    If (shr<=pmmin+parp(111)) Then
      mint(51) = 1
      Return
    End If
  End If
  140 i = ns
  Do jt = 1, 2
    If (jt==ilep) Goto 190
    If (jt==1) ipu = ipu1
    If (jt==2) ipu = ipu2
    Call pyspli(mint(10+jt), mint(12+jt), kflch(jt), kflsp(jt))
    i = i + 1
    is(jt) = i
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
    k(i, 1) = 3
    k(i, 2) = kflsp(jt)
    k(i, 3) = mint(83) + jt
    p(i, 5) = ulmass(k(i,2))
    kfls = (3-kchg(lucomp(kflsp(jt)),2)*isign(1,kflsp(jt)))/2
    k(i, kfls+3) = ipu
    k(ipu, 6-kfls) = mod(k(ipu,6-kfls), mstu(5)) + mstu(5)*i
    If (kflch(jt)==0) Then
      p(i, 1) = -p(mint(83)+jt+2, 1)
      p(i, 2) = -p(mint(83)+jt+2, 2)
      pms(jt) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
    Else
      Call luptdi(1, p(i,1), p(i,2))
      pms(jt+2) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
      i = i + 1
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
      k(i, 1) = 1
      k(i, 2) = kflch(jt)
      k(i, 3) = mint(83) + jt
      p(i, 5) = ulmass(k(i,2))
      p(i, 1) = -p(mint(83)+jt+2, 1) - p(i-1, 1)
      p(i, 2) = -p(mint(83)+jt+2, 2) - p(i-1, 2)
      pms(jt+4) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
      imb = 1
      If (mod(mint(10+jt)/1000,10)/=0) imb = 2
      If (iabs(kflch(jt))<=10 .Or. kflch(jt)==21) Then
        chik = parp(92+2*imb)
        If (mstp(92)<=1) Then
          If (imb==1) chi(jt) = rlu(0)
          If (imb==2) chi(jt) = 1. - sqrt(rlu(0))
        Else If (mstp(92)==2) Then
          chi(jt) = 1. - rlu(0)**(1./(1.+chik))
        Else If (mstp(92)==3) Then
          cut = 2.*0.3/vint(1)
          170 chi(jt) = rlu(0)**2
          If ((chi(jt)**2/(chi(jt)**2+cut**2))**0.25*(1.-chi(jt))**chik<rlu(0)) Goto 170
        Else
          cut = 2.*0.3/vint(1)
          cutr = (1.+sqrt(1.+cut**2))/cut
          180 chir = cut*cutr**rlu(0)
          chi(jt) = (chir**2-cut**2)/(2.*chir)
          If ((1.-chi(jt))**chik<rlu(0)) Goto 180
        End If
      Else
        If (mstp(92)<=1) Then
          If (imb==1) chi(jt) = rlu(0)
          If (imb==2) chi(jt) = 1. - sqrt(rlu(0))
        Else
          chi(jt) = 1. - rlu(0)**(1./(1.+parp(93+2*imb)))
        End If
        If (mod(kflch(jt)/1000,10)/=0) chi(jt) = 1. - chi(jt)
      End If
      pms(jt) = pms(jt+4)/chi(jt) + pms(jt+2)/(1.-chi(jt))
      kfls = kchg(lucomp(kflch(jt)), 2)*isign(1, kflch(jt))
      If (kfls/=0) Then
        k(i, 1) = 3
        kfls = (3-kfls)/2
        k(i, kfls+3) = ipu
        k(ipu, 6-kfls) = mod(k(ipu,6-kfls), mstu(5)) + mstu(5)*i
      End If
    End If
  190 End Do
  If (shr<=sqrt(pms(1))+sqrt(pms(2))) Goto 140
  n = i
  Do jt = 1, 2
    If (jt==ilep) Goto 200
    pe = 0.5*(shr+(pms(jt)-pms(3-jt))/shr)
    pz = sqrt(pe**2-pms(jt))
    If (kflch(jt)==0) Then
      p(is(jt), 4) = pe
      p(is(jt), 3) = pz*(-1)**(jt-1)
    Else
      pw1 = chi(jt)*(pe+pz)
      p(is(jt)+1, 4) = 0.5*(pw1+pms(jt+4)/pw1)
      p(is(jt)+1, 3) = 0.5*(pw1-pms(jt+4)/pw1)*(-1)**(jt-1)
      p(is(jt), 4) = pe - p(is(jt)+1, 4)
      p(is(jt), 3) = pz*(-1)**(jt-1) - p(is(jt)+1, 3)
    End If
  200 End Do
  If (ilep<=0) Then
    Call ludbrb(ns+1, n, 0., 0., 0D0, 0D0, -dble(pzh/(vint(1)-peh)))
  Else
    nmax = max(ip, mint(52))
    pef = shr - pe
    pzf = pz*(-1)**(ilep-1)
    pt2 = p(ilepr, 1)**2 + p(ilepr, 2)**2
    phipt = ulangl(p(ilepr,1), p(ilepr,2))
    Call ludbrb(mint(84)+1, nmax, 0., -phipt, 0D0, 0D0, 0D0)
    rqp = p(iq, 3)*(pt2+pei**2) - p(iq, 4)*pei*pzi
    sinth = p(iq, 4)*sqrt(pt2*(pt2+pei**2)/(rqp**2+pt2*p(iq,4)**2*pzi**2))*sign(1., -rqp)
    Call ludbrb(mint(84)+1, nmax, asin(sinth), 0., 0D0, 0D0, 0D0)
    betax = (-pei*pzi*sinth+sqrt(pt2*(pt2+pei**2-(pzi*sinth)**2)))/(pt2+pei**2)
    Call ludbrb(mint(84)+1, nmax, 0., 0., dble(betax), 0D0, 0D0)
    Call ludbrb(mint(84)+1, nmax, 0., phipt, 0D0, 0D0, 0D0)
    pem = p(iq, 4) + p(ip, 4)
    pzm = p(iq, 3) + p(ip, 3)
    betaz = (-pem*pzm+pzf*sqrt(pzf**2+pem**2-pzm**2))/(pzf**2+pem**2)
    Call ludbrb(mint(84)+1, nmax, 0., 0., 0D0, 0D0, dble(betaz))
    Call ludbrb(i1, i2, asin(sinth), 0., dble(betax), 0D0, 0D0)
    Call ludbrb(i1, i2, 0., phipt, 0D0, 0D0, dble(betaz))
  End If
  Return
End Subroutine pyremn
