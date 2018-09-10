Subroutine pyklim(ilim)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  isub = mint(1)
  If (isub==96) Goto 110
  sqm3 = vint(63)
  sqm4 = vint(64)
  If (ilim/=1) Then
    tau = vint(21)
    rm3 = sqm3/(tau*vint(2))
    rm4 = sqm4/(tau*vint(2))
    be34 = sqrt((1.-rm3-rm4)**2-4.*rm3*rm4)
  End If
  pthmin = ckin(3)
  If (min(sqm3,sqm4)<ckin(6)**2) pthmin = max(ckin(3), ckin(5))
  If (ilim==0) Then
    yst = vint(22)
    cth = vint(23)
    taup = vint(26)
    If (iset(isub)<=2) Then
      x1 = sqrt(tau)*exp(yst)
      x2 = sqrt(tau)*exp(-yst)
    Else
      x1 = sqrt(taup)*exp(yst)
      x2 = sqrt(taup)*exp(-yst)
    End If
    xf = x1 - x2
    If (tau*vint(2)<ckin(1)**2) mint(51) = 1
    If (ckin(2)>=0. .And. tau*vint(2)>ckin(2)**2) mint(51) = 1
    If (x1<ckin(21) .Or. x1>ckin(22)) mint(51) = 1
    If (x2<ckin(23) .Or. x2>ckin(24)) mint(51) = 1
    If (xf<ckin(25) .Or. xf>ckin(26)) mint(51) = 1
    If (yst<ckin(7) .Or. yst>ckin(8)) mint(51) = 1
    If (iset(isub)==2 .Or. iset(isub)==4) Then
      pth = 0.5*be34*sqrt(tau*vint(2)*(1.-cth**2))
      y3 = yst + 0.5*log((1.+rm3-rm4+be34*cth)/(1.+rm3-rm4-be34*cth))
      y4 = yst + 0.5*log((1.+rm4-rm3-be34*cth)/(1.+rm4-rm3+be34*cth))
      ylarge = max(y3, y4)
      ysmall = min(y3, y4)
      etalar = 10.
      etasma = -10.
      sth = sqrt(1.-cth**2)
      If (sth<1.E-6) Goto 100
      expet3 = ((1.+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth+sqrt(((1.+rm3-rm4)*cosh(yst)+be34*sinh(yst)*cth)**2-4.*rm3))/(be34*sth)
      expet4 = ((1.-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth+sqrt(((1.-rm3+rm4)*cosh(yst)-be34*sinh(yst)*cth)**2-4.*rm4))/(be34*sth)
      eta3 = log(min(1.E10,max(1.E-10,expet3)))
      eta4 = log(min(1.E10,max(1.E-10,expet4)))
      etalar = max(eta3, eta4)
      etasma = min(eta3, eta4)
      100 cts3 = ((1.+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth)/sqrt(((1.+rm3-rm4)*cosh(yst)+be34*sinh(yst)*cth)**2-4.*rm3)
      cts4 = ((1.-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth)/sqrt(((1.-rm3+rm4)*cosh(yst)-be34*sinh(yst)*cth)**2-4.*rm4)
      ctslar = max(cts3, cts4)
      ctssma = min(cts3, cts4)
      If (pth<pthmin) mint(51) = 1
      If (ckin(4)>=0. .And. pth>ckin(4)) mint(51) = 1
      If (ylarge<ckin(9) .Or. ylarge>ckin(10)) mint(51) = 1
      If (ysmall<ckin(11) .Or. ysmall>ckin(12)) mint(51) = 1
      If (etalar<ckin(13) .Or. etalar>ckin(14)) mint(51) = 1
      If (etasma<ckin(15) .Or. etasma>ckin(16)) mint(51) = 1
      If (ctslar<ckin(17) .Or. ctslar>ckin(18)) mint(51) = 1
      If (ctssma<ckin(19) .Or. ctssma>ckin(20)) mint(51) = 1
      If (cth<ckin(27) .Or. cth>ckin(28)) mint(51) = 1
    End If
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      If (taup*vint(2)<ckin(31)**2) mint(51) = 1
      If (ckin(32)>=0. .And. taup*vint(2)>ckin(32)**2) mint(51) = 1
    End If
  Else If (ilim==1) Then
    taumn0 = 0.
    taumx0 = 1.
    taumn1 = ckin(1)**2/vint(2)
    taumx1 = 1.
    If (ckin(2)>=0.) taumx1 = ckin(2)**2/vint(2)
    tm3 = sqrt(sqm3+pthmin**2)
    tm4 = sqrt(sqm4+pthmin**2)
    ydcosh = 1.
    If (ckin(9)>ckin(12)) ydcosh = cosh(ckin(9)-ckin(12))
    taumn2 = (tm3**2+2.*tm3*tm4*ydcosh+tm4**2)/vint(2)
    taumx2 = 1.
    cth2mn = min(ckin(27)**2, ckin(28)**2)
    cth2mx = max(ckin(27)**2, ckin(28)**2)
    taumn3 = 0.
    If (ckin(27)*ckin(28)>0.) taumn3 = (sqrt(sqm3+pthmin**2/(1.-cth2mn))+sqrt(sqm4+pthmin**2/(1.-cth2mn)))**2/vint(2)
    taumx3 = 1.
    If (ckin(4)>=0. .And. cth2mx<1.) taumx3 = (sqrt(sqm3+ckin(4)**2/(1.-cth2mx))+sqrt(sqm4+ckin(4)**2/(1.-cth2mx)))**2/vint(2)
    taumn4 = ckin(21)*ckin(23)
    taumx4 = ckin(22)*ckin(24)
    taumn5 = 0.
    taumx5 = max(1.-ckin(25), 1.+ckin(26))
    vint(11) = max(taumn0, taumn1, taumn2, taumn3, taumn4, taumn5)
    vint(31) = min(taumx0, taumx1, taumx2, taumx3, taumx4, taumx5)
    If (mint(43)==1 .And. (iset(isub)==1 .Or. iset(isub)==2)) Then
      vint(11) = 0.99999
      vint(31) = 1.00001
    End If
    If (vint(31)<=vint(11)) mint(51) = 1
  Else If (ilim==2) Then
    If (iset(isub)==3 .Or. iset(isub)==4) tau = vint(26)
    taurt = sqrt(tau)
    ystmn0 = log(taurt)
    ystmx0 = -ystmn0
    ystmn1 = ckin(7)
    ystmx1 = ckin(8)
    ystmn2 = log(max(tau,ckin(21))/taurt)
    ystmx2 = log(max(tau,ckin(22))/taurt)
    ystmn3 = -log(max(tau,ckin(24))/taurt)
    ystmx3 = -log(max(tau,ckin(23))/taurt)
    yepmn4 = 0.5*abs(ckin(25))/taurt
    ystmn4 = sign(log(sqrt(1.+yepmn4**2)+yepmn4), ckin(25))
    yepmx4 = 0.5*abs(ckin(26))/taurt
    ystmx4 = sign(log(sqrt(1.+yepmx4**2)+yepmx4), ckin(26))
    yepsmn = (rm3-rm4)*sinh(ckin(9)-ckin(11))
    yepsmx = (rm3-rm4)*sinh(ckin(10)-ckin(12))
    ydifmn = abs(log(sqrt(1.+yepsmn**2)-yepsmn))
    ydifmx = abs(log(sqrt(1.+yepsmx**2)-yepsmx))
    ystmn5 = 0.5*(ckin(9)+ckin(11)-ydifmn)
    ystmx5 = 0.5*(ckin(10)+ckin(12)+ydifmx)
    cthlim = sqrt(1.-4.*pthmin**2/(be34*tau*vint(2)))
    rzmn = be34*max(ckin(27), -cthlim)
    rzmx = be34*min(ckin(28), cthlim)
    yex3mx = (1.+rm3-rm4+rzmx)/max(1E-10, 1.+rm3-rm4-rzmx)
    yex4mx = (1.+rm4-rm3-rzmn)/max(1E-10, 1.+rm4-rm3+rzmn)
    yex3mn = max(1E-10, 1.+rm3-rm4+rzmn)/(1.+rm3-rm4-rzmn)
    yex4mn = max(1E-10, 1.+rm4-rm3-rzmx)/(1.+rm4-rm3+rzmx)
    ystmn6 = ckin(9) - 0.5*log(max(yex3mx,yex4mx))
    ystmx6 = ckin(12) - 0.5*log(min(yex3mn,yex4mn))
    vint(12) = max(ystmn0, ystmn1, ystmn2, ystmn3, ystmn4, ystmn5, ystmn6)
    vint(32) = min(ystmx0, ystmx1, ystmx2, ystmx3, ystmx4, ystmx5, ystmx6)
    If (mint(43)==1) Then
      vint(12) = -0.00001
      vint(32) = 0.00001
    Else If (mint(43)==2) Then
      vint(12) = 0.99999*ystmx0
      vint(32) = 1.00001*ystmx0
    Else If (mint(43)==3) Then
      vint(12) = -1.00001*ystmx0
      vint(32) = -0.99999*ystmx0
    End If
    If (vint(32)<=vint(12)) mint(51) = 1
  Else If (ilim==3) Then
    yst = vint(22)
    ctnmn0 = -1.
    ctnmx0 = 0.
    ctpmn0 = 0.
    ctpmx0 = 1.
    ctnmn1 = min(0., ckin(27))
    ctnmx1 = min(0., ckin(28))
    ctpmn1 = max(0., ckin(27))
    ctpmx1 = max(0., ckin(28))
    ctnmn2 = -sqrt(1.-4.*pthmin**2/(be34**2*tau*vint(2)))
    ctpmx2 = -ctnmn2
    ctnmx2 = 0.
    ctpmn2 = 0.
    If (ckin(4)>=0.) Then
      ctnmx2 = -sqrt(max(0.,1.-4.*ckin(4)**2/(be34**2*tau*vint(2))))
      ctpmn2 = -ctnmx2
    End If
    ctnmn3 = min(0., max((1.+rm3-rm4)/be34*tanh(ckin(11)-yst),-(1.-rm3+rm4)/be34*tanh(ckin(10)-yst)))
    ctnmx3 = min(0., (1.+rm3-rm4)/be34*tanh(ckin(12)-yst), -(1.-rm3+rm4)/be34*tanh(ckin(9)-yst))
    ctpmn3 = max(0., (1.+rm3-rm4)/be34*tanh(ckin(9)-yst), -(1.-rm3+rm4)/be34*tanh(ckin(12)-yst))
    ctpmx3 = max(0., min((1.+rm3-rm4)/be34*tanh(ckin(10)-yst),-(1.-rm3+rm4)/be34*tanh(ckin(11)-yst)))
    vint(13) = max(ctnmn0, ctnmn1, ctnmn2, ctnmn3)
    vint(33) = min(ctnmx0, ctnmx1, ctnmx2, ctnmx3)
    vint(14) = max(ctpmn0, ctpmn1, ctpmn2, ctpmn3)
    vint(34) = min(ctpmx0, ctpmx1, ctpmx2, ctpmx3)
    If (vint(33)<=vint(13) .And. vint(34)<=vint(14)) mint(51) = 1
  Else If (ilim==4) Then
    tapmn0 = tau
    tapmx0 = 1.
    tapmn1 = ckin(31)**2/vint(2)
    tapmx1 = 1.
    If (ckin(32)>=0.) tapmx1 = ckin(32)**2/vint(2)
    vint(16) = max(tapmn0, tapmn1)
    vint(36) = min(tapmx0, tapmx1)
    If (mint(43)==1) Then
      vint(16) = 0.99999
      vint(36) = 1.00001
    End If
    If (vint(36)<=vint(16)) mint(51) = 1
  End If
  Return
  110 If (ilim==0) Then
  Else If (ilim==1) Then
    If (mstp(82)<=1) vint(11) = 4.*parp(81)**2/vint(2)
    If (mstp(82)>=2) vint(11) = parp(82)**2/vint(2)
    vint(31) = 1.
  Else If (ilim==2) Then
    vint(12) = 0.5*log(vint(21))
    vint(32) = -vint(12)
  Else If (ilim==3) Then
    If (mstp(82)<=1) st2eff = 4.*parp(81)**2/(vint(21)*vint(2))
    If (mstp(82)>=2) st2eff = 0.01*parp(82)**2/(vint(21)*vint(2))
    vint(13) = -sqrt(max(0.,1.-st2eff))
    vint(33) = 0.
    vint(14) = 0.
    vint(34) = -vint(13)
  End If
  Return
End Subroutine pyklim
