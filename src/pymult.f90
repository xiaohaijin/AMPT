Subroutine pymult(mmul)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension nmul(20), sigm(20), kstr(500, 2)
  Save xt2, xt2fac, xc2, xts, irbin, rbin, nmul, sigm
  If (mmul==1) Then
    If (mstp(122)>=1) Write (mstu(11), 1000) mstp(82)
    isub = 96
    mint(1) = 96
    vint(63) = 0.
    vint(64) = 0.
    vint(143) = 1.
    vint(144) = 1.
    100 sigsum = 0.
    Do ixt2 = 1, 20
      nmul(ixt2) = mstp(83)
      sigm(ixt2) = 0.
      Do itry = 1, mstp(83)
        rsca = 0.05*((21-ixt2)-rlu(0))
        xt2 = vint(149)*(1.+vint(149))/(vint(149)+rsca) - vint(149)
        xt2 = max(0.01*vint(149), xt2)
        vint(25) = xt2
        If (rlu(0)<=coef(isub,1)) Then
          taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
          tau = xt2*(1.+taup)**2/(4.*taup)
        Else
          tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
        End If
        vint(21) = tau
        Call pyklim(2)
        ryst = rlu(0)
        myst = 1
        If (ryst>coef(isub,7)) myst = 2
        If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
        Call pykmap(2, myst, rlu(0))
        vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))
        vint(71) = 0.5*vint(1)*sqrt(xt2)
        Call pysigh(nchn, sigs)
        sigm(ixt2) = sigm(ixt2) + sigs
      End Do
      sigsum = sigsum + sigm(ixt2)
    End Do
    sigsum = sigsum/(20.*mstp(83))
    If (sigsum<1.1*vint(106)) Then
      If (mstp(122)>=1) Write (mstu(11), 1100) parp(82), sigsum
      parp(82) = 0.9*parp(82)
      vint(149) = 4.*parp(82)**2/vint(2)
      Goto 100
    End If
    If (mstp(122)>=1) Write (mstu(11), 1200) parp(82), sigsum
    yke = sigsum/vint(106)
    so = 0.5
    xi = 0.
    yi = 0.
    xk = 0.5
    iit = 0
    130 If (iit==0) Then
      xk = 2.*xk
    Else If (iit==1) Then
      xk = 0.5*xk
    Else
      xk = xi + (yke-yi)*(xf-xi)/(yf-yi)
    End If
    If (mstp(82)==2) Then
      sp = 0.5*paru(1)*(1.-exp(-xk))
      sop = sp/paru(1)
    Else
      If (mstp(82)==3) deltab = 0.02
      If (mstp(82)==4) deltab = min(0.01, 0.05*parp(84))
      sp = 0.
      sop = 0.
      b = -0.5*deltab
      140 b = b + deltab
      If (mstp(82)==3) Then
        ov = exp(-b**2)/paru(2)
      Else
        cq2 = parp(84)**2
        ov = ((1.-parp(83))**2*exp(-min(100.,b**2))+2.*parp(83)*(1.-parp(83))*2./(1.+cq2)*exp(-min(100.,b**2*2./(1.+cq2)))+parp(83)**2/cq2*exp(-min(100.,b**2/cq2)))/paru(2)
      End If
      pacc = 1. - exp(-min(100.,paru(1)*xk*ov))
      sp = sp + paru(2)*b*deltab*pacc
      sop = sop + paru(2)*b*deltab*ov*pacc
      If (b<1. .Or. b*pacc>1E-6) Goto 140
    End If
    yk = paru(1)*xk*so/sp
    If (yk<yke) Then
      xi = xk
      yi = yk
      If (iit==1) iit = 2
    Else
      xf = xk
      yf = yk
      If (iit==0) iit = 1
    End If
    If (abs(yk-yke)>=1E-5*yke) Goto 130
    vint(145) = sigsum
    vint(146) = sop/so
    vint(147) = sop/sp
  Else If (mmul==2) Then
    If (mstp(82)<=0) Then
    Else If (mstp(82)==1) Then
      xt2 = 1.
      xt2fac = xsec(96, 1)/vint(106)*vint(149)/(1.-vint(149))
    Else If (mstp(82)==2) Then
      xt2 = 1.
      xt2fac = vint(146)*xsec(96, 1)/vint(106)*vint(149)*(1.+vint(149))
    Else
      xc2 = 4.*ckin(3)**2/vint(2)
      If (ckin(3)<=ckin(5) .Or. mint(82)>=2) xc2 = 0.
    End If
  Else If (mmul==3) Then
    isub = mint(1)
    If (mstp(82)<=0) Then
      xt2 = 0.
    Else If (mstp(82)==1) Then
      xt2 = xt2fac*xt2/(xt2fac-xt2*log(rlu(0)))
    Else If (mstp(82)==2) Then
      If (xt2<1. .And. exp(-xt2fac*xt2/(vint(149)*(xt2+vint(149))))>rlu(0)) xt2 = 1.
      If (xt2>=1.) Then
        xt2 = (1.+vint(149))*xt2fac/(xt2fac-(1.+vint(149))*log(1.-rlu(0)*(1.-exp(-xt2fac/(vint(149)*(1.+vint(149))))))) - vint(149)
      Else
        xt2 = -xt2fac/log(exp(-xt2fac/(xt2+vint(149)))+rlu(0)*(exp(-xt2fac/vint(149))-exp(-xt2fac/(xt2+vint(149))))) - vint(149)
      End If
      xt2 = max(0.01*vint(149), xt2)
    Else
      xt2 = (xc2+vint(149))*(1.+vint(149))/(1.+vint(149)-rlu(0)*(1.-xc2)) - vint(149)
      xt2 = max(0.01*vint(149), xt2)
    End If
    vint(25) = xt2
    If (mstp(82)<=1 .And. xt2<vint(149)) Then
      If (mint(82)==1) ngen(0, 1) = ngen(0, 1) - 1
      If (mint(82)==1) ngen(isub, 1) = ngen(isub, 1) - 1
      isub = 95
      mint(1) = isub
      vint(21) = 0.01*vint(149)
      vint(22) = 0.
      vint(23) = 0.
      vint(25) = 0.01*vint(149)
    Else
      If (rlu(0)<=coef(isub,1)) Then
        taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
        tau = xt2*(1.+taup)**2/(4.*taup)
      Else
        tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
      End If
      vint(21) = tau
      Call pyklim(2)
      ryst = rlu(0)
      myst = 1
      If (ryst>coef(isub,7)) myst = 2
      If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
      Call pykmap(2, myst, rlu(0))
      vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))
    End If
    vint(71) = 0.5*vint(1)*sqrt(vint(25))
  Else If (mmul==4) Then
    isub = mint(1)
    xts = vint(25)
    If (iset(isub)==1) xts = vint(21)
    If (iset(isub)==2) xts = (4.*vint(48)+2.*vint(63)+2.*vint(64))/vint(2)
    If (iset(isub)==3 .Or. iset(isub)==4) xts = vint(26)
    rbin = max(0.000001, min(0.999999,xts*(1.+vint(149))/(xts+vint(149))))
    irbin = int(1.+20.*rbin)
    If (isub==96) nmul(irbin) = nmul(irbin) + 1
    If (isub==96) sigm(irbin) = sigm(irbin) + vint(153)
  Else If (mmul==5) Then
    If (mstp(82)==3) Then
      vint(148) = rlu(0)/(paru(2)*vint(147))
    Else
      rtype = rlu(0)
      cq2 = parp(84)**2
      If (rtype<(1.-parp(83))**2) Then
        b2 = -log(rlu(0))
      Else If (rtype<1.-parp(83)**2) Then
        b2 = -0.5*(1.+cq2)*log(rlu(0))
      Else
        b2 = -cq2*log(rlu(0))
      End If
      vint(148) = ((1.-parp(83))**2*exp(-min(100.,b2))+2.*parp(83)*(1.-parp(83))*2./(1.+cq2)*exp(-min(100.,b2*2./(1.+cq2)))+parp(83)**2/cq2*exp(-min(100.,b2/cq2)))/(paru(2)*vint(147))
    End If
    rncor = (irbin-20.*rbin)*nmul(irbin)
    sigcor = (irbin-20.*rbin)*sigm(irbin)
    Do ibin = irbin + 1, 20
      rncor = rncor + nmul(ibin)
      sigcor = sigcor + sigm(ibin)
    End Do
    sigabv = (sigcor/rncor)*vint(149)*(1.-xts)/(xts+vint(149))
    vint(150) = exp(-min(100.,vint(146)*vint(148)*sigabv/vint(106)))
  Else If (mmul==6) Then
    isub = mint(1)
    nmax = mint(84) + 4
    If (iset(isub)==1) nmax = mint(84) + 2
    nstr = 0
    Do i = mint(84) + 1, nmax
      kcs = kchg(lucomp(k(i,2)), 2)*isign(1, k(i,2))
      If (kcs==0) Goto 170
      Do j = 1, 4
        If (kcs==1 .And. (j==2 .Or. j==4)) Goto 160
        If (kcs==-1 .And. (j==1 .Or. j==3)) Goto 160
        If (j<=2) Then
          ist = mod(k(i,j+3)/mstu(5), mstu(5))
        Else
          ist = mod(k(i,j+1), mstu(5))
        End If
        If (ist<mint(84) .Or. ist>i) Goto 160
        If (kchg(lucomp(k(ist,2)),2)==0) Goto 160
        nstr = nstr + 1
        If (j==1 .Or. j==4) Then
          kstr(nstr, 1) = i
          kstr(nstr, 2) = ist
        Else
          kstr(nstr, 1) = ist
          kstr(nstr, 2) = i
        End If
      160 End Do
    170 End Do
    xt2 = vint(25)
    If (iset(isub)==1) xt2 = vint(21)
    If (iset(isub)==2) xt2 = (4.*vint(48)+2.*vint(63)+2.*vint(64))/vint(2)
    If (iset(isub)==3 .Or. iset(isub)==4) xt2 = vint(26)
    isub = 96
    mint(1) = 96
    If (mstp(82)<=1) Then
      xt2fac = xsec(isub, 1)*vint(149)/((1.-vint(149))*vint(106))
    Else
      xt2fac = vint(146)*vint(148)*xsec(isub, 1)/vint(106)*vint(149)*(1.+vint(149))
    End If
    vint(63) = 0.
    vint(64) = 0.
    vint(151) = 0.
    vint(152) = 0.
    vint(143) = 1. - vint(141)
    vint(144) = 1. - vint(142)
    180 If (mstp(82)<=1) Then
      xt2 = xt2fac*xt2/(xt2fac-xt2*log(rlu(0)))
      If (xt2<vint(149)) Goto 220
    Else
      If (xt2<=0.01*vint(149)) Goto 220
      xt2 = xt2fac*(xt2+vint(149))/(xt2fac-(xt2+vint(149))*log(rlu(0))) - vint(149)
      If (xt2<=0.) Goto 220
      xt2 = max(0.01*vint(149), xt2)
    End If
    vint(25) = xt2
    If (rlu(0)<=coef(isub,1)) Then
      taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
      tau = xt2*(1.+taup)**2/(4.*taup)
    Else
      tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
    End If
    vint(21) = tau
    Call pyklim(2)
    ryst = rlu(0)
    myst = 1
    If (ryst>coef(isub,7)) myst = 2
    If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
    Call pykmap(2, myst, rlu(0))
    vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))
    x1m = sqrt(tau)*exp(vint(22))
    x2m = sqrt(tau)*exp(-vint(22))
    If (vint(143)-x1m<0.01 .Or. vint(144)-x2m<0.01) Goto 180
    vint(71) = 0.5*vint(1)*sqrt(xt2)
    Call pysigh(nchn, sigs)
    If (sigs<xsec(isub,1)*rlu(0)) Goto 180
    Do i = n + 1, n + 2
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
    End Do
    rflav = rlu(0)
    pt = 0.5*vint(1)*sqrt(xt2)
    phi = paru(2)*rlu(0)
    cth = vint(23)
    k(n+1, 1) = 3
    k(n+1, 2) = 21
    If (rflav>=max(parp(85),parp(86))) k(n+1, 2) = 1 + int((2.+parj(2))*rlu(0))
    p(n+1, 1) = pt*cos(phi)
    p(n+1, 2) = pt*sin(phi)
    p(n+1, 3) = 0.25*vint(1)*(vint(41)*(1.+cth)-vint(42)*(1.-cth))
    p(n+1, 4) = 0.25*vint(1)*(vint(41)*(1.+cth)+vint(42)*(1.-cth))
    p(n+1, 5) = 0.
    k(n+2, 1) = 3
    k(n+2, 2) = 21
    If (k(n+1,2)/=21) k(n+2, 2) = -k(n+1, 2)
    p(n+2, 1) = -p(n+1, 1)
    p(n+2, 2) = -p(n+1, 2)
    p(n+2, 3) = 0.25*vint(1)*(vint(41)*(1.-cth)-vint(42)*(1.+cth))
    p(n+2, 4) = 0.25*vint(1)*(vint(41)*(1.-cth)+vint(42)*(1.+cth))
    p(n+2, 5) = 0.
    If (rflav<parp(85) .And. nstr>=1) Then
      Do i = n + 1, n + 2
        dmin = 1E8
        Do istr = 1, nstr
          i1 = kstr(istr, 1)
          i2 = kstr(istr, 2)
          dist = (p(i,4)*p(i1,4)-p(i,1)*p(i1,1)-p(i,2)*p(i1,2)-p(i,3)*p(i1,3))*(p(i,4)*p(i2,4)-p(i,1)*p(i2,1)-p(i,2)*p(i2,2)-p(i,3)*p(i2,3))/max(1., p(i1,4)*p(i2,4)-p(i1,1)*p(i2,1)-p(i1,2)*p(i2,2)-p(i1,3)*p(i2,3))
          If (istr==1 .Or. dist<dmin) Then
            dmin = dist
            ist1 = i1
            ist2 = i2
            istm = istr
          End If
        End Do
        If (k(ist1,4)/mstu(5)==ist2) k(ist1, 4) = mstu(5)*i + mod(k(ist1,4), mstu(5))
        If (mod(k(ist1,5),mstu(5))==ist2) k(ist1, 5) = mstu(5)*(k(ist1,5)/mstu(5)) + i
        k(i, 5) = mstu(5)*ist1
        k(i, 4) = mstu(5)*ist2
        If (k(ist2,5)/mstu(5)==ist1) k(ist2, 5) = mstu(5)*i + mod(k(ist2,5), mstu(5))
        If (mod(k(ist2,4),mstu(5))==ist1) k(ist2, 4) = mstu(5)*(k(ist2,4)/mstu(5)) + i
        kstr(istm, 2) = i
        kstr(nstr+1, 1) = i
        kstr(nstr+1, 2) = ist2
        nstr = nstr + 1
      End Do
    Else If (k(n+1,2)==21) Then
      k(n+1, 4) = mstu(5)*(n+2)
      k(n+1, 5) = mstu(5)*(n+2)
      k(n+2, 4) = mstu(5)*(n+1)
      k(n+2, 5) = mstu(5)*(n+1)
      kstr(nstr+1, 1) = n + 1
      kstr(nstr+1, 2) = n + 2
      kstr(nstr+2, 1) = n + 2
      kstr(nstr+2, 2) = n + 1
      nstr = nstr + 2
    Else
      k(n+1, 4) = mstu(5)*(n+2)
      k(n+2, 5) = mstu(5)*(n+1)
      kstr(nstr+1, 1) = n + 1
      kstr(nstr+1, 2) = n + 2
      nstr = nstr + 1
    End If
    n = n + 2
    If (n>mstu(4)-mstu(32)-10) Then
      Call luerrm(11, '(PYMULT:) no more memory left in LUJETS')
      If (mstu(21)>=1) Return
    End If
    mint(31) = mint(31) + 1
    vint(151) = vint(151) + vint(41)
    vint(152) = vint(152) + vint(42)
    vint(143) = vint(143) - vint(41)
    vint(144) = vint(144) - vint(42)
    If (mint(31)<240) Goto 180
    220 Continue
  End If
  Return
  1000 Format (/1X, '****** PYMULT: initialization of multiple inter', 'actions for MSTP(82) =', I2, ' ******')
  1100 Format (8X, 'pT0 =', F5.2, ' GeV gives sigma(parton-parton) =', 1P, E9.2, ' mb: rejected')
  1200 Format (8X, 'pT0 =', F5.2, ' GeV gives sigma(parton-parton) =', 1P, E9.2, ' mb: accepted')
End Subroutine pymult
