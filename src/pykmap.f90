Subroutine pykmap(ivar, mvar, vvar)
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  isub = mint(1)
  If (ivar==1) Then
    taumin = vint(11)
    taumax = vint(31)
    If (mvar==3 .Or. mvar==4) Then
      taure = vint(73)
      gamre = vint(74)
    Else If (mvar==5 .Or. mvar==6) Then
      taure = vint(75)
      gamre = vint(76)
    End If
    If (mint(43)==1 .And. (iset(isub)==1 .Or. iset(isub)==2)) Then
      tau = 1.
    Else If (mvar==1) Then
      tau = taumin*(taumax/taumin)**vvar
    Else If (mvar==2) Then
      tau = taumax*taumin/(taumin+(taumax-taumin)*vvar)
    Else If (mvar==3 .Or. mvar==5) Then
      ratgen = (taure+taumax)/(taure+taumin)*taumin/taumax
      tau = taure*taumin/((taure+taumin)*ratgen**vvar-taumin)
    Else
      aupp = atan((taumax-taure)/gamre)
      alow = atan((taumin-taure)/gamre)
      tau = taure + gamre*tan(alow+(aupp-alow)*vvar)
    End If
    vint(21) = min(taumax, max(taumin,tau))
  Else If (ivar==2) Then
    ystmin = vint(12)
    ystmax = vint(32)
    If (mint(43)==1) Then
      yst = 0.
    Else If (mint(43)==2) Then
      If (iset(isub)<=2) yst = -0.5*log(vint(21))
      If (iset(isub)>=3) yst = -0.5*log(vint(26))
    Else If (mint(43)==3) Then
      If (iset(isub)<=2) yst = 0.5*log(vint(21))
      If (iset(isub)>=3) yst = 0.5*log(vint(26))
    Else If (mvar==1) Then
      yst = ystmin + (ystmax-ystmin)*sqrt(vvar)
    Else If (mvar==2) Then
      yst = ystmax - (ystmax-ystmin)*sqrt(1.-vvar)
    Else
      aupp = atan(exp(ystmax))
      alow = atan(exp(ystmin))
      yst = log(tan(alow+(aupp-alow)*vvar))
    End If
    vint(22) = min(ystmax, max(ystmin,yst))
  Else If (ivar==3) Then
    rm34 = 2.*vint(63)*vint(64)/(vint(21)*vint(2))**2
    rsqm = 1. + rm34
    If (2.*vint(71)**2/(vint(21)*vint(2))<0.0001) rm34 = max(rm34, 2.*vint(71)**2/(vint(21)*vint(2)))
    ctnmin = vint(13)
    ctnmax = vint(33)
    ctpmin = vint(14)
    ctpmax = vint(34)
    If (mvar==1) Then
      aneg = ctnmax - ctnmin
      apos = ctpmax - ctpmin
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = ctnmin + (ctnmax-ctnmin)*vctn
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = ctpmin + (ctpmax-ctpmin)*vctp
      End If
    Else If (mvar==2) Then
      rmnmin = max(rm34, rsqm-ctnmin)
      rmnmax = max(rm34, rsqm-ctnmax)
      rmpmin = max(rm34, rsqm-ctpmin)
      rmpmax = max(rm34, rsqm-ctpmax)
      aneg = log(rmnmin/rmnmax)
      apos = log(rmpmin/rmpmax)
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rsqm - rmnmin*(rmnmax/rmnmin)**vctn
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rsqm - rmpmin*(rmpmax/rmpmin)**vctp
      End If
    Else If (mvar==3) Then
      rmnmin = max(rm34, rsqm+ctnmin)
      rmnmax = max(rm34, rsqm+ctnmax)
      rmpmin = max(rm34, rsqm+ctpmin)
      rmpmax = max(rm34, rsqm+ctpmax)
      aneg = log(rmnmax/rmnmin)
      apos = log(rmpmax/rmpmin)
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rmnmin*(rmnmax/rmnmin)**vctn - rsqm
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rmpmin*(rmpmax/rmpmin)**vctp - rsqm
      End If
    Else If (mvar==4) Then
      rmnmin = max(rm34, rsqm-ctnmin)
      rmnmax = max(rm34, rsqm-ctnmax)
      rmpmin = max(rm34, rsqm-ctpmin)
      rmpmax = max(rm34, rsqm-ctpmax)
      aneg = 1./rmnmax - 1./rmnmin
      apos = 1./rmpmax - 1./rmpmin
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rsqm - 1./(1./rmnmin+aneg*vctn)
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rsqm - 1./(1./rmpmin+apos*vctp)
      End If
    Else If (mvar==5) Then
      rmnmin = max(rm34, rsqm+ctnmin)
      rmnmax = max(rm34, rsqm+ctnmax)
      rmpmin = max(rm34, rsqm+ctpmin)
      rmpmax = max(rm34, rsqm+ctpmax)
      aneg = 1./rmnmin - 1./rmnmax
      apos = 1./rmpmin - 1./rmpmax
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = 1./(1./rmnmin-aneg*vctn) - rsqm
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = 1./(1./rmpmin-apos*vctp) - rsqm
      End If
    End If
    If (cth<0.) cth = min(ctnmax, max(ctnmin,cth))
    If (cth>0.) cth = min(ctpmax, max(ctpmin,cth))
    vint(23) = cth
  Else If (ivar==4) Then
    tau = vint(11)
    taupmn = vint(16)
    taupmx = vint(36)
    If (mint(43)==1) Then
      taup = 1.
    Else If (mvar==1) Then
      taup = taupmn*(taupmx/taupmn)**vvar
    Else
      aupp = (1.-tau/taupmx)**4
      alow = (1.-tau/taupmn)**4
      taup = tau/(1.-(alow+(aupp-alow)*vvar)**0.25)
    End If
    vint(26) = min(taupmx, max(taupmn,taup))
  End If
  Return
End Subroutine pykmap
