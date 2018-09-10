Subroutine pywidt(kflr, rmas, wdtp, wdte)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Dimension wdtp(0:40), wdte(0:40, 0:5)
  kfla = iabs(kflr)
  sqm = rmas**2
  as = ulalps(sqm)
  aem = paru(101)
  xw = paru(102)
  radc = 1. + as/paru(1)
  Do i = 0, 40
    wdtp(i) = 0.
    Do j = 0, 5
      wdte(i, j) = 0.
    End Do
  End Do
  If (kfla==21) Then
    Do i = 1, mdcy(21, 3)
      idc = i + mdcy(21, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 110
      If (i<=8) Then
        wdtp(i) = (1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    110 End Do
  Else If (kfla==23) Then
    If (mint(61)==1) Then
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei)
      vi = ai - 4.*ei*xw
      sqmz = pmas(23, 1)**2
      gzmz = pmas(23, 2)*pmas(23, 1)
      ggi = ei**2
      gzi = ei*vi/(8.*xw*(1.-xw))*sqm*(sqm-sqmz)/((sqm-sqmz)**2+gzmz**2)
      zzi = (vi**2+ai**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmz)**2+gzmz**2)
      If (mstp(43)==1) Then
        gzi = 0.
        zzi = 0.
      Else If (mstp(43)==2) Then
        ggi = 0.
        gzi = 0.
      End If
    Else If (mint(61)==2) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(114) = 0.
    End If
    Do i = 1, mdcy(23, 3)
      idc = i + mdcy(23, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 120
      If (i<=8) Then
        ef = kchg(i, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        If (mint(61)==0) Then
          wdtp(i) = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==1) Then
          wdtp(i) = 3.*((ggi*ef**2+gzi*ef*vf+zzi*vf**2)*(1.+2.*rm1)+zzi*af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==2) Then
          ggf = 3.*ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzf = 3.*ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          zzf = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        End If
        wid2 = 1.
      Else If (i<=16) Then
        ef = kchg(i+2, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        wdtp(i) = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        If (mint(61)==0) Then
          wdtp(i) = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = ((ggi*ef**2+gzi*ef*vf+zzi*vf**2)*(1.+2.*rm1)+zzi*af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = 1.
      Else
        cf = 2.*(1.-2.*xw)
        If (mint(61)==0) Then
          wdtp(i) = 0.25*cf**2*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = 0.25*(ggi+gzi*cf+zzi*cf**2)*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = 0.25*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = 0.25*cf*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = 0.25*cf**2*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = wids(37, 1)
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
        If (mint(61)==2) Then
          vint(111) = vint(111) + ggf*wid2
          vint(112) = vint(112) + gzf*wid2
          vint(114) = vint(114) + zzf*wid2
        End If
      End If
    120 End Do
    If (mstp(43)==1) Then
      vint(112) = 0.
      vint(114) = 0.
    Else If (mstp(43)==2) Then
      vint(111) = 0.
      vint(112) = 0.
    End If
  Else If (kfla==24) Then
    Do i = 1, mdcy(24, 3)
      idc = i + mdcy(24, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 130
      If (i<=16) Then
        wdtp(i) = 3.*(2.-rm1-rm2-(rm1-rm2)**2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))*vckm((i-1)/4+1, mod(i-1,4)+1)*radc
        wid2 = 1.
      Else
        wdtp(i) = (2.-rm1-rm2-(rm1-rm2)**2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    130 End Do
  Else If (kfla==25) Then
    Do i = 1, mdcy(25, 3)
      idc = i + mdcy(25, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 170
      If (i<=8) Then
        wdtp(i) = 3.*rm1*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
        wid2 = 1.
      Else If (i<=12) Then
        wdtp(i) = rm1*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        wid2 = 1.
      Else If (i==13) Then
        etare = 0.
        etaim = 0.
        Do j = 1, 2*mstp(1)
          eps = (2.*pmas(j,1)/rmas)**2
          If (eps<=1.) Then
            If (eps>1.E-4) Then
              root = sqrt(1.-eps)
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
          End If
          etare = etare + 0.5*eps*(1.+(eps-1.)*phire)
          etaim = etaim + 0.5*eps*(eps-1.)*phiim
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (as/paru(1))**2*eta2
        wid2 = 1.
      Else If (i==14) Then
        etare = 0.
        etaim = 0.
        Do j = 1, 3*mstp(1) + 1
          If (j<=2*mstp(1)) Then
            ej = kchg(j, 1)/3.
            eps = (2.*pmas(j,1)/rmas)**2
          Else If (j<=3*mstp(1)) Then
            jl = 2*(j-2*mstp(1)) - 1
            ej = kchg(10+jl, 1)/3.
            eps = (2.*pmas(10+jl,1)/rmas)**2
          Else
            eps = (2.*pmas(24,1)/rmas)**2
          End If
          If (eps<=1.) Then
            If (eps>1.E-4) Then
              root = sqrt(1.-eps)
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
          End If
          If (j<=2*mstp(1)) Then
            etare = etare + 0.5*3.*ej**2*eps*(1.+(eps-1.)*phire)
            etaim = etaim + 0.5*3.*ej**2*eps*(eps-1.)*phiim
          Else If (j<=3*mstp(1)) Then
            etare = etare + 0.5*ej**2*eps*(1.+(eps-1.)*phire)
            etaim = etaim + 0.5*ej**2*eps*(eps-1.)*phiim
          Else
            etare = etare - 0.5 - 0.75*eps*(1.+(eps-2.)*phire)
            etaim = etaim + 0.75*eps*(eps-2.)*phiim
          End If
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (aem/paru(1))**2*0.5*eta2
        wid2 = 1.
      Else If (i==15) Then
        etare = 0.
        etaim = 0.
        Do j = 1, 3*mstp(1) + 1
          If (j<=2*mstp(1)) Then
            ej = kchg(j, 1)/3.
            aj = sign(1., ej+0.1)
            vj = aj - 4.*ej*xw
            eps = (2.*pmas(j,1)/rmas)**2
            epsp = (2.*pmas(j,1)/pmas(23,1))**2
          Else If (j<=3*mstp(1)) Then
            jl = 2*(j-2*mstp(1)) - 1
            ej = kchg(10+jl, 1)/3.
            aj = sign(1., ej+0.1)
            vj = ai - 4.*ej*xw
            eps = (2.*pmas(10+jl,1)/rmas)**2
            epsp = (2.*pmas(10+jl,1)/pmas(23,1))**2
          Else
            eps = (2.*pmas(24,1)/rmas)**2
            epsp = (2.*pmas(24,1)/pmas(23,1))**2
          End If
          If (eps<=1.) Then
            root = sqrt(1.-eps)
            If (eps>1.E-4) Then
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
            psire = -(1.+0.5*root*rln)
            psiim = 0.5*paru(1)*root
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
            psire = -(1.+sqrt(eps-1.)*asin(1./sqrt(eps)))
            psiim = 0.
          End If
          If (epsp<=1.) Then
            root = sqrt(1.-epsp)
            If (epsp>1.E-4) Then
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./epsp-2.)
            End If
            phirep = 0.25*(rln**2-paru(1)**2)
            phiimp = 0.5*paru(1)*rln
            psirep = -(1.+0.5*root*rln)
            psiimp = 0.5*paru(1)*root
          Else
            phirep = -(asin(1./sqrt(epsp)))**2
            phiimp = 0.
            psirep = -(1.+sqrt(epsp-1.)*asin(1./sqrt(epsp)))
            psiimp = 0.
          End If
          fxyre = eps*epsp/(8.*(eps-epsp))*(1.-eps*epsp/(eps-epsp)*(phire-phirep)+2.*eps/(eps-epsp)*(psire-psirep))
          fxyim = eps*epsp/(8.*(eps-epsp))*(-eps*epsp/(eps-epsp)*(phiim-phiimp)+2.*eps/(eps-epsp)*(psiim-psiimp))
          f1re = eps*epsp/(2.*(eps-epsp))*(phire-phirep)
          f1im = eps*epsp/(2.*(eps-epsp))*(phiim-phiimp)
          If (j<=2*mstp(1)) Then
            etare = etare - 3.*ej*vj*(fxyre-0.25*f1re)
            etaim = etaim - 3.*ej*vj*(fxyim-0.25*f1im)
          Else If (j<=3*mstp(1)) Then
            etare = etare - ej*vj*(fxyre-0.25*f1re)
            etaim = etaim - ej*vj*(fxyim-0.25*f1im)
          Else
            etare = etare - sqrt(1.-xw)*(((1.+2./eps)*xw/sqrt(1.-xw)-(5.+2./eps))*fxyre+(3.-xw/sqrt(1.-xw))*f1re)
            etaim = etaim - sqrt(1.-xw)*(((1.+2./eps)*xw/sqrt(1.-xw)-(5.+2./eps))*fxyim+(3.-xw/sqrt(1.-xw))*f1im)
          End If
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (aem/paru(1))**2*(1.-(pmas(23,1)/rmas)**2)**3/xw*eta2
        wid2 = wids(23, 2)
      Else
        wdtp(i) = (1.-4.*rm1+12.*rm1**2)*sqrt(max(0.,1.-4.*rm1))/(2.*(18-i))
        wid2 = wids(7+i, 1)
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    170 End Do
  Else If (kfla==32) Then
    If (mint(61)==1) Then
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei)
      vi = ai - 4.*ei*xw
      sqmz = pmas(23, 1)**2
      gzmz = pmas(23, 2)*pmas(23, 1)
      api = sign(1., ei)
      vpi = api - 4.*ei*xw
      sqmzp = pmas(32, 1)**2
      gzpmzp = pmas(32, 2)*pmas(32, 1)
      ggi = ei**2
      gzi = ei*vi/(8.*xw*(1.-xw))*sqm*(sqm-sqmz)/((sqm-sqmz)**2+gzmz**2)
      gzpi = ei*vpi/(8.*xw*(1.-xw))*sqm*(sqm-sqmzp)/((sqm-sqmzp)**2+gzpmzp**2)
      zzi = (vi**2+ai**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmz)**2+gzmz**2)
      zzpi = 2.*(vi*vpi+ai*api)/(16.*xw*(1.-xw))**2*sqm**2*((sqm-sqmz)*(sqm-sqmzp)+gzmz*gzpmzp)/(((sqm-sqmz)**2+gzmz**2)*((sqm-sqmzp)**2+gzpmzp**2))
      zpzpi = (vpi**2+api**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmzp)**2+gzpmzp**2)
      If (mstp(44)==1) Then
        gzi = 0.
        gzpi = 0.
        zzi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==2) Then
        ggi = 0.
        gzi = 0.
        gzpi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==3) Then
        ggi = 0.
        gzi = 0.
        gzpi = 0.
        zzi = 0.
        zzpi = 0.
      Else If (mstp(44)==4) Then
        gzpi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==5) Then
        gzi = 0.
        zzi = 0.
        zzpi = 0.
      Else If (mstp(44)==6) Then
        ggi = 0.
        gzi = 0.
        gzpi = 0.
      End If
    Else If (mint(61)==2) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
      vint(116) = 0.
    End If
    Do i = 1, mdcy(32, 3)
      idc = i + mdcy(32, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 180
      If (i<=8) Then
        ef = kchg(i, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        apf = sign(1., ef+0.1)
        vpf = apf - 4.*ef*xw
        If (mint(61)==0) Then
          wdtp(i) = 3.*(vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==1) Then
          wdtp(i) = 3.*((ggi*ef**2+gzi*ef*vf+gzpi*ef*vpf+zzi*vf**2+zzpi*vf*vpf+zpzpi*vpf**2)*(1.+2.*rm1)+(zzi*af**2+zzpi*af*apf+zpzpi*apf**2)*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==2) Then
          ggf = 3.*ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzf = 3.*ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzpf = 3.*ef*vpf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          zzf = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
          zzpf = 3.*(vf*vpf*(1.+2.*rm1)+af*apf*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
          zpzpf = 3.*(vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        End If
        wid2 = 1.
      Else
        ef = kchg(i+2, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        If (i<=10) Then
          vpf = paru(127-2*mod(i,2))
          apf = paru(128-2*mod(i,2))
        Else If (i<=12) Then
          vpf = parj(186-2*mod(i,2))
          apf = parj(187-2*mod(i,2))
        Else
          vpf = parj(194-2*mod(i,2))
          apf = parj(195-2*mod(i,2))
        End If
        If (mint(61)==0) Then
          wdtp(i) = (vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = ((ggi*ef**2+gzi*ef*vf+gzpi*ef*vpf+zzi*vf**2+zzpi*vf*vpf+zpzpi*vpf**2)*(1.+2.*rm1)+(zzi*af**2+zzpi*af*apf+zpzpi*apf**2)*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzpf = ef*vpf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
          zzpf = (vf*vpf*(1.+2.*rm1)+af*apf*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
          zpzpf = (vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
        If (mint(61)==2) Then
          vint(111) = vint(111) + ggf
          vint(112) = vint(112) + gzf
          vint(113) = vint(113) + gzpf
          vint(114) = vint(114) + zzf
          vint(115) = vint(115) + zzpf
          vint(116) = vint(116) + zpzpf
        End If
      End If
    180 End Do
    If (mstp(44)==1) Then
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==2) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==3) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
    Else If (mstp(44)==4) Then
      vint(113) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==5) Then
      vint(112) = 0.
      vint(114) = 0.
      vint(115) = 0.
    Else If (mstp(44)==6) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
    End If
  Else If (kfla==37) Then
    Do i = 1, mdcy(37, 3)
      idc = i + mdcy(37, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 190
      If (i<=4) Then
        wdtp(i) = 3.*((rm1*paru(121)+rm2/paru(121))*(1.-rm1-rm2)-4.*rm1*rm2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))*radc
        wid2 = 1.
      Else
        wdtp(i) = ((rm1*paru(121)+rm2/paru(121))*(1.-rm1-rm2)-4.*rm1*rm2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    190 End Do
  Else If (kfla==40) Then
    Do i = 1, mdcy(40, 3)
      idc = i + mdcy(40, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 200
      If (i<=4) Then
        wdtp(i) = 3.*radc
        wid2 = 1.
      Else
        wdtp(i) = 1.
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    200 End Do
  End If
  mint(61) = 0
  Return
End Subroutine pywidt
