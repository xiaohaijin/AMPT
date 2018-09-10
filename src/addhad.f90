Subroutine addhad
  Parameter (maxstr=150001, maxr=1, xmd=1.8756)
  Double Precision smearp, smearh
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /smearz/smearp, smearh
  Common /rndf77/nseed
  Common /para8/idpert, npertd, idxsec
  Save
  np0 = iaint2(1)
  Do i = 1, np0
    dpertp(i) = 1.
  End Do
  nadd = 0
  tau0 = arpar1(1)
  Do i = np0 + 1, np0 + nadd
    itypar(i) = 42
    dpertp(i) = 1./float(nadd)
    gxar(i) = 5.*(1.-2.*ranart(nseed))
    gyar(i) = 5.*(1.-2.*ranart(nseed))
    gzar(i) = 2.*(1.-2.*ranart(nseed))
    ftar(i) = 0.
    pxar(i) = 1.
    pyar(i) = 0.
    pzar(i) = 1.
    xmar(i) = xmd
    pear(i) = sqrt(pxar(i)**2+pyar(i)**2+pzar(i)**2+xmar(i)**2)
    rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
    vx = pxar(i)/pear(i)
    vy = pyar(i)/pear(i)
    taui = ftar(i) + tau0
    ftar(i) = taui*cosh(rap)
    gxar(i) = gxar(i) + vx*tau0*cosh(rap)
    gyar(i) = gyar(i) + vy*tau0*cosh(rap)
    gzar(i) = taui*sinh(rap) + gzar(i)
    zsmear = sngl(smearh)*(2.*ranart(nseed)-1.)
    gzar(i) = gzar(i) + zsmear
  End Do
  iaint2(1) = iaint2(1) + nadd
  If (nadd>=1 .And. idpert/=1 .And. idpert/=2) Then
    Write (16, *) 'IDPERT must be 1 or 2 to add initial hadrons,      set NPERTD to 0 if you do not need perturbative deuterons'
    Stop
  End If
  If (iaint2(1)>maxstr) Then
    Write (16, *) 'Too many initial hadrons, array size is exceeded!'
    Stop
  End If
  Return
End Subroutine addhad
