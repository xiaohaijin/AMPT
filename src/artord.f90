Subroutine artord
  Parameter (maxstr=150001, maxr=1)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Dimension dptemp(maxstr)
  Dimension ityp0(maxstr), gx0(maxstr), gy0(maxstr), gz0(maxstr), ft0(maxstr), px0(maxstr), py0(maxstr), pz0(maxstr), ee0(maxstr), xm0(maxstr)
  Dimension indx(maxstr)
  External arindx
  Save
  npar = 0
  np = iaint2(1)
  Do i = 1, np
    ityp0(i) = itypar(i)
    gx0(i) = gxar(i)
    gy0(i) = gyar(i)
    gz0(i) = gzar(i)
    ft0(i) = ftar(i)
    px0(i) = pxar(i)
    py0(i) = pyar(i)
    pz0(i) = pzar(i)
    ee0(i) = pear(i)
    xm0(i) = xmar(i)
    dptemp(i) = dpertp(i)
  End Do
  Call arindx(maxstr, np, ft0, indx)
  Do i = 1, np
    npar = npar + 1
    itypar(npar) = ityp0(indx(i))
    gxar(npar) = gx0(indx(i))
    gyar(npar) = gy0(indx(i))
    gzar(npar) = gz0(indx(i))
    ftar(npar) = ft0(indx(i))
    pxar(npar) = px0(indx(i))
    pyar(npar) = py0(indx(i))
    pzar(npar) = pz0(indx(i))
    pear(npar) = ee0(indx(i))
    xmar(npar) = xm0(indx(i))
    dpertp(npar) = dptemp(indx(i))
  End Do
  iaint2(1) = npar
  Return
End Subroutine artord
