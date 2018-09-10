Subroutine arini2(k)
  Parameter (maxstr=150001, maxr=1)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  multi1(k) = iaint2(1)
  Do i = 1, multi1(k)
    ityp1(i, k) = itypar(i)
    gx1(i, k) = gxar(i)
    gy1(i, k) = gyar(i)
    gz1(i, k) = gzar(i)
    ft1(i, k) = ftar(i)
    px1(i, k) = pxar(i)
    py1(i, k) = pyar(i)
    pz1(i, k) = pzar(i)
    ee1(i, k) = pear(i)
    xm1(i, k) = xmar(i)
    dpp1(i, k) = dpertp(i)
  End Do
  Do ip = 1, maxstr
    tfdcy(ip) = ntmax*dt
    tft(ip) = ntmax*dt
  End Do
  Do irun = 1, maxr
    Do ip = 1, maxstr
      tfdpi(ip, irun) = ntmax*dt
    End Do
  End Do
  Return
End Subroutine arini2
