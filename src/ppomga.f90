Subroutine ppomga(srt, iseed, px, py, pz, dm1, pnx, pny, pnz, dm2, ppx, ppy, ppz, icou1)
  Common /table/xarray(0:1000), earray(0:1000)
  Common /rndf77/nseed
  Save
  ntrym = 0
  icou1 = 0
  pi = 3.1415926
  amn = 938.925/1000.
  amp = 782./1000.
  dm1 = amn
  dm2 = amn
  v = 0.43
  w = -0.84
  ptmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2
  scheck = ptmax2
  If (scheck<0) Then
     Write (99, *) 'scheck30: ', scheck
     scheck = 0.
  End If
  ptmax = sqrt(scheck)*1./3.
7 pt = ptr(ptmax, iseed)
  pzmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2 - pt**2
  ntrym = ntrym + 1
  If ((pzmax2<0.) .And. ntrym<=100) Then
     Goto 7
  Else
     pzmax2 = 1.E-09
  End If
  pzmax = sqrt(pzmax2)
  xmax = 2.*pzmax/srt
  ntryx = 0
  fmax00 = 1.056
  x00 = 0.26
  If (abs(xmax)>0.26) Then
     f00 = fmax00
  Else
     f00 = 1. + v*abs(xmax) + w*xmax**2
  End If
9 x = xmax*(1.-2.*ranart(nseed))
  ntryx = ntryx + 1
  xratio = (1.+v*abs(x)+w*x**2)/f00
  If (xratio<ranart(nseed) .And. ntryx<=50) Goto 9
  pz = 0.5*srt*x
  fai = 2.*pi*ranart(nseed)
  px = pt*cos(fai)
  py = pt*sin(fai)
  ek = sqrt(dm1**2+pt**2+pz**2)
  eln = srt - ek
  If (eln<=0) Then
     icou1 = -1
     Return
  End If
  bx = -px/eln
  by = -py/eln
  bz = -pz/eln
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck31: ', scheck
     Stop
  End If
  ga = 1./sqrt(scheck)
  elnc = eln/ga
  pn2 = ((elnc**2+dm2**2-amp**2)/(2.*elnc))**2 - dm2**2
  If (pn2<=0) Then
     icou1 = -1
     Return
  End If
  pn = sqrt(pn2)
  xptr = 0.33*pn
  pnt = ptr(xptr, iseed)
  fain = 2.*pi*ranart(nseed)
  pnx = pnt*cos(fain)
  pny = pnt*sin(fain)
  sig = 1
  If (x>0) sig = -1
  scheck = pn**2 - pnt**2
  If (scheck<0) Then
     Write (99, *) 'scheck32: ', scheck
     scheck = 0.
  End If
  pnz = sig*sqrt(scheck)
  en = sqrt(dm2**2+pnx**2+pny**2+pnz**2)
  ppx = -pnx
  ppy = -pny
  ppz = -pnz
  ep = sqrt(amp**2+ppx**2+ppy**2+ppz**2)
  pbeta = pnx*bx + pny*by + pnz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  pnx = bx*trans0 + pnx
  pny = by*trans0 + pny
  pnz = bz*trans0 + pnz
  If (ep==0.) ep = 1.E-09
  pbeta = ppx*bx + ppy*by + ppz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+ep)
  ppx = bx*trans0 + ppx
  ppy = by*trans0 + ppy
  ppz = bz*trans0 + ppz
  Return
End Subroutine ppomga
