Function pyspen(xrein, ximin, ireim)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension b(0:14)
  Data b/1.000000E+00, -5.000000E-01, 1.666667E-01, 0.000000E+00, -3.333333E-02, 0.000000E+00, 2.380952E-02, 0.000000E+00, -3.333333E-02, 0.000000E+00, 7.575757E-02, 0.000000E+00, -2.531135E-01, 0.000000E+00, 1.166667E+00/
  xre = xrein
  xim = ximin
  If (abs(1.-xre)<1.E-6 .And. abs(xim)<1.E-6) Then
    If (ireim==1) pyspen = paru(1)**2/6.
    If (ireim==2) pyspen = 0.
    Return
  End If
  xmod = sqrt(xre**2+xim**2)
  If (xmod<1.E-6) Then
    If (ireim==1) pyspen = 0.
    If (ireim==2) pyspen = 0.
    Return
  End If
  xarg = sign(acos(xre/xmod), xim)
  sp0re = 0.
  sp0im = 0.
  sgn = 1.
  If (xmod>1.) Then
    algxre = log(xmod)
    algxim = xarg - sign(paru(1), xarg)
    sp0re = -paru(1)**2/6. - (algxre**2-algxim**2)/2.
    sp0im = -algxre*algxim
    sgn = -1.
    xmod = 1./xmod
    xarg = -xarg
    xre = xmod*cos(xarg)
    xim = xmod*sin(xarg)
  End If
  If (xre>0.5) Then
    algxre = log(xmod)
    algxim = xarg
    xre = 1. - xre
    xim = -xim
    xmod = sqrt(xre**2+xim**2)
    xarg = sign(acos(xre/xmod), xim)
    algyre = log(xmod)
    algyim = xarg
    sp0re = sp0re + sgn*(paru(1)**2/6.-(algxre*algyre-algxim*algyim))
    sp0im = sp0im - sgn*(algxre*algyim+algxim*algyre)
    sgn = -sgn
  End If
  xre = 1. - xre
  xim = -xim
  xmod = sqrt(xre**2+xim**2)
  xarg = sign(acos(xre/xmod), xim)
  zre = -log(xmod)
  zim = -xarg
  spre = 0.
  spim = 0.
  savere = 1.
  saveim = 0.
  Do i = 0, 14
    termre = (savere*zre-saveim*zim)/float(i+1)
    termim = (savere*zim+saveim*zre)/float(i+1)
    savere = termre
    saveim = termim
    spre = spre + b(i)*termre
    spim = spim + b(i)*termim
  End Do
  If (ireim==1) pyspen = sp0re + sgn*spre
  If (ireim==2) pyspen = sp0im + sgn*spim
  Return
End Function pyspen
