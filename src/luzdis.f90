Subroutine luzdis(kfl1, kfl2, pr, z)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  kfla = iabs(kfl1)
  kflb = iabs(kfl2)
  kflh = kfla
  If (kfla>=10) kflh = mod(kfla/1000, 10)
  If (mstj(11)==1 .Or. (mstj(11)==3 .And. kflh<=3)) Then
    fa = parj(41)
    If (mstj(91)==1) fa = parj(43)
    If (kflb>=10) fa = fa + parj(45)
    fb = parj(42)*pr
    If (mstj(91)==1) fb = parj(44)*pr
    fc = 1.
    If (kfla>=10) fc = fc - parj(45)
    If (kflb>=10) fc = fc + parj(45)
    mc = 1
    If (abs(fc-1.)>0.01) mc = 2
    If (fa<0.02) Then
      ma = 1
      zmax = 1.
      If (fc>fb) zmax = fb/fc
    Else If (abs(fc-fa)<0.01) Then
      ma = 2
      zmax = fb/(fb+fc)
    Else
      ma = 3
      zmax = 0.5*(fb+fc-sqrt((fb-fc)**2+4.*fa*fb))/(fc-fa)
      If (zmax>0.99 .And. fb>100.) zmax = 1. - fa/fb
    End If
    mmax = 2
    If (zmax<0.1) Then
      mmax = 1
      zdiv = 2.75*zmax
      If (mc==1) Then
        fint = 1. - log(zdiv)
      Else
        zdivc = zdiv**(1.-fc)
        fint = 1. + (1.-1./zdivc)/(fc-1.)
      End If
    Else If (zmax>0.85 .And. fb>1.) Then
      mmax = 3
      fscb = sqrt(4.+(fc/fb)**2)
      zdiv = fscb - 1./zmax - (fc/fb)*log(zmax*0.5*(fscb+fc/fb))
      If (ma>=2) zdiv = zdiv + (fa/fb)*log(1.-zmax)
      zdiv = min(zmax, max(0.,zdiv))
      fint = 1. + fb*(1.-zdiv)
    End If
    100 z = rlu(0)
    fpre = 1.
    If (mmax==1) Then
      If (fint*rlu(0)<=1.) Then
        z = zdiv*z
      Else If (mc==1) Then
        z = zdiv**z
        fpre = zdiv/z
      Else
        z = 1./(zdivc+z*(1.-zdivc))**(1./(1.-fc))
        fpre = (zdiv/z)**fc
      End If
    Else If (mmax==3) Then
      If (fint*rlu(0)<=1.) Then
        z = zdiv + log(z)/fb
        fpre = exp(fb*(z-zdiv))
      Else
        z = zdiv + z*(1.-zdiv)
      End If
    End If
    If (z<=fb/(50.+fb) .Or. z>=1.) Goto 100
    fval = (zmax/z)**fc*exp(fb*(1./zmax-1./z))
    If (ma>=2) fval = ((1.-z)/(1.-zmax))**fa*fval
    If (fval<rlu(0)*fpre) Goto 100
  Else
    fc = parj(50+max(1,kflh))
    If (mstj(91)==1) fc = parj(59)
    110 z = rlu(0)
    If (fc>=0. .And. fc<=1.) Then
      If (fc>rlu(0)) z = 1. - z**(1./3.)
    Else If (fc>-1.) Then
      If (-4.*fc*z*(1.-z)**2<rlu(0)*((1.-z)**2-fc*z)**2) Goto 110
    Else
      If (fc>0.) z = 1. - z**(1./fc)
      If (fc<0.) z = z**(-1./fc)
    End If
  End If
  Return
End Subroutine luzdis
