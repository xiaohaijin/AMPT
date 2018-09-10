Function pyi3au(be, eps, ireim)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  If (eps<1.) ga = 0.5*(1.+sqrt(1.-eps))
  If (eps<0.) Then
    f3re = pyspen((ga-1.)/(ga+be-1.), 0., 1) - pyspen(ga/(ga+be-1.), 0., 1) + pyspen((be-ga)/be, 0., 1) - pyspen((be-ga)/(be-1.), 0., 1) + (log(be)**2-log(be-1.)**2)/2. + log(ga)*log((ga+be-1.)/be) + log(ga-1.)*log((be-1.)/(ga+be-1.))
    f3im = 0.
  Else If (eps<1.) Then
    f3re = pyspen((ga-1.)/(ga+be-1.), 0., 1) - pyspen(ga/(ga+be-1.), 0., 1) + pyspen(ga/(ga-be), 0., 1) - pyspen((ga-1.)/(ga-be), 0., 1) + log(ga/(1.-ga))*log((ga+be-1.)/(be-ga))
    f3im = -paru(1)*log((ga+be-1.)/(be-ga))
  Else
    rsq = eps/(eps-1.+(2.*be-1.)**2)
    rcthe = rsq*(1.-2.*be/eps)
    rsthe = sqrt(rsq-rcthe**2)
    rcphi = rsq*(1.+2.*(be-1.)/eps)
    rsphi = sqrt(rsq-rcphi**2)
    r = sqrt(rsq)
    the = acos(rcthe/r)
    phi = acos(rcphi/r)
    f3re = pyspen(rcthe, rsthe, 1) + pyspen(rcthe, -rsthe, 1) - pyspen(rcphi, rsphi, 1) - pyspen(rcphi, -rsphi, 1) + (phi-the)*(phi+the-paru(1))
    f3im = pyspen(rcthe, rsthe, 2) + pyspen(rcthe, -rsthe, 2) - pyspen(rcphi, rsphi, 2) - pyspen(rcphi, -rsphi, 2)
  End If
  If (ireim==1) pyi3au = 2./(2.*be-1.)*f3re
  If (ireim==2) pyi3au = 2./(2.*be-1.)*f3im
  Return
End Function pyi3au
