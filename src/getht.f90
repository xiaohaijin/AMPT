Subroutine getht(iscat, jscat, pp2, that)
  Implicit Double Precision (A-H, O-Z)
  Parameter (hbarc=0.197327054D0)
  Parameter (maxptn=400001)
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /anim/nevent, isoft, isflag, izpc
  External ran1
  Common /rndm3/iseedp
  Save
  iseed = iseedp
  xmu2 = (hbarc*xmu)**2
  xmp2 = xmp**2
  xm2 = xmu2 + xmp2
  rx = ran1(iseed)
  that = xm2*(1D0+1D0/((1D0-xm2/(4D0*pp2+xm2))*rx-1D0))
  If (izpc==100) that = -4D0*pp2*rx
  Return
End Subroutine getht
