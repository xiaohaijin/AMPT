Subroutine energy(e, temp)
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /rndm3/iseedp
  Save
  iseed = iseedp
1000 Continue
  e = ran1(iseed)
  e = e*ran1(iseed)
  e = e*ran1(iseed)
  If (e<=0D0) Goto 1000
  e = -temp*log(e)
  If (ran1(iseed)>exp((e-dsqrt(e**2+xmp**2))/temp)) Then
     Goto 1000
  End If
  Return
End Subroutine energy
