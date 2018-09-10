Subroutine posit1(x, y, r0)
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Common /rndm3/iseedp
  Save
  iseed = iseedp
10 x = 2D0*ran1(iseed) - 1D0
  y = 2D0*ran1(iseed) - 1D0
  If (x**2+y**2>1D0) Goto 10
  x = x*r0
  y = y*r0
  Return
End Subroutine posit1
