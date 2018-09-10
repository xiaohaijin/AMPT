Subroutine momntm(px, py, pz, e)
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Parameter (pi=3.14159265358979D0)
  Common /rndm3/iseedp
  Save
  iseed = iseedp
  cost = 2D0*ran1(iseed) - 1D0
  sint = dsqrt(1D0-cost**2)
  phi = 2D0*pi*ran1(iseed)
  px = e*sint*cos(phi)
  py = e*sint*sin(phi)
  pz = e*cost
  Return
End Subroutine momntm
