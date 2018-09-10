Subroutine posit3(x, y, z)
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /rndm3/iseedp
  Save
  iseed = iseedp
  x = 2D0*ran1(iseed) - 1D0
  y = 2D0*ran1(iseed) - 1D0
  z = 2D0*ran1(iseed) - 1D0
  x = x*5D0*size1
  y = y*5D0*size2
  z = z*5D0*size3
  Return
End Subroutine posit3
