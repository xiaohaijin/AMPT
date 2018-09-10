Double Precision Function ftime1(iseed)
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Parameter (hbarc=0.197327054D0)
  Common /par1/formt
  Save
  aa = hbarc/formt
  ftime1 = aa*dsqrt(1D0/ran1(iseed)-1D0)
  Return
End Function ftime1
