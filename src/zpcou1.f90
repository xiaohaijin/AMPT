Subroutine zpcou1
  Implicit Double Precision (A-H, O-Z)
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /ana1/ts(12)
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
  Save
  dpt = 0.5D0
  dy2 = 1D0
  dy1 = 0.5D0
  dy = 0.2D0
  ntotal = nevnt*nsbrun
  Return
End Subroutine zpcou1
