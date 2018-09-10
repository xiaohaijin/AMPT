Subroutine pbarfs(srt, npion, iseed)
  Parameter (pimass=0.140, pi=3.1415926)
  Dimension factor(6), pnpi(6)
  Common /rndf77/nseed
  Save
  factor(2) = 1.
  factor(3) = 1.17E-01
  factor(4) = 3.27E-03
  factor(5) = 3.58E-05
  factor(6) = 1.93E-07
  ene = (srt/pimass)**3/(6.*pi**2)
  Do n = 2, 6
    pnpi(n) = ene**n*factor(n)
  End Do
  pmax = max(pnpi(2), pnpi(3), pnpi(4), pnpi(5), pnpi(6))
  ntry = 0
  10 npion = 2 + int(5*ranart(nseed))
  If (npion>6) Goto 10
  thisp = pnpi(npion)/pmax
  ntry = ntry + 1
  If ((thisp<ranart(nseed)) .And. (ntry<=20)) Goto 10
  Return
End Subroutine pbarfs
