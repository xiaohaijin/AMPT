Subroutine dmelangle(pxn, pyn, pzn, pfinal)
  Parameter (pi=3.1415926)
  Common /rndf77/nseed
  Save
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pzn = pfinal*c1
  pxn = pfinal*s1*ct1
  pyn = pfinal*s1*st1
  Return
End Subroutine dmelangle
