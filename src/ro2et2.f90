Subroutine ro2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (etam=0.5475, arho=0.77)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /rndf77/nseed
  Save
  If (lb(i1)>=25 .And. lb(i1)<=27 .And. lb(i2)>=25 .And. lb(i2)<=27) Then
    iblock = 1895
    ei1 = etam
    ei2 = etam
    lbb1 = 0
    lbb2 = 0
  Else If (lb(i1)==0 .And. lb(i2)==0) Then
    iblock = 1896
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 25 + int(3*ranart(nseed))
    ei1 = arho
    ei2 = arho
  End If
  Return
End Subroutine ro2et2
