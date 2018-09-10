Subroutine opioet(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957, etam=0.5475, aomega=0.782)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /rndf77/nseed
  Save
  If ((lb(i1)>=3 .And. lb(i1)<=5 .And. lb(i2)==28) .Or. (lb(i2)>=3 .And. lb(i2)<=5 .And. lb(i1)==28)) Then
    iblock = 1890
    ei1 = aomega
    ei2 = etam
    lbb1 = 28
    lbb2 = 0
  Else If ((lb(i1)==28 .And. lb(i2)==0) .Or. (lb(i1)==0 .And. lb(i2)==28)) Then
    iblock = 1891
    lbb1 = 28
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = aomega
    ei2 = ap2
    If (lbb2==4) ei2 = ap1
  End If
  Return
End Subroutine opioet
