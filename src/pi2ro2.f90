Subroutine pi2ro2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /rndf77/nseed
  Save
  If ((lb(i1)>=3 .And. lb(i1)<=5) .And. (lb(i2)>=3 .And. lb(i2)<=5)) Then
    iblock = 1850
    ei1 = 0.77
    ei2 = 0.77
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 25 + int(3*ranart(nseed))
  Else If ((lb(i1)>=25 .And. lb(i1)<=27) .And. (lb(i2)>=25 .And. lb(i2)<=27)) Then
    iblock = 1851
    lbb1 = 3 + int(3*ranart(nseed))
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = ap2
    ei2 = ap2
    If (lbb1==4) ei1 = ap1
    If (lbb2==4) ei2 = ap1
  End If
  Return
End Subroutine pi2ro2
