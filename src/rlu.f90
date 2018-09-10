Function rlu(idum)
  Common /ludatr/mrlu(6), rrlu(100)
  Save /ludatr/
  Equivalence (mrlu1, mrlu(1)), (mrlu2, mrlu(2)), (mrlu3, mrlu(3)), (mrlu4, mrlu(4)), (mrlu5, mrlu(5)), (mrlu6, mrlu(6)), (rrlu98, rrlu(98)), (rrlu99, rrlu(99)), (rrlu00, rrlu(100))
  If (mrlu2==0) Then
    ij = mod(mrlu1/30082, 31329)
    kl = mod(mrlu1, 30082)
    i = mod(ij/177, 177) + 2
    j = mod(ij, 177) + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl, 169)
    Do ii = 1, 97
      s = 0.
      t = 0.5
      Do jj = 1, 24
        m = mod(mod(i*j,179)*k, 179)
        i = j
        j = k
        k = m
        l = mod(53*l+1, 169)
        If (mod(l*m,64)>=32) s = s + t
        t = 0.5*t
      End Do
      rrlu(ii) = s
    End Do
    twom24 = 1.
    Do i24 = 1, 24
      twom24 = 0.5*twom24
    End Do
    rrlu98 = 362436.*twom24
    rrlu99 = 7654321.*twom24
    rrlu00 = 16777213.*twom24
    mrlu2 = 1
    mrlu3 = 0
    mrlu4 = 97
    mrlu5 = 33
  End If
  130 runi = rrlu(mrlu4) - rrlu(mrlu5)
  If (runi<0.) runi = runi + 1.
  rrlu(mrlu4) = runi
  mrlu4 = mrlu4 - 1
  If (mrlu4==0) mrlu4 = 97
  mrlu5 = mrlu5 - 1
  If (mrlu5==0) mrlu5 = 97
  rrlu98 = rrlu98 - rrlu99
  If (rrlu98<0.) rrlu98 = rrlu98 + rrlu00
  runi = runi - rrlu98
  If (runi<0.) runi = runi + 1.
  If (runi<=0 .Or. runi>=1.) Goto 130
  mrlu3 = mrlu3 + 1
  If (mrlu3==1000000000) Then
    mrlu2 = mrlu2 + 1
    mrlu3 = 0
  End If
  rlu = runi
  Return
End Function rlu
