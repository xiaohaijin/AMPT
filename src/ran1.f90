Double Precision Function ran1(idum)
  Implicit Double Precision (A-H, O-Z)
  Dimension r(97)
  Common /rndm1/number
  Parameter (m1=259200, ia1=7141, ic1=54773, rm1=1D0/m1)
  Parameter (m2=134456, ia2=8121, ic2=28411, rm2=1D0/m2)
  Parameter (m3=243000, ia3=4561, ic3=51349)
  Save
  Data iff/0/
  If (idum<0 .Or. iff==0) Then
     iff = 1
     ix1 = mod(ic1-idum, m1)
     ix1 = mod(ia1*ix1+ic1, m1)
     ix2 = mod(ix1, m2)
     ix1 = mod(ia1*ix1+ic1, m1)
     ix3 = mod(ix1, m3)
     Do j = 1, 97
        ix1 = mod(ia1*ix1+ic1, m1)
        ix2 = mod(ia2*ix2+ic2, m2)
        r(j) = (dble(ix1)+dble(ix2)*rm2)*rm1
     End Do
     idum = 1
  End If
  ix1 = mod(ia1*ix1+ic1, m1)
  ix2 = mod(ia2*ix2+ic2, m2)
  ix3 = mod(ia3*ix3+ic3, m3)
  j = 1 + (97*ix3)/m3
  If (j>97 .Or. j<1) Print *, 'In zpc ran1, j<1 or j>97', j
  ran1 = r(j)
  r(j) = (dble(ix1)+dble(ix2)*rm2)*rm1
  number = number + 1
  Return
End Function ran1
