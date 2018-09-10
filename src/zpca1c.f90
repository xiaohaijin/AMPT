Subroutine zpca1c(p0, p1, p2, p3, ian)
  Implicit Double Precision (A-H, O-Z)
  Common /ana3/em(4, 4, 12)
  Dimension en(4)
  Save
  en(1) = p0
  en(2) = p1
  en(3) = p2
  en(4) = p3
  Do i = 1, 4
     Do j = 1, 4
        em(i, j, ian) = em(i, j, ian) + en(i)*en(j)/p0
     End Do
  End Do
  Return
End Subroutine zpca1c
