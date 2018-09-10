Subroutine aran9(qran, ndim)
  Dimension qran(10)
  Common /sedvax/num1
  Save
  Do i = 1, ndim
    qran(i) = ranart(num1)
  End Do
  Return
End Subroutine aran9
