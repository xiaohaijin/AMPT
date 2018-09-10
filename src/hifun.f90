Subroutine hifun(i, xmin, xmax, fhb)
  Common /hijhb/rr(10, 201), xx(10, 201)
  External fhb
  Save
  fnorm = gauss1(fhb, xmin, xmax, 0.001)
  Do j = 1, 201
    xx(i, j) = xmin + (xmax-xmin)*(j-1)/200.0
    xdd = xx(i, j)
    rr(i, j) = gauss1(fhb, xmin, xdd, 0.001)/fnorm
  End Do
  Return
End Subroutine hifun
