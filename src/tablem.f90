Subroutine tablem
  Common /table/xarray(0:1000), earray(0:1000)
  Save
  ptmax = 2.01
  anorm = ptdis(ptmax)
  Do l = 0, 200
     x = 0.01*float(l+1)
     rr = ptdis(x)/anorm
     earray(l) = rr
     xarray(l) = x
  End Do
  Return
End Subroutine tablem
