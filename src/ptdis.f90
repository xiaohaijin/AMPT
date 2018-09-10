Real Function ptdis(x)
  Save
  b = 3.78
  c = 0.47
  d = 3.60
  ptdis = 1./(2.*b)*(1.-exp(-b*x**2)) - c/d*x*exp(-d*x) - c/d**2*(exp(-d*x)-1.)
  Return
End Function ptdis
