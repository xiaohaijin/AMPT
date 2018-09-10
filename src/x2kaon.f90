Real Function x2kaon(srt)
  Save
  smin = 2.8639
  x2kaon = 0.0000001
  If (srt<smin) Return
  sigma1 = 2.8
  sigma2 = 7.7
  sigma3 = 3.9
  x = srt**2/smin**2 + 0.0000001
  f1 = (1.+1./sqrt(x))*alog(x) - 4.*(1.-1./sqrt(x))
  f2 = 1. - (1./sqrt(x))*(1.+alog(sqrt(x)))
  f3 = ((x-1.)/x**2)**3.5
  x2kaon = (1.-1./x)**3*(sigma1*f1+sigma2*f2) + sigma3*f3
  Return
End Function x2kaon
