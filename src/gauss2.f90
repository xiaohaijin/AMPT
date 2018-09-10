Function gauss2(f, a, b, eps)
  External f
  Dimension w(12), x(12)
  Save
  Data const/1.0E-12/
  Data w/0.1012285, .2223810, .3137067, .3623838, .0271525, .0622535, 0.0951585, .1246290, .1495960, .1691565, .1826034, .1894506/
  Data x/0.9602899, .7966665, .5255324, .1834346, .9894009, .9445750, 0.8656312, .7554044, .6178762, .4580168, .2816036, .0950125/
  delta = const*abs(a-b)
  gauss2 = 0.0
  aa = a
  5 y = b - aa
  If (abs(y)<=delta) Return
  2 bb = aa + y
  c1 = 0.5*(aa+bb)
  c2 = c1 - aa
  s8 = 0.0
  s16 = 0.0
  Do i = 1, 4
    u = x(i)*c2
    s8 = s8 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  Do i = 5, 12
    u = x(i)*c2
    s16 = s16 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  s8 = s8*c2
  s16 = s16*c2
  If (abs(s16-s8)>eps*(1.+abs(s16))) Goto 4
  gauss2 = gauss2 + s16
  aa = bb
  Goto 5
  4 y = 0.5*y
  If (abs(y)>delta) Goto 2
  Write (6, 7)
  gauss2 = 0.0
  Return
  7 Format (1X, 'GAUSS2....TOO HIGH ACURACY REQUIRED')
End Function gauss2
