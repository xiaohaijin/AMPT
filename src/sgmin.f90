Function sgmin(n)
  Save
  ga = 0.
  If (n<=2) Goto 20
  Do i = 1, n - 1
    z = i
    ga = ga + alog(z)
  End Do
  20 sgmin = ga
  Return
End Function sgmin
