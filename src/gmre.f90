Function gmre(x)
  Implicit Double Precision (A-H, O-Z)
  Save
  z = x
  If (x>3.0D0) Goto 10
  z = x + 3.D0
  10 gmre = 0.5D0*dlog(2.D0*3.14159265D0/z) + z*dlog(z) - z + dlog(1.D0+1.D0/12.D0/z+1.D0/288.D0/z**2-139.D0/51840.D0/z**3-571.D0/2488320.D0/z**4)
  If (z==x) Goto 20
  gmre = gmre - dlog(z-1.D0) - dlog(z-2.D0) - dlog(z-3.D0)
  20 Continue
  Return
End Function gmre
