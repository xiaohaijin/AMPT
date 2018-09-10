Subroutine boosti
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para6/centy
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /lor/enenew, pxnew, pynew, pznew
  Save
  External lorenz
  bex = 0D0
  bey = 0D0
  bez = -tanh(centy)
  Do i = 1, mul
     px1 = gx0(i)
     py1 = gy0(i)
     pz1 = gz0(i)
     e1 = ft0(i)
     Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
     gx0(i) = pxnew
     gy0(i) = pynew
     gz0(i) = pznew
     ft0(i) = enenew
     px1 = px0(i)
     py1 = py0(i)
     pz1 = pz0(i)
     e1 = e0(i)
     Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
     px0(i) = pxnew
     py0(i) = pynew
     pz0(i) = pznew
     e0(i) = enenew
  End Do
  Return
End Subroutine boosti
