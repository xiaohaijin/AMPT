Subroutine cropro(vx1, vy1, vz1, vx2, vy2, vz2)
  Implicit Double Precision (A-H, O-Z)
  Common /cprod/vx3, vy3, vz3
  Save
  vx3 = vy1*vz2 - vz1*vy2
  vy3 = vz1*vx2 - vx1*vz2
  vz3 = vx1*vy2 - vy1*vx2
  Return
End Subroutine cropro
