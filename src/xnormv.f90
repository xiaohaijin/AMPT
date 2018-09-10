Subroutine xnormv(vx, vy, vz)
  Implicit Double Precision (A-H, O-Z)
  Save
  vv = dsqrt(vx**2+vy**2+vz**2)
  vx = vx/vv
  vy = vy/vv
  vz = vz/vv
  Return
End Subroutine xnormv
