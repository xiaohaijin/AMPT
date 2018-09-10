Subroutine zprota(xn1, xn2, xn3, theta, v1, v2, v3)
  Implicit Double Precision (A-H, O-Z)
  Save
  vx = v1
  vy = v2
  vz = v3
  c = cos(theta)
  omc = 1D0 - c
  s = sin(theta)
  a11 = xn1**2*omc + c
  a12 = xn1*xn2*omc - s*xn3
  a13 = xn1*xn3*omc + s*xn2
  a21 = xn1*xn2*omc + s*xn3
  a22 = xn2**2*omc + c
  a23 = xn2*xn3*omc - s*xn1
  a31 = xn1*xn3*omc - s*xn2
  a32 = xn3*xn2*omc + s*xn1
  a33 = xn3**2*omc + c
  v1 = vx*a11 + vy*a12 + vz*a13
  v2 = vx*a21 + vy*a22 + vz*a23
  v3 = vx*a31 + vy*a32 + vz*a33
  Return
End Subroutine zprota
