Subroutine chkout(l, t, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Save
  m1 = 11
  m2 = 11
  m3 = 11
  Call chkcel(l, m1, m2, m3, t, tmin, nc)
  Do i = 1, 10
     Do j = 1, 10
        Do k = 1, 10
           If (i==1 .Or. i==10 .Or. j==1 .Or. j==10 .Or. k==1 .Or. k==10) Call chkcel(l, i, j, k, t, tmin, nc)
        End Do
     End Do
  End Do
  Return
End Subroutine chkout
