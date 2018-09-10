Subroutine dchout(l, ii, t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  External integ
  tt = ft(l)
  td = t - size
  x1 = gx(l) + vx(l)*(t-tt)
  x2 = gy(l) + vy(l)*(t-tt)
  x3 = gz(l) + vz(l)*(t-tt)
  If (td<=0D0) Then
     i1 = integ(x1/size1) + 6
     i2 = integ(x2/size2) + 6
     i3 = integ(x3/size3) + 6
     If (integ(x1/size1)==x1/size1 .And. vx(l)<0D0) i1 = i1 - 1
     If (integ(x2/size2)==x2/size2 .And. vy(l)<0D0) i2 = i2 - 1
     If (integ(x3/size3)==x3/size3 .And. vz(l)<0D0) i3 = i3 - 1
  Else
     i1 = integ(x1/(size1+v1*td)) + 6
     i2 = integ(x2/(size2+v2*td)) + 6
     i3 = integ(x3/(size3+v3*td)) + 6
     If (integ(x1/(size1+v1*td))==x1/(size1+v1*td) .And. vx(l)<(i1-6)*v1) i1 = i1 - 1
     If (integ(x2/(size2+v2*td))==x2/(size2+v2*td) .And. vy(l)<(i2-6)*v2) i2 = i2 - 1
     If (integ(x3/(size3+v3*td))==x3/(size3+v3*td) .And. vz(l)<(i3-6)*v3) i3 = i3 - 1
  End If
  If (ii==1) Then
     i = 9
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  If (ii==2) Then
     i = 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  If (ii==3) Then
     j = 9
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  If (ii==4) Then
     j = 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  If (ii==5) Then
     k = 9
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  If (ii==6) Then
     k = 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If
  Return
End Subroutine dchout
