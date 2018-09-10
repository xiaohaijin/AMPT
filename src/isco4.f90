Subroutine isco4(i, j, allok, tm, t1, t2)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  Logical allok
  allok = last(i) /= j .Or. last(j) /= i
  icels1 = icels(i)
  ii1 = icels1/10000
  jj1 = (icels1-ii1*10000)/100
  kk1 = icels1 - ii1*10000 - jj1*100
  icels2 = icels(j)
  ii2 = icels2/10000
  jj2 = (icels2-ii2*10000)/100
  kk2 = icels2 - ii2*10000 - jj2*100
  i1 = i
  i2 = j
  p4 = ft(i2) - ft(i1)
  p1 = gx(i2) - gx(i1)
  p2 = gy(i2) - gy(i1)
  p3 = gz(i2) - gz(i1)
  If (ii1-ii2>5) Then
     p1 = p1 + 10D0*size1
  Else If (ii1-ii2<-5) Then
     p1 = p1 - 10D0*size1
  End If
  If (jj1-jj2>5) Then
     p2 = p2 + 10D0*size2
  Else If (jj1-jj2<-5) Then
     p2 = p2 - 10D0*size2
  End If
  If (kk1-kk2>5) Then
     p3 = p3 + 10D0*size3
  Else If (kk1-kk2<-5) Then
     p3 = p3 - 10D0*size3
  End If
  q4 = e(i1)
  q1 = px(i1)
  q2 = py(i1)
  q3 = pz(i1)
  r4 = e(i2)
  r1 = px(i2)
  r2 = py(i2)
  r3 = pz(i2)
  a = p4*q4 - p1*q1 - p2*q2 - p3*q3
  b = p4*r4 - p1*r1 - p2*r2 - p3*r3
  c = q4*q4 - q1*q1 - q2*q2 - q3*q3
  d = r4*r4 - r1*r1 - r2*r2 - r3*r3
  ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
  f = p4*p4 - p1*p1 - p2*p2 - p3*p3
  h = a + b
  If (h>0D0) Then
     g = a
     a = -b
     b = -g
     g = c
     c = d
     d = g
     i1 = j
     i2 = i
  End If
  If (allok) Then
     vp = a*d - b*ee
     allok = allok .And. vp < 0D0
  End If
  If (allok) Then
     dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
     allok = allok .And. dm2 < cutof2
  End If
  If (allok) Then
     tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
     tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
     tm = 0.5D0*(tc1+tc2)
     allok = allok .And. tm > ft(i) .And. tm > ft(j)
  End If
  If (allok) Then
     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2
     allok = allok .And. rts2 > rscut2
  End If
  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tm
     t2 = tm
  Else
     t1 = tm
     t2 = tm
  End If
  Return
End Subroutine isco4
