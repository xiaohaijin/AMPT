Subroutine isco11(i, j, allok, tm, t1, t2)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /aurec1/jxa, jya, jza
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  Logical allok, allokp
  allok = last(i) /= j .Or. last(j) /= i
  tm = tlarge
  If (allok) Then
     Do ii = -1, 1
        Do jj = -1, 1
           Do kk = -1, 1
              allokp = .True.
              i1 = i
              i2 = j
              p4 = ft(j) - ft(i)
              p1 = gx(j) - gx(i)
              p2 = gy(j) - gy(i)
              p3 = gz(j) - gz(i)
              p1 = p1 + ii*10D0*size1
              p2 = p2 + jj*10D0*size2
              p3 = p3 + kk*10D0*size3
              q4 = e(i)
              q1 = px(i)
              q2 = py(i)
              q3 = pz(i)
              r4 = e(j)
              r1 = px(j)
              r2 = py(j)
              r3 = pz(j)
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
              If (allokp) Then
                 vp = a*d - b*ee
                 allokp = allokp .And. vp < 0D0
              End If
              If (allokp) Then
                 dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
                 allokp = allokp .And. dm2 < cutof2
              End If
              If (allokp) Then
                 tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
                 tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
                 If (iordsc==20) Then
                    tmp = min(tc1, tc2)
                 Else If (iordsc==21) Then
                    tmp = 0.5D0*(tc1+tc2)
                 Else
                    tmp = max(tc1, tc2)
                 End If
                 allokp = allokp .And. tmp > ft(i) .And. tmp > ft(j)
              End If
              If (allokp .And. tmp<tm) Then
                 tm = tmp
                 jxa = ii
                 jya = jj
                 jza = kk
                 ha = h
                 tc1a = tc1
                 tc2a = tc2
              End If
           End Do
        End Do
     End Do
     If (tm==tlarge) Then
        allok = .False.
     End If
  End If
  If (allok) Then
     q4 = e(i1)
     q1 = px(i1)
     q2 = py(i1)
     q3 = pz(i1)
     r4 = e(i2)
     r1 = px(i2)
     r2 = py(i2)
     r3 = pz(i2)
     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2
     allok = allok .And. rts2 > rscut2
  End If
  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (ha>0D0) Then
     t1 = tc2a
     t2 = tc1a
  Else
     t1 = tc1a
     t2 = tc2a
  End If
  Return
End Subroutine isco11
