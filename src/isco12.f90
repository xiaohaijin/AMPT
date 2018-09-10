Subroutine isco12(i, j, allok, tm, t1, t2)
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
  If (ft(i)>=ft(j)) Then
     i1 = j
     i2 = i
     isign = -1
  Else
     i1 = i
     i2 = j
     isign = 1
  End If
  If (allok) Then
     tm = tlarge
     t1 = ft(i1)
     vx1 = vx(i1)
     vy1 = vy(i1)
     vz1 = vz(i1)
     t2 = ft(i2)
     dvx = vx(i2) - vx1
     dvy = vy(i2) - vy1
     dvz = vz(i2) - vz1
     dt = t2 - t1
     Do ii = -1, 1
        Do jj = -1, 1
           Do kk = -1, 1
              allokp = .True.
              dx = gx(i2) - gx(i1) - vx1*dt
              dy = gy(i2) - gy(i1) - vy1*dt
              dz = gz(i2) - gz(i1) - vz1*dt
              dx = dx + ii*10D0*size1
              dy = dy + jj*10D0*size2
              dz = dz + kk*10D0*size3
              vp = dvx*dx + dvy*dy + dvz*dz
              allokp = allokp .And. vp < 0D0
              If (allokp) Then
                 v2 = dvx*dvx + dvy*dvy + dvz*dvz
                 If (v2==0D0) Then
                    tmp = tlarge
                 Else
                    tmp = t2 - vp/v2
                 End If
                 allokp = allokp .And. tmp > t1 .And. tmp > t2
              End If
              If (allokp) Then
                 dgx = dx - dvx*t2
                 dgy = dy - dvy*t2
                 dgz = dz - dvz*t2
                 dm2 = -v2*tmp**2 + dgx*dgx + dgy*dgy + dgz*dgz
                 allokp = allokp .And. dm2 < cutof2
              End If
              If (allokp .And. tmp<tm) Then
                 tm = tmp
                 jxa = isign*ii
                 jya = isign*jj
                 jza = isign*kk
              End If
           End Do
        End Do
     End Do
     If (tm==tlarge) Then
        allok = .False.
     End If
  End If
  If (allok) Then
     e1 = e(i1)
     px1 = px(i1)
     py1 = py(i1)
     pz1 = pz(i1)
     e2 = e(i2)
     px2 = px(i2)
     py2 = py(i2)
     pz2 = pz(i2)
     rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
     allok = allok .And. rts2 > rscut2
  End If
  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else
     t1 = tm
     t2 = tm
  End If
  Return
End Subroutine isco12
