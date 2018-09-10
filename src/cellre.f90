Subroutine cellre(i, t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /aurec1/jxa, jya, jza
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist2/icell, icel(10, 10, 10)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  Logical good
  External integ
  t0 = t
1000 Continue
  If (iconfg==3 .Or. iconfg==5) Then
     k = mod(icsta(i), 10)
     If (k==1) Then
        gx(i) = gx(i) - 10D0*size1
        dgxa(i) = dgxa(i) + 10D0*size1
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgxa(ii) = dgxa(ii) - 10D0*size1
           End If
        End Do
     End If
     If (k==2) Then
        gx(i) = gx(i) + 10D0*size1
        dgxa(i) = dgxa(i) - 10D0*size1
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgxa(ii) = dgxa(ii) + 10D0*size1
           End If
        End Do
     End If
     If (k==3) Then
        gy(i) = gy(i) - 10D0*size2
        dgya(i) = dgya(i) + 10D0*size2
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgya(ii) = dgya(ii) - 10D0*size2
           End If
        End Do
     End If
     If (k==4) Then
        gy(i) = gy(i) + 10D0*size2
        dgya(i) = dgya(i) - 10D0*size2
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgya(ii) = dgya(ii) + 10D0*size2
           End If
        End Do
     End If
     If (iconfg==5) Then
        If (k==5) Then
           gz(i) = gz(i) - 10D0*size3
           dgza(i) = dgza(i) + 10D0*size3
           Do ii = 1, ichkpt
              If (next(ii)==i) Then
                 dgza(ii) = dgza(ii) - 10D0*size3
              End If
           End Do
        End If
        If (k==6) Then
           gz(i) = gz(i) + 10D0*size3
           dgza(i) = dgza(i) - 10D0*size3
           Do ii = 1, ichkpt
              If (next(ii)==i) Then
                 dgza(ii) = dgza(ii) + 10D0*size3
              End If
           End Do
        End If
     End If
  Else
     icels0 = icels(i)
     i1 = icels0/10000
     i2 = (icels0-i1*10000)/100
     i3 = icels0 - i1*10000 - i2*100
     If (i1>=1 .And. i1<=10 .And. i2>=1 .And. i2<=10 .And. i3>=1 .And. i3<=10) Then
        If (icel(i1,i2,i3)==i) icel(i1, i2, i3) = nic(i)
        Call oldcre(i)
        k = mod(icsta(i), 10)
        If (iconfg==1) Then
           good = (i1==1 .And. k==2) .Or. (i1==10 .And. k==1) .Or. (i2==1 .And. k==4) .Or. (i2==10 .And. k==3) .Or. (i3==1 .And. k==6) .Or. (i3==10 .And. k==5)
        End If
        If (iconfg==2) Then
           good = (i3==1 .And. k==6) .Or. (i3==10 .And. k==5)
        End If
        If (good) Then
           Call newcre(i, icell)
           icels(i) = 111111
        Else
           If (k==1) i1 = i1 + 1
           If (k==2) i1 = i1 - 1
           If (k==3) i2 = i2 + 1
           If (k==4) i2 = i2 - 1
           If (k==5) i3 = i3 + 1
           If (k==6) i3 = i3 - 1
           If (iconfg==2 .Or. iconfg==4) Then
              If (i1==0) Then
                 i1 = 10
                 gx(i) = gx(i) + 10D0*size1
              End If
              If (i1==11) Then
                 i1 = 1
                 gx(i) = gx(i) - 10D0*size1
              End If
              If (i2==0) Then
                 i2 = 10
                 gy(i) = gy(i) + 10D0*size2
              End If
              If (i2==11) Then
                 i2 = 1
                 gy(i) = gy(i) - 10D0*size2
              End If
              If (iconfg==4) Then
                 If (i3==0) Then
                    i3 = 10
                    gz(i) = gz(i) + 10D0*size3
                 End If
                 If (i3==11) Then
                    i3 = 1
                    gz(i) = gz(i) - 10D0*size3
                 End If
              End If
           End If
           j = icel(i1, i2, i3)
           Call newcre(i, j)
           icel(i1, i2, i3) = j
           icels(i) = i1*10000 + i2*100 + i3
        End If
     Else
        If (icell==i) icell = nic(i)
        Call oldcre(i)
        k = mod(icsta(i), 10)
        ddt = t - ft(i)
        dtt = t - size
        If (dtt<=0D0) Then
           i1 = integ((gx(i)+vx(i)*ddt)/size1) + 6
           i2 = integ((gy(i)+vy(i)*ddt)/size2) + 6
           i3 = integ((gz(i)+vz(i)*ddt)/size3) + 6
        Else
           i1 = integ((gx(i)+vx(i)*ddt)/(size1+v1*dtt)) + 6
           i2 = integ((gy(i)+vy(i)*ddt)/(size2+v2*dtt)) + 6
           i3 = integ((gz(i)+vz(i)*ddt)/(size3+v3*dtt)) + 6
        End If
        If (k==1) i1 = 1
        If (k==2) i1 = 10
        If (k==3) i2 = 1
        If (k==4) i2 = 10
        If (k==5) i3 = 1
        If (k==6) i3 = 10
        j = icel(i1, i2, i3)
        Call newcre(i, j)
        icel(i1, i2, i3) = j
        icels(i) = i1*10000 + i2*100 + i3
     End If
  End If
  If (next(i)/=0) Then
     otmp = ot(next(i))
     ctmp = ct(next(i))
  End If
  If (i1==11 .And. i2==11 .And. i3==11) Then
     Call dchout(i, k, t)
  Else
     If (iconfg==1) Then
        Call dchin1(i, k, i1, i2, i3, t)
     Else If (iconfg==2) Then
        Call dchin2(i, k, i1, i2, i3, t)
     Else If (iconfg==4) Then
        Call dchin3(i, k, i1, i2, i3, t)
     End If
  End If
  If (icsta(i)/10==11) Then
     ot(next(i)) = otmp
     ct(next(i)) = ctmp
     next(next(i)) = i
     Call wallc(i, i1, i2, i3, t0, tmin1)
     If (tmin1<ct(i)) Then
        icsta(i) = icsta(i) + 10
        t0 = tmin1
        Goto 1000
     End If
  End If
  Return
End Subroutine cellre
