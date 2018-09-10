Subroutine celasn
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist2/icell, icel(10, 10, 10)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Save
  External integ
  i = ifmpt
  tt = ft(i)
  td = tt - size
  If (iconfg==1 .And. (size1==0D0 .Or. size2==0D0 .Or. size3==0D0)) Then
     i1 = 11
     i2 = 11
     i3 = 11
  Else If (iconfg==4 .Or. td<=0D0) Then
     i1 = integ(gx(i)/size1) + 6
     i2 = integ(gy(i)/size2) + 6
     i3 = integ(gz(i)/size3) + 6
     If (integ(gx(i)/size1)==gx(i)/size1 .And. vx(i)<0D0) i1 = i1 - 1
     If (integ(gy(i)/size2)==gy(i)/size2 .And. vy(i)<0D0) i2 = i2 - 1
     If (integ(gz(i)/size3)==gz(i)/size3 .And. vz(i)<0D0) i3 = i3 - 1
  Else
     i1 = integ(gx(i)/(size1+v1*td)) + 6
     i2 = integ(gy(i)/(size2+v2*td)) + 6
     i3 = integ(gz(i)/(size3+v3*td)) + 6
     If (integ(gx(i)/(size1+v1*td))==gx(i)/(size1+v1*td) .And. vx(i)<(i1-6)*v1) i1 = i1 - 1
     If (integ(gy(i)/(size2+v2*td))==gy(i)/(size2+v2*td) .And. vy(i)<(i2-6)*v2) i2 = i2 - 1
     If (integ(gz(i)/(size3+v3*td))==gz(i)/(size3+v3*td) .And. vz(i)<(i3-6)*v3) i3 = i3 - 1
  End If
  If (i1<=0 .Or. i1>=11 .Or. i2<=0 .Or. i2>=11 .Or. i3<=0 .Or. i3>=11) Then
     i1 = 11
     i2 = 11
     i3 = 11
  End If
  If (i1==11) Then
     j = icell
     Call newcre(i, j)
     icell = j
     icels(i) = 111111
  Else
     j = icel(i1, i2, i3)
     Call newcre(i, j)
     icel(i1, i2, i3) = j
     icels(i) = i1*10000 + i2*100 + i3
  End If
  Return
End Subroutine celasn
