Subroutine wallc1(i, i1, i2, i3, t, tmin)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  x1p = gx(i)
  x2p = gy(i)
  x3p = gz(i)
  tf = ft(i)
  v1p = vx(i)
  v2p = vy(i)
  v3p = vz(i)
  If (t<size .And. tf<size) Then
     If (v1p>0D0) Then
        t1 = ((dble(i1)-5D0)*size1-x1p)/v1p + tf
     Else If (v1p<0D0) Then
        t1 = ((dble(i1)-6D0)*size1-x1p)/v1p + tf
     Else
        t1 = tlarge
     End If
     If (v2p>0D0) Then
        t2 = ((dble(i2)-5D0)*size2-x2p)/v2p + tf
     Else If (v2p<0D0) Then
        t2 = ((dble(i2)-6D0)*size2-x2p)/v2p + tf
     Else
        t2 = tlarge
     End If
     If (v3p>0D0) Then
        t3 = ((dble(i3)-5D0)*size3-x3p)/v3p + tf
     Else If (v3p<0D0) Then
        t3 = ((dble(i3)-6D0)*size3-x3p)/v3p + tf
     Else
        t3 = tlarge
     End If
     tmin = min(t1, t2, t3)
     If (tmin==t1) Then
        If (v1p>0D0) Then
           icsta(i) = 101
        Else
           icsta(i) = 102
        End If
     End If
     If (tmin==t2) Then
        If (v2p>0D0) Then
           icsta(i) = 103
        Else
           icsta(i) = 104
        End If
     End If
     If (tmin==t3) Then
        If (v3p>0D0) Then
           icsta(i) = 105
        Else
           icsta(i) = 106
        End If
     End If
     If (tmin<=size) Return
  End If
  If (v1p>(i1-5)*v1) Then
     t1 = ((i1-5)*(size1-v1*size)+v1p*tf-x1p)/(v1p-(i1-5)*v1)
  Else If (v1p<(i1-6)*v1) Then
     t1 = ((i1-6)*(size1-v1*size)+v1p*tf-x1p)/(v1p-(i1-6)*v1)
  Else
     t1 = tlarge
  End If
  If (v2p>(i2-5)*v2) Then
     t2 = ((i2-5)*(size2-v2*size)+v2p*tf-x2p)/(v2p-(i2-5)*v2)
  Else If (v2p<(i2-6)*v2) Then
     t2 = ((i2-6)*(size2-v2*size)+v2p*tf-x2p)/(v2p-(i2-6)*v2)
  Else
     t2 = tlarge
  End If
  If (v3p>(i3-5)*v3) Then
     t3 = ((i3-5)*(size3-v3*size)+v3p*tf-x3p)/(v3p-(i3-5)*v3)
  Else If (v3p<(i3-6)*v3) Then
     t3 = ((i3-6)*(size3-v3*size)+v3p*tf-x3p)/(v3p-(i3-6)*v3)
  Else
     t3 = tlarge
  End If
  tmin = min(t1, t2, t3)
  If (tmin==t1) Then
     If (v1p>(i1-5)*v1) Then
        icsta(i) = 101
     Else
        icsta(i) = 102
     End If
  End If
  If (tmin==t2) Then
     If (v2p>(i2-5)*v2) Then
        icsta(i) = 103
     Else
        icsta(i) = 104
     End If
  End If
  If (tmin==t3) Then
     If (v3p>(i3-5)*v3) Then
        icsta(i) = 105
     Else
        icsta(i) = 106
     End If
  End If
  Return
End Subroutine wallc1
