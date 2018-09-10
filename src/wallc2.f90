Subroutine wallc2(i, i1, i2, i3, t, tmin)
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
  If (v1p>0D0) Then
     t1 = (5D0*size1-x1p)/v1p + tf
  Else If (v1p<0D0) Then
     t1 = (-5D0*size1-x1p)/v1p + tf
  Else
     t1 = tlarge
  End If
  If (v2p>0D0) Then
     t2 = (5D0*size2-x2p)/v2p + tf
  Else If (v2p<0D0) Then
     t2 = (-5D0*size2-x2p)/v2p + tf
  Else
     t2 = tlarge
  End If
  If (iconfg==5) Then
     If (v3p>0D0) Then
        t3 = (5D0*size3-x3p)/v3p + tf
     Else If (v3p<0D0) Then
        t3 = (-5D0*size3-x3p)/v3p + tf
     Else
        t3 = tlarge
     End If
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
  Return
End Subroutine wallc2
