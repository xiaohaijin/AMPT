Subroutine wallcb(i, t, tmin)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  If (size1==0D0 .Or. size2==0D0 .Or. size3==0D0) Return
  x1p = gx(i)
  x2p = gy(i)
  x3p = gz(i)
  v1p = vx(i)
  v2p = vy(i)
  v3p = vz(i)
  tf = ft(i)
  If (t<size .And. tf<size) Then
     If (x1p<-5D0*size1 .And. v1p>0D0) Then
        t1 = (-5D0*size1-x1p)/v1p + tf
     Else If (x1p>5D0*size1 .And. v1p<0D0) Then
        t1 = -(x1p-5D0*size1)/v1p + tf
     Else
        t1 = tlarge
     End If
     If (t1/=tlarge) Then
        x2pp = x2p + v2p*(t1-tf)
        x3pp = x3p + v3p*(t1-tf)
        If (x2pp<=-5D0*size2 .Or. x2pp>=5D0*size2 .Or. x3pp<=-5D0*size3 .Or. x3pp>=5D0*size3) t1 = tlarge
     End If
     If (x2p<-5D0*size2 .And. v2p>0D0) Then
        t2 = (-5D0*size2-x2p)/v2p + tf
     Else If (x2p>5D0*size2 .And. v2p<0D0) Then
        t2 = -(x2p-5D0*size2)/v2p + tf
     Else
        t2 = tlarge
     End If
     If (t2/=tlarge) Then
        x1pp = x1p + v1p*(t2-tf)
        x3pp = x3p + v3p*(t2-tf)
        If (x1pp<=-5D0*size1 .Or. x1pp>=5D0*size1 .Or. x3pp<=-5D0*size3 .Or. x3pp>=5D0*size3) t2 = tlarge
     End If
     If (x3p<-5D0*size3 .And. v3p>0D0) Then
        t3 = (-5D0*size3-x3p)/v3p + tf
     Else If (x3p>5D0*size3 .And. v3p<0D0) Then
        t3 = -(x3p-5D0*size3)/v3p + tf
     Else
        t3 = tlarge
     End If
     If (t3/=tlarge) Then
        x1pp = x1p + v1p*(t3-tf)
        x2pp = x2p + v2p*(t3-tf)
        If (x1pp<=-5D0*size1 .Or. x1pp>=5D0*size1 .Or. x2pp<=-5D0*size2 .Or. x2pp>=5D0*size2) t3 = tlarge
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
  x1q = x1p + v1p*(t-tf)
  x2q = x2p + v2p*(t-tf)
  x3q = x3p + v3p*(t-tf)
  If (x1q<-5D0*(size1+v1*(t-size)) .And. v1p>-5D0*v1) Then
     t1 = (-5D0*(size1-v1*size)+v1p*tf-x1p)/(v1p-(-5D0)*v1)
     icsta1 = 101
  Else If (x1q>5D0*(size1+v1*(t-size)) .And. v1p<5D0*v1) Then
     t1 = (5D0*(size1-v1*size)+v1p*tf-x1p)/(v1p-5D0*v1)
     icsta1 = 102
  Else
     t1 = tlarge
  End If
  If (t1/=tlarge) Then
     x2pp = x2p + v2p*(t1-tf)
     x3pp = x3p + v3p*(t1-tf)
     If (x2pp<=-5D0*(size2+v2*(t1-size)) .Or. x2pp>=5D0*(size2+v2*(t1-size)) .Or. x3pp<=-5D0*(size3+v3*(t1-size)) .Or. x3pp>=5D0*(size3+v3*(t1-size))) t1 = tlarge
  End If
  If (x2q<-5D0*(size2+v2*(t-size)) .And. v2p>-5D0*v2) Then
     t2 = (-5D0*(size2-v2*size)+v2p*tf-x2p)/(v2p-(-5D0)*v2)
     icsta2 = 103
  Else If (x2q>5D0*(size2+v2*(t-size)) .And. v2p<5D0*v2) Then
     t2 = (5D0*(size2-v2*size)+v2p*tf-x2p)/(v2p-5D0*v2)
     icsta2 = 104
  Else
     t2 = tlarge
  End If
  If (t2/=tlarge) Then
     x1pp = x1p + v1p*(t2-tf)
     x3pp = x3p + v3p*(t2-tf)
     If (x1pp<=-5D0*(size1+v1*(t2-size)) .Or. x1pp>=5D0*(size1+v1*(t2-size)) .Or. x3pp<=-5D0*(size3+v3*(t2-size)) .Or. x3pp>=5D0*(size3+v3*(t2-size))) t2 = tlarge
  End If
  If (x3q<-5D0*(size3+v3*(t-size)) .And. v3p>-5D0*v3) Then
     t3 = (-5D0*(size3-v3*size)+v3p*tf-x3p)/(v3p-(-5D0)*v3)
     icsta3 = 105
  Else If (x3q>5D0*(size3+v3*(t-size)) .And. v3p<5D0*v3) Then
     t3 = (5D0*(size3-v3*size)+v3p*tf-x3p)/(v3p-5D0*v3)
     icsta3 = 106
  Else
     t3 = tlarge
  End If
  If (t3/=tlarge) Then
     x2pp = x2p + v2p*(t3-tf)
     x1pp = x1p + v1p*(t3-tf)
     If (x2pp<=-5D0*(size2+v2*(t3-size)) .Or. x2pp>=5D0*(size2+v2*(t3-size)) .Or. x1pp<=-5D0*(size1+v1*(t3-size)) .Or. x1pp>=5D0*(size1+v1*(t3-size))) t3 = tlarge
  End If
  tmin = min(t1, t2, t3)
  If (tmin==t1) Then
     icsta(i) = icsta1
  Else If (tmin==t2) Then
     icsta(i) = icsta2
  Else If (tmin==t3) Then
     icsta(i) = icsta3
  End If
  Return
End Subroutine wallcb
