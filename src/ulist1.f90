Subroutine ulist1(l, t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  icels0 = icels(l)
  i1 = icels0/10000
  i2 = (icels0-i1*10000)/100
  i3 = icels0 - i1*10000 - i2*100
  k = mod(icsta(l), 10)
  Call wallc(l, i1, i2, i3, t, tmin1)
  tmin = tmin1
  nc = 0
  If (i1==11 .And. i2==11 .And. i3==11) Then
     Call chkout(l, t, tmin, nc)
  Else
     If (iconfg==1) Then
        Call chkin1(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==2) Then
        Call chkin2(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==4) Then
        Call chkin3(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==3 .Or. iconfg==5) Then
        Call chkcel(l, i1, i2, i3, t, tmin, nc)
     End If
  End If
  Call fixtim(l, t, tmin1, tmin, nc)
  Return
End Subroutine ulist1
