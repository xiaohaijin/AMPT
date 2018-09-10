Subroutine reor(t, tmin, j, last0)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Save
  icels0 = icels(j)
  i1 = icels0/10000
  i2 = (icels0-i1*10000)/100
  i3 = icels0 - i1*10000 - i2*100
  Call wallc(j, i1, i2, i3, t, tmin1)
  If (tmin<=tmin1) Then
     nc = last0
  Else
     tmin = tmin1
     nc = 0
  End If
  If (iconfg==3 .Or. iconfg==5) Then
     Call chcell(j, i1, i2, i3, last0, t, tmin, nc)
  Else
     If (i1==11 .And. i2==11 .And. i3==11) Then
        Call chout(j, last0, t, tmin, nc)
     Else
        If (iconfg==1) Then
           Call chin1(j, i1, i2, i3, last0, t, tmin, nc)
        Else If (iconfg==2) Then
           Call chin2(j, i1, i2, i3, last0, t, tmin, nc)
        Else If (iconfg==4) Then
           Call chin3(j, i1, i2, i3, last0, t, tmin, nc)
        End If
     End If
  End If
  Call fixtim(j, t, tmin1, tmin, nc)
  Return
End Subroutine reor
