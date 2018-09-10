Subroutine fixtim(l, t, tmin1, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  k = nc
  If (tmin<tmin1) Then
     ot(l) = tmin
     If (ct(l)<tmin1) Then
        icsta(l) = 0
     Else
        icsta(l) = icsta(l) + 10
     End If
     next(l) = k
  Else If (tmin==tmin1) Then
     ot(l) = tmin
     If (nc==0) Then
        next(l) = 0
     Else
        icsta(l) = icsta(l) + 10
        next(l) = k
     End If
  Else
     ot(l) = tmin1
     next(l) = 0
  End If
  Return
End Subroutine fixtim
