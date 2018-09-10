Subroutine wallc(i, i1, i2, i3, t, tmin)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  tmin = tlarge
  If (iconfg<=2 .Or. iconfg==4) Then
     If ((i1>=1 .And. i1<=10) .Or. (i2>=1 .And. i2<=10) .Or. (i3>=1 .And. i3<=10)) Then
        Call wallc1(i, i1, i2, i3, t, tmin)
     Else
        Call wallcb(i, t, tmin)
     End If
  Else If (iconfg==3 .Or. iconfg==5) Then
     Call wallc2(i, i1, i2, i3, t, tmin)
  End If
  Return
End Subroutine wallc
