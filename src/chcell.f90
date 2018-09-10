Subroutine chcell(il, i1, i2, i3, last0, t, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist2/icell, icel(10, 10, 10)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Save
  If (iconfg==3 .Or. iconfg==5) Then
     jj = ichkpt
     Do j = 1, jj
        jmintm = j
        If (j/=il .And. j/=last0) Call mintm(il, jmintm, tmin, nc)
     End Do
     Return
  End If
  If (i1==11 .And. i2==11 .And. i3==11) Then
     l = icell
  Else
     l = icel(i1, i2, i3)
  End If
  If (l==0) Return
  j = nic(l)
  If (j==0) Then
     If (l==il .Or. l==last0) Return
     Call mintm(il, l, tmin, nc)
  Else
     If (l/=il .And. l/=last0) Call mintm(il, l, tmin, nc)
     Do While (j/=l)
        If (j/=il .And. j/=last0) Call mintm(il, j, tmin, nc)
        j = nic(j)
     End Do
  End If
  Return
End Subroutine chcell
