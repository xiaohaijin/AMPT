Subroutine chkcel(il, i1, i2, i3, t, tmin, nc)
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
        Call ck(j, ick)
        jud2 = j
        If (ick==1) Call ud2(jud2, il, t, tmin, nc)
     End Do
     Return
  End If
  If (i1==11 .And. i2==11 .And. i3==11) Then
     l = icell
  Else
     l = icel(i1, i2, i3)
  End If
  If (l==0) Then
     Return
  End If
  j = nic(l)
  If (j==0) Then
     Call ck(l, ick)
     If (ick==1) Call ud2(l, il, t, tmin, nc)
  Else
     Call ck(l, ick)
     If (ick==1) Call ud2(l, il, t, tmin, nc)
     Do While (j/=l)
        Call ck(j, ick)
        If (ick==1) Call ud2(j, il, t, tmin, nc)
        j = nic(j)
     End Do
  End If
  Return
End Subroutine chkcel
