Subroutine ulist(t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Save
  If (ictype==1 .Or. ictype==2 .Or. ictype==5 .Or. ictype==6) Then
     l = ifmpt
     Call ulist1(l, t)
  End If
  If (ictype/=1) Then
     l = iscat
     Call ulist1(l, t)
     If (jscat/=0) Then
        l = jscat
        Call ulist1(l, t)
     End If
  End If
  Return
End Subroutine ulist
