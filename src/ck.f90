Subroutine ck(l, ick)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Save
  ick = 1
  If (ictype==1) Then
     If (l==ifmpt) ick = 0
  Else If (ictype==0 .Or. ictype==3 .Or. ictype==4) Then
     If (l==iscat .Or. l==jscat) ick = 0
  Else
     If (l==iscat .Or. l==jscat .Or. l==ifmpt) ick = 0
  End If
  Return
End Subroutine ck
