Subroutine oldcre(i)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Save
  If (nic(i)==0) Return
  j = nic(i)
  If (nic(j)==i) Then
     nic(j) = 0
     Return
  End If
  Do While (nic(j)/=i)
     j = nic(j)
  End Do
  nic(j) = nic(i)
  Return
End Subroutine oldcre
