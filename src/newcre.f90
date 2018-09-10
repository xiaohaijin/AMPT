Subroutine newcre(i, k)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Save
  If (k==0) Then
     k = i
     nic(i) = 0
  Else If (nic(k)==0) Then
     nic(k) = i
     nic(i) = k
  Else
     j = k
     Do While (nic(j)/=k)
        j = nic(j)
     End Do
     nic(j) = i
     nic(i) = k
  End If
  Return
End Subroutine newcre
