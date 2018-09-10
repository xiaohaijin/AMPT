Subroutine zpca1
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Save
  If (mod(ictype,2)==0) Then
     Call zpca1a(iscat)
     Call zpca1a(jscat)
  End If
  Return
End Subroutine zpca1
