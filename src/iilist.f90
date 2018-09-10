Subroutine iilist
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist2/icell, icel(10, 10, 10)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /ilist6/t, iopern, icolln
  Save
  iscat = maxptn
  jscat = maxptn
  Do i = 1, mul
     next(i) = 0
     last(i) = 0
     icsta(i) = 0
     nic(i) = 0
     icels(i) = 0
  End Do
  icell = 0
  Do i1 = 1, 10
     Do i2 = 1, 10
        Do i3 = 1, 10
           icel(i1, i2, i3) = 0
        End Do
     End Do
  End Do
  ichkpt = 0
  ifmpt = 1
  Do i = 1, mul
     ct(i) = tlarge
     ot(i) = tlarge
  End Do
  iopern = 0
  icolln = 0
  t = 0.D0
  Return
End Subroutine iilist
