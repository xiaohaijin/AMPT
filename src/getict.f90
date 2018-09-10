Subroutine getict(t1)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  t1 = tlarge
  iscat = 0
  jscat = 0
  Do i = 1, ichkpt
     If (ot(i)<t1) Then
        t1 = ot(i)
        iscat = i
     End If
  End Do
  If (iscat/=0) jscat = next(iscat)
  If (iscat/=0 .And. jscat/=0) Then
     If (icsta(iscat)==0 .And. icsta(jscat)==0) Then
        ictype = 0
     Else
        ictype = 4
     End If
  Else If (iscat/=0 .Or. jscat/=0) Then
     ictype = 3
  End If
  If (ifmpt<=mul) Then
     If (ft(ifmpt)<t1) Then
        ictype = 1
        t1 = ft(ifmpt)
     Else If (ft(ifmpt)==t1) Then
        If (ictype==0) ictype = 2
        If (ictype==3) ictype = 5
        If (ictype==4) ictype = 6
     End If
  End If
  Return
End Subroutine getict
