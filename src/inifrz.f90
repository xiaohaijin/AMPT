Subroutine inifrz
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  Save
  step1 = 0.1D0
  step2 = 1D0
  step3 = 10D0
  step4 = 100D0
  Do it = 1, 101
     tfrz(it) = 0D0 + dble(it-1)*step1
  End Do
  Do it = 102, 191
     tfrz(it) = 10D0 + dble(it-101)*step2
  End Do
  Do it = 192, 281
     tfrz(it) = 100D0 + dble(it-191)*step3
  End Do
  Do it = 282, 301
     tfrz(it) = 1000D0 + dble(it-281)*step4
  End Do
  tfrz(302) = tlarge
  Return
End Subroutine inifrz
