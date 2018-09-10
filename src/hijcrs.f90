Subroutine hijcrs
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /njet/n, ipcrs
  External fhin, ftot, fnjet, ftotjt, ftotrg
  Save
  If (hint1(1)>=10.0) Call crsjet
  aphx1 = hipr1(6)*(ihnt2(1)**0.3333333-1.0)
  aphx2 = hipr1(6)*(ihnt2(3)**0.3333333-1.0)
  hint1(11) = hint1(14) - aphx1*hint1(15) - aphx2*hint1(16) + aphx1*aphx2*hint1(17)
  hint1(10) = gauss1(ftotjt, 0.0, 20.0, 0.01)
  hint1(12) = gauss1(fhin, 0.0, 20.0, 0.01)
  hint1(13) = gauss1(ftot, 0.0, 20.0, 0.01)
  hint1(60) = hint1(61) - aphx1*hint1(62) - aphx2*hint1(63) + aphx1*aphx2*hint1(64)
  hint1(59) = gauss1(ftotrg, 0.0, 20.0, 0.01)
  If (hint1(59)==0.0) hint1(59) = hint1(60)
  If (hint1(1)>=10.0) Then
    Do i = 0, 20
      n = i
      hint1(80+i) = gauss1(fnjet, 0.0, 20.0, 0.01)/hint1(12)
    End Do
  End If
  hint1(10) = hint1(10)*hipr1(31)
  hint1(12) = hint1(12)*hipr1(31)
  hint1(13) = hint1(13)*hipr1(31)
  hint1(59) = hint1(59)*hipr1(31)
  If (ihpr2(13)/=0) Then
    hipr1(33) = 1.36*(1.0+36.0/hint1(1)**2)*alog(0.6+0.1*hint1(1)**2)
    hipr1(33) = hipr1(33)/hint1(12)
  End If
  Return
End Subroutine hijcrs
