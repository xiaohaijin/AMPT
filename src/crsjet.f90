Subroutine crsjet
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
  Common /njet/n, ipcrs
  Common /bveg1/xl(10), xu(10), acc, ndim, ncall, itmx, nprn
  Common /bveg2/xi(50, 10), si, si2, swgt, schi, ndo, it
  Common /bveg3/f, ti, tsi
  Common /sedvax/num1
  External fjet, fjetrg
  Save
  ndim = 3
  ipcrs = 0
  Call vegas(fjet, avgi, sd, chi2a)
  hint1(14) = sngl(avgi)/2.5682
  If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
    ipcrs = 1
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(15) = sngl(avgi)/2.5682
  End If
  If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
    ipcrs = 2
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(16) = sngl(avgi)/2.5682
  End If
  If (ihpr2(6)==1 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
    ipcrs = 3
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(17) = sngl(avgi)/2.5682
  End If
  If (ihpr2(3)/=0) Then
    ipcrs = 0
    Call vegas(fjetrg, avgi, sd, chi2a)
    hint1(61) = sngl(avgi)/2.5682
    If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
      ipcrs = 1
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(62) = sngl(avgi)/2.5682
    End If
    If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
      ipcrs = 2
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(63) = sngl(avgi)/2.5682
    End If
    If (ihpr2(6)==1 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
      ipcrs = 3
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(64) = sngl(avgi)/2.5682
    End If
  End If
  Return
End Subroutine crsjet
