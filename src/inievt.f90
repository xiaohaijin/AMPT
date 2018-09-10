Subroutine inievt
  Implicit Double Precision (A-H, O-Z)
  Common /para1/mul
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Save
  If (ireflg==0) Call readi
  If (igeflg/=0) Call genei
  If (ibstfg/=0) Call boosti
  Return
End Subroutine inievt
