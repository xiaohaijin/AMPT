Subroutine scat(t, iscat, jscat)
  Implicit Double Precision (A-H, O-Z)
  Save
  Call newpos(t, iscat)
  Call newpos(t, jscat)
  Call newmom(t)
  Return
End Subroutine scat
