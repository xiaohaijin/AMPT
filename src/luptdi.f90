Subroutine luptdi(kfl, px, py)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  kfla = iabs(kfl)
  pt = parj(21)*sqrt(-log(max(1E-10,rlu(0))))
  If (mstj(91)==1) pt = parj(22)*pt
  If (kfla==0 .And. mstj(13)<=0) pt = 0.
  phi = paru(2)*rlu(0)
  px = pt*cos(phi)
  py = pt*sin(phi)
  Return
End Subroutine luptdi
