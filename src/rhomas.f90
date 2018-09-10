Real Function rhomas(dmax, iseed)
  Common /rndf77/nseed
  Save
  dmin = 0.28
  If (dmax<0.77) Then
     fm = frho(dmax)
  Else
     fm = 1.
  End If
  If (fm==0.) fm = 1.E-06
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  If ((ranart(nseed)>frho(dm)/fm) .And. (ntry1<=10)) Goto 10
  If (dm>1.07) Goto 10
  rhomas = dm
  Return
End Function rhomas
