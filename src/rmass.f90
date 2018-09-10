Real Function rmass(dmax, iseed)
  Common /rndf77/nseed
  Save
  dmin = 1.078
  If (dmax<1.232) Then
     fm = fdelta(dmax)
  Else
     fm = 1.
  End If
  If (fm==0.) fm = 1.E-06
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  If ((ranart(nseed)>fdelta(dm)/fm) .And. (ntry1<=10)) Goto 10
  If (dm>1.47) Goto 10
  rmass = dm
  Return
End Function rmass
