Subroutine getnst(srt)
  Parameter (pimass=0.140, pi=3.1415926)
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  s = srt**2
  nstate = 0
  wtot = 0.
  If (srt<=thresh(1)) Return
  Do i = 1, 15
    weight(i) = 0.
    If (srt>thresh(i)) nstate = i
  End Do
  Do i = 1, nstate
    pf2 = (s-(ppbm(i,1)+ppbm(i,2))**2)*(s-(ppbm(i,1)-ppbm(i,2))**2)/4/s
    weight(i) = pf2*niso(i)
    wtot = wtot + weight(i)
  End Do
  ene = (srt/pimass)**3/(6.*pi**2)
  fsum = factr2(2) + factr2(3)*ene + factr2(4)*ene**2 + factr2(5)*ene**3 + factr2(6)*ene**4
  Return
End Subroutine getnst
