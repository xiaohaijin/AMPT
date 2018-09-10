Subroutine gradup(ix, iy, iz, gradxp, gradyp, gradzp)
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ss/inout(20)
  Save
  rxplus = rhop(ix+1, iy, iz)/rho0
  rxmins = rhop(ix-1, iy, iz)/rho0
  ryplus = rhop(ix, iy+1, iz)/rho0
  rymins = rhop(ix, iy-1, iz)/rho0
  rzplus = rhop(ix, iy, iz+1)/rho0
  rzmins = rhop(ix, iy, iz-1)/rho0
  gradxp = (rxplus-rxmins)/2.
  gradyp = (ryplus-rymins)/2.
  gradzp = (rzplus-rzmins)/2.
  Return
End Subroutine gradup
