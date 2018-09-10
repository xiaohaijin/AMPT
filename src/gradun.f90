Subroutine gradun(ix, iy, iz, gradxn, gradyn, gradzn)
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ss/inout(20)
  Save
  rxplus = rhon(ix+1, iy, iz)/rho0
  rxmins = rhon(ix-1, iy, iz)/rho0
  ryplus = rhon(ix, iy+1, iz)/rho0
  rymins = rhon(ix, iy-1, iz)/rho0
  rzplus = rhon(ix, iy, iz+1)/rho0
  rzmins = rhon(ix, iy, iz-1)/rho0
  gradxn = (rxplus-rxmins)/2.
  gradyn = (ryplus-rymins)/2.
  gradzn = (rzplus-rzmins)/2.
  Return
End Subroutine gradun
