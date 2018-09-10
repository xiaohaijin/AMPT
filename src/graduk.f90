Subroutine graduk(ix, iy, iz, gradxk, gradyk, gradzk)
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ss/inout(20)
  Save
  rxplus = rho(ix+1, iy, iz)
  rxmins = rho(ix-1, iy, iz)
  ryplus = rho(ix, iy+1, iz)
  rymins = rho(ix, iy-1, iz)
  rzplus = rho(ix, iy, iz+1)
  rzmins = rho(ix, iy, iz-1)
  gradxk = (rxplus-rxmins)/2.
  gradyk = (ryplus-rymins)/2.
  gradzk = (rzplus-rzmins)/2.
  Return
End Subroutine graduk
