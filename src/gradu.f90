Subroutine gradu(iopt, ix, iy, iz, gradx, grady, gradz)
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.167)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ss/inout(20)
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Save
  rxplus = rho(ix+1, iy, iz)/rho0
  rxmins = rho(ix-1, iy, iz)/rho0
  ryplus = rho(ix, iy+1, iz)/rho0
  rymins = rho(ix, iy-1, iz)/rho0
  rzplus = rho(ix, iy, iz+1)/rho0
  rzmins = rho(ix, iy, iz-1)/rho0
  den0 = rho(ix, iy, iz)/rho0
  ene0 = pel(ix, iy, iz)
  Goto (1, 2, 3, 4, 5) iopt
  If (iopt==6) Goto 6
  If (iopt==7) Goto 7
1 Continue
  gradx = -0.062*(rxplus-rxmins) + 0.03525*(rxplus**2-rxmins**2)
  grady = -0.062*(ryplus-rymins) + 0.03525*(ryplus**2-rymins**2)
  gradz = -0.062*(rzplus-rzmins) + 0.03525*(rzplus**2-rzmins**2)
  Return
2 Continue
  expnt = 1.3333333
  gradx = -0.109*(rxplus-rxmins) + 0.082*(rxplus**expnt-rxmins**expnt)
  grady = -0.109*(ryplus-rymins) + 0.082*(ryplus**expnt-rymins**expnt)
  gradz = -0.109*(rzplus-rzmins) + 0.082*(rzplus**expnt-rzmins**expnt)
  Return
3 Continue
  expnt = 1.1666667
  acoef = 0.178
  gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
  grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
  gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
  Return
4 Continue
  eh = 4.
  eqgp = 7.
  acoef = 0.178
  expnt = 1.1666667
  denr = rho(ix, iy, iz)/rho0
  If (denr<=eh .Or. denr>=eqgp) Then
     gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
     grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
     gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
  Else
     acoef1 = 0.178
     acoef2 = 0.0
     expnt2 = 2./3.
     gradx = -acoef1*(rxplus**expnt-rxmins**expnt) - acoef2*(rxplus**expnt2-rxmins**expnt2)
     grady = -acoef1*(ryplus**expnt-rymins**expnt) - acoef2*(ryplus**expnt2-rymins**expnt2)
     gradz = -acoef1*(rzplus**expnt-rzmins**expnt) - acoef2*(rzplus**expnt2-rzmins**expnt2)
  End If
  Return
5 Continue
  expnt = 2.77
  gradx = -0.0516*(rxplus-rxmins) + 0.02498*(rxplus**expnt-rxmins**expnt)
  grady = -0.0516*(ryplus-rymins) + 0.02498*(ryplus**expnt-rymins**expnt)
  gradz = -0.0516*(rzplus-rzmins) + 0.02498*(rzplus**expnt-rzmins**expnt)
  Return
6 Continue
  If (ene0<=0.5) Then
     gradx = -0.062*(rxplus-rxmins) + 0.03525*(rxplus**2-rxmins**2)
     grady = -0.062*(ryplus-rymins) + 0.03525*(ryplus**2-rymins**2)
     gradz = -0.062*(rzplus-rzmins) + 0.03525*(rzplus**2-rzmins**2)
     Return
  End If
  If (ene0>0.5 .And. ene0<=1.5) Then
     ef = 36./1000.
     gradx = -0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = -0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = -0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
  If (ene0>1.5) Then
     ef = 36./1000.
     cf0 = 0.8
     gradx = 0.5*cf0*(rxplus**0.333-rxmins**0.333) - 0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = 0.5*cf0*(ryplus**0.333-rymins**0.333) - 0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = 0.5*cf0*(rzplus**0.333-rzmins**0.333) - 0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
7 Continue
  If (den0<=4.5) Then
     expnt = 1.1666667
     acoef = 0.178
     gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
     grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
     gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
     Return
  End If
  If (den0>4.5 .And. den0<=5.1) Then
     ef = 36./1000.
     gradx = -0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = -0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = -0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
  If (den0>5.1) Then
     ef = 36./1000.
     cf0 = 0.8
     gradx = 0.5*cf0*(rxplus**0.333-rxmins**0.333) - 0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = 0.5*cf0*(ryplus**0.333-rymins**0.333) - 0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = 0.5*cf0*(rzplus**0.333-rzmins**0.333) - 0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
End Subroutine gradu
