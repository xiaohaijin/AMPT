Subroutine dens(ipot, mass, num, nesc)
  Parameter (maxstr=150001, maxr=1)
  Parameter (maxx=20, maxz=24)
  Dimension pxl(-maxx:maxx, -maxx:maxx, -maxz:maxz), pyl(-maxx:maxx, -maxx:maxx, -maxz:maxz), pzl(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ddpi/pirho(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ss/inout(20)
  Common /rr/massr(0:maxr)
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /bbb/bxx(-maxx:maxx, -maxx:maxx, -maxz:maxz), byy(-maxx:maxx, -maxx:maxx, -maxz:maxz), bzz(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./
  Do iz = -maxz, maxz
     Do iy = -maxx, maxx
        Do ix = -maxx, maxx
           rho(ix, iy, iz) = 0.0
           rhon(ix, iy, iz) = 0.0
           rhop(ix, iy, iz) = 0.0
           pirho(ix, iy, iz) = 0.0
           pxl(ix, iy, iz) = 0.0
           pyl(ix, iy, iz) = 0.0
           pzl(ix, iy, iz) = 0.0
           pel(ix, iy, iz) = 0.0
           bxx(ix, iy, iz) = 0.0
           byy(ix, iy, iz) = 0.0
           bzz(ix, iy, iz) = 0.0
        End Do
     End Do
  End Do
  nesc = 0
  big = 1.0/(3.0*float(num))
  small = 1.0/(9.0*float(num))
  msum = 0
  Do irun = 1, num
     msum = msum + massr(irun-1)
     Do j = 1, massr(irun)
        i = j + msum
        ix = nint(r(1,i))
        iy = nint(r(2,i))
        iz = nint(r(3,i))
        If (ix<=-maxx .Or. ix>=maxx .Or. iy<=-maxx .Or. iy>=maxx .Or. iz<=-maxz .Or. iz>=maxz) Then
           nesc = nesc + 1
        Else
           If (j>mass) Goto 30
           rho(ix, iy, iz) = rho(ix, iy, iz) + big
           rho(ix+1, iy, iz) = rho(ix+1, iy, iz) + small
           rho(ix-1, iy, iz) = rho(ix-1, iy, iz) + small
           rho(ix, iy+1, iz) = rho(ix, iy+1, iz) + small
           rho(ix, iy-1, iz) = rho(ix, iy-1, iz) + small
           rho(ix, iy, iz+1) = rho(ix, iy, iz+1) + small
           rho(ix, iy, iz-1) = rho(ix, iy, iz-1) + small
           If (zet(lb(i))/=0) Then
              rhop(ix, iy, iz) = rhop(ix, iy, iz) + big
              rhop(ix+1, iy, iz) = rhop(ix+1, iy, iz) + small
              rhop(ix-1, iy, iz) = rhop(ix-1, iy, iz) + small
              rhop(ix, iy+1, iz) = rhop(ix, iy+1, iz) + small
              rhop(ix, iy-1, iz) = rhop(ix, iy-1, iz) + small
              rhop(ix, iy, iz+1) = rhop(ix, iy, iz+1) + small
              rhop(ix, iy, iz-1) = rhop(ix, iy, iz-1) + small
              Goto 40
           End If
           If (zet(lb(i))==0) Then
              rhon(ix, iy, iz) = rhon(ix, iy, iz) + big
              rhon(ix+1, iy, iz) = rhon(ix+1, iy, iz) + small
              rhon(ix-1, iy, iz) = rhon(ix-1, iy, iz) + small
              rhon(ix, iy+1, iz) = rhon(ix, iy+1, iz) + small
              rhon(ix, iy-1, iz) = rhon(ix, iy-1, iz) + small
              rhon(ix, iy, iz+1) = rhon(ix, iy, iz+1) + small
              rhon(ix, iy, iz-1) = rhon(ix, iy, iz-1) + small
              Goto 40
           End If
30         pirho(ix, iy, iz) = pirho(ix, iy, iz) + big
           pirho(ix+1, iy, iz) = pirho(ix+1, iy, iz) + small
           pirho(ix-1, iy, iz) = pirho(ix-1, iy, iz) + small
           pirho(ix, iy+1, iz) = pirho(ix, iy+1, iz) + small
           pirho(ix, iy-1, iz) = pirho(ix, iy-1, iz) + small
           pirho(ix, iy, iz+1) = pirho(ix, iy, iz+1) + small
           pirho(ix, iy, iz-1) = pirho(ix, iy, iz-1) + small
40         pxl(ix, iy, iz) = pxl(ix, iy, iz) + p(1, i)*big
           pxl(ix+1, iy, iz) = pxl(ix+1, iy, iz) + p(1, i)*small
           pxl(ix-1, iy, iz) = pxl(ix-1, iy, iz) + p(1, i)*small
           pxl(ix, iy+1, iz) = pxl(ix, iy+1, iz) + p(1, i)*small
           pxl(ix, iy-1, iz) = pxl(ix, iy-1, iz) + p(1, i)*small
           pxl(ix, iy, iz+1) = pxl(ix, iy, iz+1) + p(1, i)*small
           pxl(ix, iy, iz-1) = pxl(ix, iy, iz-1) + p(1, i)*small
           pyl(ix, iy, iz) = pyl(ix, iy, iz) + p(2, i)*big
           pyl(ix+1, iy, iz) = pyl(ix+1, iy, iz) + p(2, i)*small
           pyl(ix-1, iy, iz) = pyl(ix-1, iy, iz) + p(2, i)*small
           pyl(ix, iy+1, iz) = pyl(ix, iy+1, iz) + p(2, i)*small
           pyl(ix, iy-1, iz) = pyl(ix, iy-1, iz) + p(2, i)*small
           pyl(ix, iy, iz+1) = pyl(ix, iy, iz+1) + p(2, i)*small
           pyl(ix, iy, iz-1) = pyl(ix, iy, iz-1) + p(2, i)*small
           pzl(ix, iy, iz) = pzl(ix, iy, iz) + p(3, i)*big
           pzl(ix+1, iy, iz) = pzl(ix+1, iy, iz) + p(3, i)*small
           pzl(ix-1, iy, iz) = pzl(ix-1, iy, iz) + p(3, i)*small
           pzl(ix, iy+1, iz) = pzl(ix, iy+1, iz) + p(3, i)*small
           pzl(ix, iy-1, iz) = pzl(ix, iy-1, iz) + p(3, i)*small
           pzl(ix, iy, iz+1) = pzl(ix, iy, iz+1) + p(3, i)*small
           pzl(ix, iy, iz-1) = pzl(ix, iy, iz-1) + p(3, i)*small
           pel(ix, iy, iz) = pel(ix, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*big
           pel(ix+1, iy, iz) = pel(ix+1, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix-1, iy, iz) = pel(ix-1, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy+1, iz) = pel(ix, iy+1, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy-1, iz) = pel(ix, iy-1, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy, iz+1) = pel(ix, iy, iz+1) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy, iz-1) = pel(ix, iy, iz-1) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
        End If
     End Do
  End Do
  Do iz = -maxz, maxz
     Do iy = -maxx, maxx
        Do ix = -maxx, maxx
           If ((rho(ix,iy,iz)==0) .Or. (pel(ix,iy,iz)==0)) Goto 101
           smass2 = pel(ix, iy, iz)**2 - pxl(ix, iy, iz)**2 - pyl(ix, iy, iz)**2 - pzl(ix, iy, iz)**2
           If (smass2<=0) smass2 = 1.E-06
           smass = sqrt(smass2)
           If (smass==0.) smass = 1.E-06
           gamma = pel(ix, iy, iz)/smass
           If (gamma==0) Goto 101
           bxx(ix, iy, iz) = pxl(ix, iy, iz)/pel(ix, iy, iz)
           byy(ix, iy, iz) = pyl(ix, iy, iz)/pel(ix, iy, iz)
           bzz(ix, iy, iz) = pzl(ix, iy, iz)/pel(ix, iy, iz)
           rho(ix, iy, iz) = rho(ix, iy, iz)/gamma
           rhon(ix, iy, iz) = rhon(ix, iy, iz)/gamma
           rhop(ix, iy, iz) = rhop(ix, iy, iz)/gamma
           pirho(ix, iy, iz) = pirho(ix, iy, iz)/gamma
           pel(ix, iy, iz) = pel(ix, iy, iz)/(gamma**2)
           rho0 = 0.163
           If (ipot==0) Then
              u = 0
              Goto 70
           End If
           If (ipot==1 .Or. ipot==6) Then
              a = -0.1236
              b = 0.0704
              s = 2
              Goto 60
           End If
           If (ipot==2 .Or. ipot==7) Then
              a = -0.218
              b = 0.164
              s = 4./3.
           End If
           If (ipot==3) Then
              a = -0.3581
              b = 0.3048
              s = 1.167
              Goto 60
           End If
           If (ipot==4) Then
              denr = rho(ix, iy, iz)/rho0
              b = 0.3048
              s = 1.167
              If (denr<=4 .Or. denr>7) Then
                 a = -0.3581
              Else
                 a = -b*denr**(1./6.) - 2.*0.036/3.*denr**(-0.333)
              End If
              Goto 60
           End If
60         u = 0.5*a*rho(ix, iy, iz)**2/rho0 + b/(1+s)*(rho(ix,iy,iz)/rho0)**s*rho(ix, iy, iz)
70         pel(ix, iy, iz) = pel(ix, iy, iz) + u
101     End Do
     End Do
  End Do
  Return
End Subroutine dens
