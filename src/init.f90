Subroutine init(minnum, maxnum, num, radius, x0, z0, p0, gamma, iseed, mass, iopt)
  Parameter (maxstr=150001, amu=0.9383)
  Parameter (maxx=20, maxz=24)
  Parameter (pi=3.1415926)
  Real ptot(3)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ss/inout(20)
  Common /rndf77/nseed
  Save
  If (p0/=0.) Then
     sign = p0/abs(p0)
  Else
     sign = 0.
  End If
  scheck = gamma**2 - 1.
  If (scheck<0) Then
     Write (99, *) 'scheck10: ', scheck
     scheck = 0.
  End If
  beta = sign*sqrt(scheck)/gamma
  If (minnum==1) Then
     idnum = 1
  Else
     idnum = -1
  End If
  Do irun = 1, num
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
        id(i) = idnum
        e(i) = amu
     End Do
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
200     Continue
        x = 1.0 - 2.0*ranart(nseed)
        y = 1.0 - 2.0*ranart(nseed)
        z = 1.0 - 2.0*ranart(nseed)
        If ((x*x+y*y+z*z)>1.0) Goto 200
        r(1, i) = x*radius
        r(2, i) = y*radius
        r(3, i) = z*radius
     End Do
  End Do
  If (iopt/=3) Then
     rhow0 = 0.168
     Do irun = 1, num
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
500        Continue
           px = 1.0 - 2.0*ranart(nseed)
           py = 1.0 - 2.0*ranart(nseed)
           pz = 1.0 - 2.0*ranart(nseed)
           If (px*px+py*py+pz*pz>1.0) Goto 500
           rdist = sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
           rhows = rhow0/(1.0+exp((rdist-radius)/0.55))
           pfermi = 0.197*(1.5*pi*pi*rhows)**(1./3.)
           If (iopt==2) pfermi = 0.27
           If (iopt==4) pfermi = 0.
           p(1, i) = pfermi*px
           p(2, i) = pfermi*py
           p(3, i) = pfermi*pz
        End Do
        Do idir = 1, 3
           ptot(idir) = 0.0
        End Do
        npart = 0
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           npart = npart + 1
           Do idir = 1, 3
              ptot(idir) = ptot(idir) + p(idir, i)
           End Do
        End Do
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           Do idir = 1, 3
              p(idir, i) = p(idir, i) - ptot(idir)/float(npart)
           End Do
           If ((iopt==1) .Or. (iopt==2)) Then
              epart = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+amu**2)
              p(3, i) = gamma*(p(3,i)+beta*epart)
           Else
              p(3, i) = p(3, i) + p0
           End If
        End Do
     End Do
  Else
     Do irun = 1, num
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           p(1, i) = 0.0
           p(2, i) = 0.0
           p(3, i) = p0
        End Do
     End Do
  End If
  Do irun = 1, num
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
        r(1, i) = r(1, i) + x0
        r(3, i) = (r(3,i)+z0)/gamma
     End Do
  End Do
  Return
End Subroutine init
