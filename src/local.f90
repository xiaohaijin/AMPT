Subroutine local(t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (r0=1D0)
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /coal/dpcoal, drcoal, ecritl
  Save
  Do it = 1, 301
     If (t>=tfrz(it) .And. t<tfrz(it+1)) Then
        If (it==itlast) Then
           Return
        Else
           itlast = it
           Goto 50
        End If
     End If
  End Do
  Write (1, *) 'local time out of range in LOCAL, stop', t, it
  Stop
50 Continue
  Do ip = 1, mul
     If (ifrz(ip)==1) Goto 200
     If (it==301) Then
        etcrit = 1D6
        Goto 150
     Else
        etcrit = (ecritl*2D0/3D0)
     End If
     If (t<ft5(ip)) Goto 200
     rap0 = rap(ip)
     eta0 = eta(ip)
     x0 = gx5(ip) + vx(ip)*(t-ft5(ip))
     y0 = gy5(ip) + vy(ip)*(t-ft5(ip))
     detdy = 0D0
     Do itest = 1, mul
        If (itest==ip .Or. t<ft5(itest)) Goto 100
        ettest = eta(itest)
        xtest = gx5(itest) + vx(itest)*(t-ft5(itest))
        ytest = gy5(itest) + vy(itest)*(t-ft5(itest))
        drt = sqrt((xtest-x0)**2+(ytest-y0)**2)
        If (dabs(ettest-eta0)<=1D0 .And. drt<=r0) detdy = detdy + dsqrt(px5(itest)**2+py5(itest)**2+xmass5(itest)**2)*0.5D0
100  End Do
     detdy = detdy*(dcosh(eta0)**2)/(t*3.1416D0*r0**2*dcosh(rap0))
150  If (detdy<=etcrit) Then
        ifrz(ip) = 1
        idfrz(ip) = ityp5(ip)
        pxfrz(ip) = px5(ip)
        pyfrz(ip) = py5(ip)
        pzfrz(ip) = pz5(ip)
        efrz(ip) = e5(ip)
        xmfrz(ip) = xmass5(ip)
        If (t>ft5(ip)) Then
           gxfrz(ip) = x0
           gyfrz(ip) = y0
           gzfrz(ip) = gz5(ip) + vz(ip)*(t-ft5(ip))
           ftfrz(ip) = t
        Else
           gxfrz(ip) = gx5(ip)
           gyfrz(ip) = gy5(ip)
           gzfrz(ip) = gz5(ip)
           ftfrz(ip) = ft5(ip)
        End If
     End If
200 End Do
  Return
End Subroutine local
