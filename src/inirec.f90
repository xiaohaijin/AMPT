Subroutine inirec
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Common /para5/iconfg, iordsc
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  Common /smearz/smearp, smearh
  Common /precpb/vxp(maxptn), vyp(maxptn), vzp(maxptn)
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  Common /anim/nevent, isoft, isflag, izpc
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  Common /rndm3/iseedp
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /arevt/iaevt, iarun, miss
  Save
  iseed = iseedp
  If (isoft==5) Then
     itlast = 0
     Call inifrz
  End If
  Do i = 1, mul
     indxi = indx(i)
     ityp(i) = ityp0(indxi)
     gx(i) = gx0(indxi)
     gy(i) = gy0(indxi)
     gz(i) = gz0(indxi)
     ft(i) = ft0(indxi)
     px(i) = px0(indxi)
     py(i) = py0(indxi)
     pz(i) = pz0(indxi)
     e(i) = e0(indxi)
     xmass(i) = xmass0(indxi)
     lstrg1(i) = lstrg0(indxi)
     lpart1(i) = lpart0(indxi)
     vxp(i) = vxp0(indxi)
     vyp(i) = vyp0(indxi)
     vzp(i) = vzp0(indxi)
     xstrg0(i) = xstrg(indxi)
     ystrg0(i) = ystrg(indxi)
     istrg0(i) = istrg(indxi)
     If (isoft==5) Then
        idfrz(i) = ityp(i)
        gxfrz(i) = gx(i)
        gyfrz(i) = gy(i)
        gzfrz(i) = gz(i)
        ftfrz(i) = ft(i)
        pxfrz(i) = px(i)
        pyfrz(i) = py(i)
        pzfrz(i) = pz(i)
        efrz(i) = e(i)
        xmfrz(i) = xmass(i)
        ifrz(i) = 0
     End If
  End Do
  Do i = 1, mul
     ityps(i) = ityp(i)
     gxs(i) = gx(i)
     gys(i) = gy(i)
     gzs(i) = gz(i)
     fts(i) = ft(i)
     pxs(i) = px(i)
     pys(i) = py(i)
     pzs(i) = pz(i)
     es(i) = e(i)
     xmasss(i) = xmass(i)
  End Do
  If (isoft==1 .And. (ioscar==2 .Or. ioscar==3)) Write (92, *) iaevt, miss, mul
  Do i = 1, mul
     energy = e(i)
     vx(i) = px(i)/energy
     vy(i) = py(i)/energy
     vz(i) = pz(i)/energy
     If (iftflg==0) Then
        formt = ft(i)
        If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
           gx(i) = gx(i) + vxp(i)*formt
           gy(i) = gy(i) + vyp(i)*formt
           gz(i) = gz(i) + vzp(i)*formt
        Else
           gx(i) = gx(i) + vx(i)*formt
           gy(i) = gy(i) + vy(i)*formt
           gz(i) = gz(i) + vz(i)*formt
        End If
     End If
     If (ioscar==2 .Or. ioscar==3) Then
        If (dmax1(abs(gx(i)),abs(gy(i)),abs(gz(i)),abs(ft(i)))<9999) Then
           Write (92, 200) ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i), ft(i), istrg0(i), xstrg0(i), ystrg0(i)
        Else
           Write (92, 201) ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i), ft(i), istrg0(i), xstrg0(i), ystrg0(i)
        End If
     End If
  End Do
  If (iconfg<=3) Then
     Do i = 1, mul
        If (ft(i)<=abs(gz(i))) Then
           eta(i) = 1000000.D0
        Else
           eta(i) = 0.5D0*log((ft(i)+gz(i))/(ft(i)-gz(i)))
        End If
        If (e(i)<=abs(pz(i))) Then
           rap(i) = 1000000.D0
        Else
           rap(i) = 0.5D0*log((e(i)+pz(i))/(e(i)-pz(i)))
        End If
        If (eta(i)<1000000.D0) Then
           tau(i) = ft(i)/cosh(eta(i))
        Else
           tau(i) = 1D-10
        End If
     End Do
     Do i = 1, mul
        etas(i) = eta(i)
        raps(i) = rap(i)
        taus(i) = tau(i)
     End Do
  End If
  Return
200 Format (I3, 2(1X,F7.2), 1X, F8.2, 1X, F6.3, 4(1X,F8.2), 1X, I5, 2(1X,F7.2))
201 Format (I3, 2(1X,F7.2), 1X, F8.2, 1X, F6.3, 4(1X,E8.2), 1X, I5, 2(1X,F7.2))
End Subroutine inirec
