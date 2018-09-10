Subroutine zpca2a
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /para5/iconfg, iordsc
  Common /para6/centy
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /ilist6/t, iopern, icolln
  Common /rndm1/number
  Common /rndm2/iff
  Common /rndm3/iseedp
  Common /ana1/ts(12)
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
  Save
  Do i = 1, ichkpt
     rapi = rap(i)
     et = dsqrt(px(i)**2+py(i)**2+xmp**2)
     Do j = 1, 24
        If (rapi>j+centy-13D0 .And. rapi<j+centy-12D0) Then
           fdetdy(j) = fdetdy(j) + et
           fdndy(j) = fdndy(j) + 1D0
        End If
     End Do
     Do j = 1, 12
        If (et>0.5D0*(j-1) .And. et<0.5D0*j) Then
           fdndpt(j) = fdndpt(j) + 1D0
        End If
     End Do
     If (iconfg==1) Then
        t1 = ft(i)
        t2 = tlarge
        ipic = 11
     Else
        t1 = tau(i)
        t2 = tlarge
        ipic = 12
     End If
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           Call zpca1b(rapi, et, ian)
        End If
     End Do
     If (iconfg==1) Then
        Call zpca1b(rapi, et, 12)
     End If
  End Do
  Do ian = 1, 12
     If (dn(ian)==0D0 .Or. dn1(ian)==0D0 .Or. dn2(ian)==0D0) Then
     End If
     detdy(ian) = detdy(ian) + det(ian)
     If (dn(ian)/=0) Then
        detdn(ian) = detdn(ian) + det(ian)/dn(ian)
     End If
     dndy(ian) = dndy(ian) + dn(ian)
     detdy1(ian) = detdy1(ian) + det1(ian)
     If (dn1(ian)/=0) Then
        detdn1(ian) = detdn1(ian) + det1(ian)/dn1(ian)
     End If
     dndy1(ian) = dndy1(ian) + dn1(ian)
     detdy2(ian) = detdy2(ian) + det2(ian)
     If (dn2(ian)/=0) Then
        detdn2(ian) = detdn2(ian) + det2(ian)/dn2(ian)
     End If
     dndy2(ian) = dndy2(ian) + dn2(ian)
  End Do
  Return
End Subroutine zpca2a
