Subroutine arini1
  Parameter (maxstr=150001)
  Double Precision smearp, smearh
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /smearz/smearp, smearh
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /anim/nevent, isoft, isflag, izpc
  Common /nzpc/nattzp
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /rndf77/nseed
  Common /para8/idpert, npertd, idxsec
  Save
  Open (91, File='../data/deuteron_processes.dat', Status='UNKNOWN')
  If (idpert==1 .Or. idpert==2) Then
    Open (90, File='../data/ampt_pert.dat', Status='UNKNOWN')
  End If
  tau0 = arpar1(1)
  np = iaint2(1)
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    If (np>nattzp) Then
      Do i = nattzp + 1, np
        If ((xmar(i)**2+pxar(i)**2+pyar(i)**2)>0.) Then
          rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
        Else
          Print *, ' IN ARINI1 mt=0'
          rap = 1000000.0*sign(1., pzar(i))
        End If
        vx = pxar(i)/pear(i)
        vy = pyar(i)/pear(i)
        ftar(i) = tau0*cosh(rap)
        gxar(i) = gxar(i) + vx*ftar(i)
        gyar(i) = gyar(i) + vy*ftar(i)
        gzar(i) = tau0*sinh(rap)
        If (pxar(i)==0 .And. pyar(i)==0 .And. (itypar(i)==2112 .Or. itypar(i)==2212)) Then
          If ((pear(i)/hint1(6)>0.99 .And. pear(i)/hint1(6)<1.01) .Or. (pear(i)/hint1(7)>0.99 .And. pear(i)/hint1(7)<1.01)) Then
            taui = 1.E-20
            ftar(i) = taui*cosh(rap)
            gzar(i) = taui*sinh(rap)
          End If
        End If
      End Do
    End If
  Else
    Do i = 1, np
      If ((xmar(i)**2+pxar(i)**2+pyar(i)**2)>0.) Then
        rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
      Else
        Print *, ' IN ARINI1 mt=0'
        rap = 1000000.0*sign(1., pzar(i))
      End If
      vx = pxar(i)/pear(i)
      vy = pyar(i)/pear(i)
      taui = ftar(i) + tau0
      ftar(i) = taui*cosh(rap)
      gxar(i) = gxar(i) + vx*tau0*cosh(rap)
      gyar(i) = gyar(i) + vy*tau0*cosh(rap)
      gzar(i) = taui*sinh(rap)
      zsmear = sngl(smearh)*(2.*ranart(nseed)-1.)
      gzar(i) = gzar(i) + zsmear
      If (pxar(i)==0 .And. pyar(i)==0 .And. (itypar(i)==2112 .Or. itypar(i)==2212)) Then
        If ((pear(i)/hint1(6)>0.99 .And. pear(i)/hint1(6)<1.01) .Or. (pear(i)/hint1(7)>0.99 .And. pear(i)/hint1(7)<1.01)) Then
          taui = 1.E-20
          ftar(i) = taui*cosh(rap)
          gzar(i) = taui*sinh(rap) + zsmear
        End If
      End If
    End Do
  End If
  Call addhad
  Return
End Subroutine arini1
