Subroutine embedhighpt
  Parameter (maxstr=150001, maxr=1, pichmass=0.140, pi0mass=0.135, pi=3.1415926, nxymax=10001)
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Common /rndf77/nseed
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /anim/nevent, isoft, isflag, izpc
  Common /arevt/iaevt, iarun, miss
  Common /xyembed/nxyjet, xyjet(nxymax, 2)
  Save
  If (iembed==1 .Or. iembed==2) Then
     xjet = xembd
     yjet = yembd
  Else If (iembed==3 .Or. iembed==4) Then
     If (nevent<=nxyjet) Then
        Read (97, *) xjet, yjet
     Else
        ixy = mod(iaevt, nxyjet)
        If (ixy==0) ixy = nxyjet
        xjet = xyjet(ixy, 1)
        yjet = xyjet(ixy, 2)
     End If
  Else
     Return
  End If
  ptq = sqrt(pxqembd**2+pyqembd**2)
  If (ptq<(pichmass/2.)) Then
     Print *, 'Embedded quark transverse momentum is too small'
     Stop
  End If
  idqembd = 1 + int(2*ranart(nseed))
  If (idqembd==1) Then
     idqsoft = -2
     idpi1 = -211
  Else If (idqembd==2) Then
     idqsoft = -1
     idpi1 = 211
  Else
     Print *, 'Wrong quark flavor embedded'
     Stop
  End If
  xmq = ulmass(idqembd)
  xmqsoft = ulmass(idqsoft)
  ptpi = ((pichmass**2+xmq**2-xmqsoft**2)*ptq-sqrt((xmq**2+ptq**2)*(pichmass**4-2.*pichmass**2*(xmq**2+xmqsoft**2)+(xmq**2-xmqsoft**2)**2)))/(2.*xmq**2)
  If (iembed==1 .Or. iembed==3) Then
     pxpi1 = ptpi*pxqembd/ptq
     pypi1 = ptpi*pyqembd/ptq
     phidecomp = acos(pxqembd/ptq)
     If (pyqembd<0) phidecomp = 2.*pi - phidecomp
  Else
     phidecomp = 2.*pi*ranart(nseed)
     pxpi1 = ptpi*cos(phidecomp)
     pypi1 = ptpi*sin(phidecomp)
  End If
  pzpi1 = 0.
  Do ipion = 1, 2
     If (ipion==1) Then
        idpi = idpi1
        pxpi = pxpi1
        pypi = pypi1
        pzpi = pzpi1
     Else If (ipion==2) Then
        idpi = -idpi1
        pxpi = -pxpi1
        pypi = -pypi1
        pzpi = -pzpi1
     End If
     natt = natt + 1
     katt(natt, 1) = idpi
     katt(natt, 2) = 40
     katt(natt, 3) = 0
     patt(natt, 1) = pxpi
     patt(natt, 2) = pypi
     patt(natt, 3) = pzpi
     patt(natt, 4) = sqrt(pxpi**2+pypi**2+pzpi**2+pichmass**2)
     eatt = eatt + patt(natt, 4)
     gxar(natt) = xjet
     gyar(natt) = yjet
     gzar(natt) = 0.
     ftar(natt) = 0.
     itypar(natt) = katt(natt, 1)
     pxar(natt) = patt(natt, 1)
     pyar(natt) = patt(natt, 2)
     pzar(natt) = patt(natt, 3)
     pear(natt) = patt(natt, 4)
     xmar(natt) = pichmass
  End Do
  If (nsembd>0) Then
     Do ipion = 1, 2
        Do ispion = 1, nsembd
           idsart = 3 + int(3*ranart(nseed))
           If (idsart==3) Then
              pimass = pichmass
              idpis = -211
           Else If (idsart==4) Then
              pimass = pi0mass
              idpis = 111
           Else
              pimass = pichmass
              idpis = 211
           End If
           natt = natt + 1
           katt(natt, 1) = idpis
           katt(natt, 2) = 40
           katt(natt, 3) = 0
           theta = tmaxembd*ranart(nseed)
           phi = 2.*pi*ranart(nseed)
           pxspi = psembd*sin(theta)*cos(phi)
           pyspi = psembd*sin(theta)*sin(phi)
           pzspi = psembd*cos(theta)
           If (ipion==1) Then
              Call rotate(pxpi1, pypi1, pzpi1, pxspi, pyspi, pzspi)
           Else
              Call rotate(-pxpi1, -pypi1, -pzpi1, pxspi, pyspi, pzspi)
           End If
           patt(natt, 1) = pxspi
           patt(natt, 2) = pyspi
           patt(natt, 3) = pzspi
           patt(natt, 4) = sqrt(psembd**2+pimass**2)
           eatt = eatt + patt(natt, 4)
           gxar(natt) = xjet
           gyar(natt) = yjet
           gzar(natt) = 0.
           ftar(natt) = 0.
           itypar(natt) = katt(natt, 1)
           pxar(natt) = patt(natt, 1)
           pyar(natt) = patt(natt, 2)
           pzar(natt) = patt(natt, 3)
           pear(natt) = patt(natt, 4)
           xmar(natt) = pimass
        End Do
     End Do
  End If
  Return
End Subroutine embedhighpt
