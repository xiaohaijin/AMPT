Subroutine luboei(nsav)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension dps(4), kfbe(9), nbe(0:9), bei(100)
  Data kfbe/211, -211, 111, 321, -321, 130, 310, 221, 331/
  If ((mstj(51)/=1 .And. mstj(51)/=2) .Or. n-nsav<=1) Return
  Do j = 1, 4
    dps(j) = 0.D0
  End Do
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 120
    Do j = 1, 4
      dps(j) = dps(j) + dble(p(i,j))
    End Do
  120 End Do
  Call ludbrb(0, 0, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))
  pecm = 0.
  Do i = 1, n
    If (k(i,1)>=1 .And. k(i,1)<=10) pecm = pecm + p(i, 4)
  End Do
  nbe(0) = n + mstu(3)
  Do ibe = 1, min(9, mstj(51))
    nbe(ibe) = nbe(ibe-1)
    Do i = nsav + 1, n
      If (k(i,2)/=kfbe(ibe)) Goto 150
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 150
      If (nbe(ibe)>=mstu(4)-mstu(32)-5) Then
        Call luerrm(11, '(LUBOEI:) no more memory left in LUJETS')
        Return
      End If
      nbe(ibe) = nbe(ibe) + 1
      k(nbe(ibe), 1) = i
      Do j = 1, 3
        p(nbe(ibe), j) = 0.
      End Do
    150 End Do
  End Do
  Do ibe = 1, min(9, mstj(51))
    If (ibe/=1 .And. ibe/=4 .And. ibe<=7) Goto 180
    If (ibe==1 .And. max(nbe(1)-nbe(0),nbe(2)-nbe(1),nbe(3)-nbe(2))<=1) Goto 180
    If (ibe==4 .And. max(nbe(4)-nbe(3),nbe(5)-nbe(4),nbe(6)-nbe(5),nbe(7)-nbe(6))<=1) Goto 180
    If (ibe>=8 .And. nbe(ibe)-nbe(ibe-1)<=1) Goto 180
    If (ibe==1) pmhq = 2.*ulmass(211)
    If (ibe==4) pmhq = 2.*ulmass(321)
    If (ibe==8) pmhq = 2.*ulmass(221)
    If (ibe==9) pmhq = 2.*ulmass(331)
    qdel = 0.1*min(pmhq, parj(93))
    If (mstj(51)==1) Then
      nbin = min(100, nint(9.*parj(93)/qdel))
      beex = exp(0.5*qdel/parj(93))
      bert = exp(-qdel/parj(93))
    Else
      nbin = min(100, nint(3.*parj(93)/qdel))
    End If
    Do ibin = 1, nbin
      qbin = qdel*(ibin-0.5)
      bei(ibin) = qdel*(qbin**2+qdel**2/12.)/sqrt(qbin**2+pmhq**2)
      If (mstj(51)==1) Then
        beex = beex*bert
        bei(ibin) = bei(ibin)*beex
      Else
        bei(ibin) = bei(ibin)*exp(-(qbin/parj(93))**2)
      End If
      If (ibin>=2) bei(ibin) = bei(ibin) + bei(ibin-1)
    End Do
    180 Do i1m = nbe(ibe-1) + 1, nbe(ibe) - 1
      i1 = k(i1m, 1)
      Do i2m = i1m + 1, nbe(ibe)
        i2 = k(i2m, 1)
        q2old = max(0., (p(i1,4)+p(i2,4))**2-(p(i1,1)+p(i2,1))**2-(p(i1,2)+p(i2,2))**2-(p(i1,3)+p(i2,3))**2-(p(i1,5)+p(i2,5))**2)
        qold = sqrt(q2old)
        If (qold<0.5*qdel) Then
          qmov = qold/3.
        Else If (qold<(nbin-0.1)*qdel) Then
          rbin = qold/qdel
          ibin = int(rbin)
          rinp = (rbin**3-ibin**3)/(3*ibin*(ibin+1)+1)
          qmov = (bei(ibin)+rinp*(bei(ibin+1)-bei(ibin)))*sqrt(q2old+pmhq**2)/q2old
        Else
          qmov = bei(nbin)*sqrt(q2old+pmhq**2)/q2old
        End If
        q2new = q2old*(qold/(qold+3.*parj(92)*qmov))**(2./3.)
        hc1 = (p(i1,4)+p(i2,4))**2 - (q2old-q2new)
        hc2 = (q2old-q2new)*(p(i1,4)-p(i2,4))**2
        ha = 0.5*(1.-sqrt(hc1*q2new/(hc1*q2old-hc2)))
        Do j = 1, 3
          pd = ha*(p(i2,j)-p(i1,j))
          p(i1m, j) = p(i1m, j) + pd
          p(i2m, j) = p(i2m, j) - pd
        End Do
      End Do
    End Do
  End Do
  Do im = nbe(0) + 1, nbe(min(9,mstj(51)))
    i = k(im, 1)
    Do j = 1, 3
      p(i, j) = p(i, j) + p(im, j)
    End Do
    p(i, 4) = sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  End Do
  pes = 0.
  pqs = 0.
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 240
    pes = pes + p(i, 4)
    pqs = pqs + p(i, 5)**2/p(i, 4)
  240 End Do
  fac = (pecm-pqs)/(pes-pqs)
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 260
    Do j = 1, 3
      p(i, j) = fac*p(i, j)
    End Do
    p(i, 4) = sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  260 End Do
  Call ludbrb(0, 0, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))
  Return
End Subroutine luboei
