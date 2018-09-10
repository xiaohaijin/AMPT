Subroutine artan1
  Parameter (maxstr=150001, maxr=1)
  Parameter (ymt1=-1.0, ymt2=1.0)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /arana1/dy1ntb(50), dy1ntp(50), dy1hm(50), dy1kp(50), dy1km(50), dy1k0s(50), dy1la(50), dy1lb(50), dy1phi(50), dm1pip(50), dm1pim(50), dmt1pr(50), dmt1pb(50), dmt1kp(50), dm1km(50), dm1k0s(50), dmt1la(50), dmt1lb(50), dy1msn(50), dy1pip(50), dy1pim(50), dy1pi0(50), dy1pr(50), dy1pb(50), dy1neg(50), dy1ch(50), de1neg(50), de1ch(50)
  Save
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      If (xm<0.01) Goto 200
      ptot = sqrt(px**2+py**2+pz**2)
      If ((px**2+py**2)>0.) Then
        eta = asinh(pz/sqrt(px**2+py**2))
      Else
        eta = 1000000.0*sign(1., pz)
        If (abs(pz)<=1E-3) eta = 0.
      End If
      xmt = sqrt(px**2+py**2+xm**2)
      dxmt = xmt - xm
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN ARTAN1 mt=0'
        y = 1000000.0*sign(1., pz)
      End If
      If (abs(y)>=10.0) Goto 100
      If (abs(eta)>=10.0) Goto 100
      iy = 1 + int((y+10.)/by)
      ieta = 1 + int((eta+10.)/by)
      If (ityp<-1000) Then
        dy1ntb(iy) = dy1ntb(iy) - 1.0
      End If
      If (ityp>1000) Then
        dy1ntb(iy) = dy1ntb(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy1ntp(iy) = dy1ntp(iy) - 1.0
      End If
      If (ityp==2212) Then
        dy1ntp(iy) = dy1ntp(iy) + 1.0
      End If
      If (ityp==-2112) Then
        dy1hm(iy) = dy1hm(iy) + 1.0
      End If
      If (luchge(ityp)/=0) Then
        dy1ch(iy) = dy1ch(iy) + 1.0
        de1ch(ieta) = de1ch(ieta) + 1.0
        If (luchge(ityp)<0) Then
          dy1neg(iy) = dy1neg(iy) + 1.0
          de1neg(ieta) = de1neg(ieta) + 1.0
        End If
      End If
      If ((ityp>=100 .And. ityp<1000) .Or. (ityp>-1000 .And. ityp<=-100)) Then
        dy1msn(iy) = dy1msn(iy) + 1.0
      End If
      If (ityp==211) Then
        dy1pip(iy) = dy1pip(iy) + 1.0
      End If
      If (ityp==-211) Then
        dy1pim(iy) = dy1pim(iy) + 1.0
      End If
      If (ityp==111) Then
        dy1pi0(iy) = dy1pi0(iy) + 1.0
      End If
      If (ityp==2212) Then
        dy1pr(iy) = dy1pr(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy1pb(iy) = dy1pb(iy) + 1.0
      End If
      If (ityp==321) Then
        dy1kp(iy) = dy1kp(iy) + 1.0
      End If
      If (ityp==-321) Then
        dy1km(iy) = dy1km(iy) + 1.0
      End If
      If (ityp==130) Then
        dy1k0s(iy) = dy1k0s(iy) + 1.0
      End If
      If (ityp==3122) Then
        dy1la(iy) = dy1la(iy) + 1.0
      End If
      If (ityp==-3122) Then
        dy1lb(iy) = dy1lb(iy) + 1.0
      End If
      If (ityp==333) Then
        dy1phi(iy) = dy1phi(iy) + 1.0
      End If
      100 If (y<ymt1 .Or. y>ymt2) Goto 200
      If (dxmt>=50.0*bmt .Or. dxmt==0) Goto 200
      imt = 1 + int(dxmt/bmt)
      If (ityp==211) Then
        dm1pip(imt) = dm1pip(imt) + 1.0/xmt
      End If
      If (ityp==-211) Then
        dm1pim(imt) = dm1pim(imt) + 1.0/xmt
      End If
      If (ityp==2212) Then
        dmt1pr(imt) = dmt1pr(imt) + 1.0/xmt
      End If
      If (ityp==-2212) Then
        dmt1pb(imt) = dmt1pb(imt) + 1.0/xmt
      End If
      If (ityp==321) Then
        dmt1kp(imt) = dmt1kp(imt) + 1.0/xmt
      End If
      If (ityp==-321) Then
        dm1km(imt) = dm1km(imt) + 1.0/xmt
      End If
      If (ityp==130) Then
        dm1k0s(imt) = dm1k0s(imt) + 1.0/xmt
      End If
      If (ityp==3122) Then
        dmt1la(imt) = dmt1la(imt) + 1.0/xmt
      End If
      If (ityp==-3122) Then
        dmt1lb(imt) = dmt1lb(imt) + 1.0/xmt
      End If
      200 Continue
    End Do
  End Do
  Return
End Subroutine artan1
