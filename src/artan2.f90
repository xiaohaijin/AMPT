Subroutine artan2
  Parameter (maxstr=150001, maxr=1)
  Parameter (ymt1=-1.0, ymt2=1.0)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /arana2/dy2ntb(50), dy2ntp(50), dy2hm(50), dy2kp(50), dy2km(50), dy2k0s(50), dy2la(50), dy2lb(50), dy2phi(50), dm2pip(50), dm2pim(50), dmt2pr(50), dmt2pb(50), dmt2kp(50), dm2km(50), dm2k0s(50), dmt2la(50), dmt2lb(50), dy2msn(50), dy2pip(50), dy2pim(50), dy2pi0(50), dy2pr(50), dy2pb(50), dy2neg(50), dy2ch(50), de2neg(50), de2ch(50)
  Save
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      xmt = sqrt(px**2+py**2+xm**2)
      If (xm<0.01) Goto 200
      dxmt = xmt - xm
      ptot = sqrt(px**2+py**2+pz**2)
      If ((px**2+py**2)>0.) Then
        eta = asinh(pz/sqrt(px**2+py**2))
      Else
        eta = 1000000.0*sign(1., pz)
        If (abs(pz)<=1E-3) eta = 0.
      End If
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN ARTAN2 mt=0'
        y = 1000000.0*sign(1., pz)
      End If
      If (abs(y)>=10.0) Goto 100
      If (abs(eta)>=10.0) Goto 100
      iy = 1 + int((y+10.)/by)
      ieta = 1 + int((eta+10.)/by)
      If (ityp<-1000) Then
        dy2ntb(iy) = dy2ntb(iy) - 1.0
      End If
      If (ityp>1000) Then
        dy2ntb(iy) = dy2ntb(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy2ntp(iy) = dy2ntp(iy) - 1.0
      End If
      If (ityp==2212) Then
        dy2ntp(iy) = dy2ntp(iy) + 1.0
      End If
      If (ityp==-2112) Then
        dy2hm(iy) = dy2hm(iy) + 1.0
      End If
      If (luchge(ityp)/=0) Then
        dy2ch(iy) = dy2ch(iy) + 1.0
        de2ch(ieta) = de2ch(ieta) + 1.0
        If (luchge(ityp)<0) Then
          dy2neg(iy) = dy2neg(iy) + 1.0
          de2neg(ieta) = de2neg(ieta) + 1.0
        End If
      End If
      If ((ityp>=100 .And. ityp<1000) .Or. (ityp>-1000 .And. ityp<=-100)) Then
        dy2msn(iy) = dy2msn(iy) + 1.0
      End If
      If (ityp==211) Then
        dy2pip(iy) = dy2pip(iy) + 1.0
      End If
      If (ityp==-211) Then
        dy2pim(iy) = dy2pim(iy) + 1.0
      End If
      If (ityp==111) Then
        dy2pi0(iy) = dy2pi0(iy) + 1.0
      End If
      If (ityp==2212) Then
        dy2pr(iy) = dy2pr(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy2pb(iy) = dy2pb(iy) + 1.0
      End If
      If (ityp==321) Then
        dy2kp(iy) = dy2kp(iy) + 1.0
      End If
      If (ityp==-321) Then
        dy2km(iy) = dy2km(iy) + 1.0
      End If
      If (ityp==130) Then
        dy2k0s(iy) = dy2k0s(iy) + 1.0
      End If
      If (ityp==3122) Then
        dy2la(iy) = dy2la(iy) + 1.0
      End If
      If (ityp==-3122) Then
        dy2lb(iy) = dy2lb(iy) + 1.0
      End If
      If (ityp==333) Then
        dy2phi(iy) = dy2phi(iy) + 1.0
      End If
      100 If (y<ymt1 .Or. y>ymt2) Goto 200
      If (dxmt>=50.0*bmt .Or. dxmt==0) Goto 200
      imt = 1 + int(dxmt/bmt)
      If (ityp==211) Then
        dm2pip(imt) = dm2pip(imt) + 1.0/xmt
      End If
      If (ityp==-211) Then
        dm2pim(imt) = dm2pim(imt) + 1.0/xmt
      End If
      If (ityp==2212) Then
        dmt2pr(imt) = dmt2pr(imt) + 1.0/xmt
      End If
      If (ityp==-2212) Then
        dmt2pb(imt) = dmt2pb(imt) + 1.0/xmt
      End If
      If (ityp==321) Then
        dmt2kp(imt) = dmt2kp(imt) + 1.0/xmt
      End If
      If (ityp==-321) Then
        dm2km(imt) = dm2km(imt) + 1.0/xmt
      End If
      If (ityp==130) Then
        dm2k0s(imt) = dm2k0s(imt) + 1.0/xmt
      End If
      If (ityp==3122) Then
        dmt2la(imt) = dmt2la(imt) + 1.0/xmt
      End If
      If (ityp==-3122) Then
        dmt2lb(imt) = dmt2lb(imt) + 1.0/xmt
      End If
      200 Continue
    End Do
  End Do
  Return
End Subroutine artan2
