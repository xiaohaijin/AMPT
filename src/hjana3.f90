Subroutine hjana3
  Parameter (maxstr=150001, maxr=1)
  Parameter (ymin=-1.0, ymax=1.0)
  Parameter (dmt=0.05, dy=0.2)
  Double Precision v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Dimension dndyh3(50), dmyh3(50), deyh3(50)
  Common /run/num
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /arout/iout
  Common /iflow/v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Save
  Data iw/0/
  iw = iw + 1
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      If (ityp>-100 .And. ityp<100) Goto 200
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      xmt = sqrt(px**2+py**2+xm**2)
      dxmt = xmt - xm
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA3 mt=0'
        y = 1000000.0*sign(1., pz)
      End If
      iy = 1 + int(y/dy)
      If (iy<1 .Or. iy>50) Goto 100
      dndyh3(iy) = dndyh3(iy) + 1.0
      deyh3(iy) = deyh3(iy) + xmt
      100 Continue
      If (y<ymin .Or. y>=ymax) Goto 200
      imt = 1 + int(dxmt/dmt)
      If (imt>50) Goto 200
      dmyh3(imt) = dmyh3(imt) + 1.0/xmt
      200 Continue
    End Do
  End Do
  Return
End Subroutine hjana3
