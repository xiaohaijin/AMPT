Subroutine iniflw(nevnt, idd)
  Parameter (maxstr=150001, maxr=1)
  Double Precision v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Common /run/num
  Common /arerc1/multi1(maxr)
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /iflow/v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Save
  If (idd==0) Then
     v2i = 0D0
     eti = 0D0
     xmulti = 0D0
     v2mi = 0D0
     s2mi = 0D0
     xmmult = 0D0
     v2bi = 0D0
     s2bi = 0D0
     xbmult = 0D0
  Else If (idd==1) Then
     Do j = 1, num
        Do i = 1, multi1(j)
           ityp = ityp1(i, j)
           If (ityp>-100 .And. ityp<100) Goto 100
           xmulti = xmulti + 1.D0
           px = px1(i, j)
           py = py1(i, j)
           xm = xm1(i, j)
           pt2 = px**2 + py**2
           xh = gx1(i, j)
           yh = gy1(i, j)
           xt2 = xh**2 + yh**2
           If (pt2>0) v2i = v2i + dble((px**2-py**2)/pt2)
           eti = eti + dble(sqrt(px**2+py**2+xm**2))
           If (ityp<-1000 .Or. ityp>1000) Then
              xbmult = xbmult + 1.D0
              If (pt2>0) v2bi = v2bi + dble((px**2-py**2)/pt2)
              If (xt2>0) s2bi = s2bi + dble((xh**2-yh**2)/xt2)
           Else
              xmmult = xmmult + 1.D0
              If (pt2>0) v2mi = v2mi + dble((px**2-py**2)/pt2)
              If (xt2>0) s2mi = s2mi + dble((xh**2-yh**2)/xt2)
           End If
100        Continue
        End Do
     End Do
  Else If (idd==2) Then
     If (xmulti/=0) v2i = v2i/xmulti
     eti = eti/dble(nevnt)
     xmulti = xmulti/dble(nevnt)
     If (xmmult/=0) Then
        v2mi = v2mi/xmmult
        s2mi = s2mi/xmmult
     End If
     xmmult = xmmult/dble(nevnt)
     If (xbmult/=0) Then
        v2bi = v2bi/xbmult
        s2bi = s2bi/xbmult
     End If
     xbmult = xbmult/dble(nevnt)
  End If
  Return
End Subroutine iniflw
