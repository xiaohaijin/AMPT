Subroutine flowh(ct)
  Parameter (maxstr=150001, maxr=1)
  Dimension tsh(31)
  Double Precision v2h, xnhadr, eth, v2h2, s2h
  Double Precision v2hp, xnhadp, v2hsum, v2h2sm, v2hevt(3)
  Double Precision pt2, v2hadr
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  Common /ebe/v2hp(3), xnhadp(3), v2hsum(3), v2h2sm(3)
  Common /lastt/itimeh, bimp
  Common /run/num
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /rr/massr(0:maxr)
  Common /anim/nevent, isoft, isflag, izpc
  Common /arevt/iaevt, iarun, miss
  Save
  Do ii = 1, 31
     tsh(ii) = float(ii-1)
  End Do
  Do ianh = 1, 30
     If ((ct+0.0001)<tsh(ianh+1) .And. (ct+0.0001)>=tsh(ianh)) Then
        If (ianh==itimeh) Goto 101
        ia = 0
        Do j = 1, num
           mult = massr(j)
           ia = ia + massr(j-1)
           Do ic = 1, mult
              i = ia + ic
              If (iabs(lb(i)-10000)<100) Goto 100
              px = p(1, i)
              py = p(2, i)
              pt2 = dble(px)**2 + dble(py)**2
              ene = sqrt(e(i)**2+sngl(pt2)+p(3,i)**2)
              rap = 0.5*alog((ene+p(3,i))/(ene-p(3,i)))
              iloop = 1
              If (abs(rap)<=1) Then
                 iloop = 2
                 If (abs(rap)<=0.5) Then
                    iloop = 3
                 End If
              End If
              Do iy = 1, iloop
                 If (pt2>0D0) Then
                    v2hadr = (dble(px)**2-dble(py)**2)/pt2
                    v2h(ianh, iy) = v2h(ianh, iy) + v2hadr
                    v2h2(ianh, iy) = v2h2(ianh, iy) + v2hadr**2
                    If (dabs(v2hadr)>1D0) Write (1, *) 'v2hadr>1', v2hadr, px, py
                 End If
                 xperp2 = r(1, i)**2 + r(2, i)**2
                 If (xperp2>0.) s2h(ianh, iy) = s2h(ianh, iy) + dble((r(1,i)**2-r(2,i)**2)/xperp2)
                 eth(ianh, iy) = eth(ianh, iy) + dble(sqrt(e(i)**2+sngl(pt2)))
                 xnhadr(ianh, iy) = xnhadr(ianh, iy) + 1D0
              End Do
100        End Do
        End Do
        itimeh = ianh
        If (ianh==30) Then
           Do iy = 1, 3
              nhadrn = idint(xnhadr(ianh,iy)-xnhadp(iy))
              If (nhadrn/=0) Then
                 v2hevt(iy) = (v2h(ianh,iy)-v2hp(iy))/dble(nhadrn)
                 v2hsum(iy) = v2hsum(iy) + v2hevt(iy)
                 v2h2sm(iy) = v2h2sm(iy) + v2hevt(iy)**2
                 v2hp(iy) = v2h(ianh, iy)
                 xnhadp(iy) = xnhadr(ianh, iy)
              End If
           End Do
           Write (88, 160) iaevt, v2hevt
        End If
        Goto 101
     End If
  End Do
101 ifanim = 0
  If (ifanim==1) Then
     ia = 0
     Do j = 1, num
        mult = massr(j)
        ia = ia + massr(j-1)
        Write (10, *) ct
        Write (10, *) mult
        Do ic = 1, mult
           i = ia + ic
           If (amax1(abs(r(1,i)),abs(r(2,i)),abs(r(3,i)))<9999) Then
              Write (10, 210) lb(i), r(1, i), r(2, i), r(3, i), p(1, i), p(2, i), p(3, i), e(i)
           Else
              Write (10, 220) lb(i), r(1, i), r(2, i), r(3, i), p(1, i), p(2, i), p(3, i), e(i)
           End If
        End Do
     End Do
     Return
  End If
  Return
160 Format (I10, 3(2X,F9.5))
210 Format (I6, 7(1X,F9.3))
220 Format (I6, 3(1X,E9.3), 4(1X,F9.3))
End Subroutine flowh
