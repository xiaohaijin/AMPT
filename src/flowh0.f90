Subroutine flowh0(nevnt, idd)
  Dimension tsh(31)
  Double Precision v2h, xnhadr, eth, v2h2, s2h
  Double Precision v2hp, xnhadp, v2hsum, v2h2sm, v2havg(3), varv2h(3)
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  Common /ebe/v2hp(3), xnhadp(3), v2hsum(3), v2h2sm(3)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /lastt/itimeh, bimp
  Save
  If (idd==0) Then
     itimeh = 0
     Do ii = 1, 31
        tsh(ii) = float(ii-1)
     End Do
     Do ii = 1, 30
        Do iy = 1, 3
           v2h(ii, iy) = 0D0
           xnhadr(ii, iy) = 0D0
           eth(ii, iy) = 0D0
           v2h2(ii, iy) = 0D0
           s2h(ii, iy) = 0D0
        End Do
     End Do
     Do iy = 1, 3
        v2hp(iy) = 0D0
        xnhadp(iy) = 0D0
        v2hsum(iy) = 0D0
        v2h2sm(iy) = 0D0
        If (iy==1) Then
           nunit = 59
        Else If (iy==2) Then
           nunit = 68
        Else
           nunit = 69
        End If
        Write (nunit, *) '   tsh,   v2h,     v2h2,     s2h, ' // ' eth,   xmulth'
     End Do
  Else If (idd==2) Then
     Do ii = 1, 30
        Do iy = 1, 3
           If (xnhadr(ii,iy)==0) Then
              xmulth = 0.
           Else If (xnhadr(ii,iy)>1) Then
              v2h(ii, iy) = v2h(ii, iy)/xnhadr(ii, iy)
              eth(ii, iy) = eth(ii, iy)/dble(nevnt)
              v2h2(ii, iy) = dsqrt((v2h2(ii,iy)/xnhadr(ii,iy)-v2h(ii,iy)**2)/(xnhadr(ii,iy)-1))
              s2h(ii, iy) = s2h(ii, iy)/xnhadr(ii, iy)
              xmulth = sngl(xnhadr(ii,iy)/nevnt)
           End If
           If (iy==1) Then
              nunit = 59
           Else If (iy==2) Then
              nunit = 68
           Else
              nunit = 69
           End If
           If (tsh(ii)<=(ntmax*dt)) Write (nunit, 200) tsh(ii), v2h(ii, iy), v2h2(ii, iy), s2h(ii, iy), eth(ii, iy), xmulth
        End Do
     End Do
     Do iy = 1, 3
        v2havg(iy) = v2hsum(iy)/dble(nevnt)
        varv2h(iy) = dsqrt(v2h2sm(iy)/dble(nevnt)-v2havg(iy)**2)
     End Do
     Write (88, 240) 'EBE v2h,v2h(y2),v2h(y1): avg=', v2havg
     Write (88, 240) 'EBE v2h,v2h(y2),v2h(y1): var=', varv2h
  End If
  Return
200 Format (2X, F5.2, 3(2X,F7.4), 2(2X,F9.2))
240 Format (A30, 3(2X,F9.5))
End Subroutine flowh0
