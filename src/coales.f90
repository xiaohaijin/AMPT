Subroutine coales
  Parameter (maxstr=150001)
  Implicit Double Precision (D)
  Double Precision gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
  Dimension iover(maxstr), dp1(2:3), dr1(2:3)
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Double Precision dpcoal, drcoal, ecritl
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Common /coal/dpcoal, drcoal, ecritl
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Save
  Do isg = 1, nsg
     iover(isg) = 0
  End Do
  Do isg = 1, nsg
     If (njsgs(isg)/=2 .Or. iover(isg)==1) Goto 150
     If (k2sgs(isg,1)<0) Then
        Write (6, *) 'Antiquark appears in quark loop; stop'
        Stop
     End If
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr0 = drlocl
     dp0 = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     Do jsg = 1, nsg
        If (jsg==isg .Or. iover(jsg)==1) Goto 120
        If (njsgs(jsg)==2) Then
           ipmin = 2
           ipmax = 2
        Else If (njsgs(jsg)==3 .And. k2sgs(jsg,1)<0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 120
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(1)*pesgs(jsg,ip)-pxp(1)*pxsgs(jsg,ip)-pyp(1)*pysgs(jsg,ip)-pzp(1)*pzsgs(jsg,ip)-pmp(1)*pmsgs(jsg,ip)))
           If (dplocl>dpcoal) Goto 120
           ftp(2) = ftsgs(jsg, ip)
           gxp(2) = gxsgs(jsg, ip)
           gyp(2) = gysgs(jsg, ip)
           gzp(2) = gzsgs(jsg, ip)
           pxp(2) = pxsgs(jsg, ip)
           pyp(2) = pysgs(jsg, ip)
           pzp(2) = pzsgs(jsg, ip)
           pmp(2) = pmsgs(jsg, ip)
           pep(2) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           If (drlocl>drcoal) Goto 120
           If ((dp0>dpcoal .Or. dr0>drcoal) .Or. (drlocl<dr0)) Then
              dp0 = dplocl
              dr0 = drlocl
              Call exchge(isg, 2, jsg, ip)
           End If
        End Do
120  End Do
     If (dp0<=dpcoal .And. dr0<=drcoal) iover(isg) = 1
150 End Do
  Do isg = 1, nsg
     If (njsgs(isg)/=2 .Or. iover(isg)==1) Goto 250
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr0 = drlocl
     dp0 = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     Do jsg = 1, nsg
        If (jsg==isg .Or. iover(jsg)==1) Goto 220
        If (njsgs(jsg)==2) Then
           ipmin = 1
           ipmax = 1
        Else If (njsgs(jsg)==3 .And. k2sgs(jsg,1)>0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 220
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(2)*pesgs(jsg,ip)-pxp(2)*pxsgs(jsg,ip)-pyp(2)*pysgs(jsg,ip)-pzp(2)*pzsgs(jsg,ip)-pmp(2)*pmsgs(jsg,ip)))
           If (dplocl>dpcoal) Goto 220
           ftp(1) = ftsgs(jsg, ip)
           gxp(1) = gxsgs(jsg, ip)
           gyp(1) = gysgs(jsg, ip)
           gzp(1) = gzsgs(jsg, ip)
           pxp(1) = pxsgs(jsg, ip)
           pyp(1) = pysgs(jsg, ip)
           pzp(1) = pzsgs(jsg, ip)
           pmp(1) = pmsgs(jsg, ip)
           pep(1) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           If (drlocl>drcoal) Goto 220
           If ((dp0>dpcoal .Or. dr0>drcoal) .Or. (drlocl<dr0)) Then
              dp0 = dplocl
              dr0 = drlocl
              Call exchge(isg, 1, jsg, ip)
           End If
        End Do
220  End Do
     If (dp0<=dpcoal .And. dr0<=drcoal) iover(isg) = 1
250 End Do
  Do isg = 1, nsg
     If (njsgs(isg)/=3 .Or. iover(isg)==1) Goto 350
     ibaryn = k2sgs(isg, 1)
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr1(2) = drlocl
     dp1(2) = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     ftp(2) = ftsgs(isg, 3)
     gxp(2) = gxsgs(isg, 3)
     gyp(2) = gysgs(isg, 3)
     gzp(2) = gzsgs(isg, 3)
     pxp(2) = pxsgs(isg, 3)
     pyp(2) = pysgs(isg, 3)
     pzp(2) = pzsgs(isg, 3)
     pmp(2) = pmsgs(isg, 3)
     pep(2) = pesgs(isg, 3)
     Call locldr(2, drlocl)
     dr1(3) = drlocl
     dp1(3) = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     Do jsg = 1, nsg
        If (jsg==isg .Or. iover(jsg)==1) Goto 320
        If (njsgs(jsg)==2) Then
           If (ibaryn>0) Then
              ipmin = 1
           Else
              ipmin = 2
           End If
           ipmax = ipmin
        Else If (njsgs(jsg)==3 .And. (ibaryn*k2sgs(jsg,1))>0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 320
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(1)*pesgs(jsg,ip)-pxp(1)*pxsgs(jsg,ip)-pyp(1)*pysgs(jsg,ip)-pzp(1)*pzsgs(jsg,ip)-pmp(1)*pmsgs(jsg,ip)))
           If (dplocl>dpcoal) Goto 320
           ftp(2) = ftsgs(jsg, ip)
           gxp(2) = gxsgs(jsg, ip)
           gyp(2) = gysgs(jsg, ip)
           gzp(2) = gzsgs(jsg, ip)
           pxp(2) = pxsgs(jsg, ip)
           pyp(2) = pysgs(jsg, ip)
           pzp(2) = pzsgs(jsg, ip)
           pmp(2) = pmsgs(jsg, ip)
           pep(2) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           If (drlocl>drcoal) Goto 320
           ipi = 0
           If (dp1(2)>dpcoal .Or. dr1(2)>drcoal) Then
              ipi = 2
              If ((dp1(3)>dpcoal .Or. dr1(3)>drcoal) .And. dr1(3)>dr1(2)) ipi = 3
           Else If (dp1(3)>dpcoal .Or. dr1(3)>drcoal) Then
              ipi = 3
           Else If (dr1(2)<dr1(3)) Then
              If (drlocl<dr1(3)) ipi = 3
           Else If (dr1(3)<=dr1(2)) Then
              If (drlocl<dr1(2)) ipi = 2
           End If
           If (ipi/=0) Then
              dp1(ipi) = dplocl
              dr1(ipi) = drlocl
              Call exchge(isg, ipi, jsg, ip)
           End If
        End Do
320  End Do
     If (dp1(2)<=dpcoal .And. dr1(2)<=drcoal .And. dp1(3)<=dpcoal .And. dr1(3)<=drcoal) iover(isg) = 1
350 End Do
  Return
End Subroutine coales
