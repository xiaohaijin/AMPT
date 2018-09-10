Subroutine minijet_out(bb, phirp)
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arevt/iaevt, iarun, miss
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
  ntrig = 0
  Do i = 1, ihnt2(1)
     Do j = 1, npj(i)
        pt = sqrt(pjpx(i,j)**2+pjpy(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  Do i = 1, ihnt2(3)
     Do j = 1, ntj(i)
        pt = sqrt(pjtx(i,j)**2+pjty(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  Do i = 1, nsg
     Do j = 1, njsg(i)
        pt = sqrt(pxsg(i,j)**2+pysg(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  If (ntrig==0) Return
  If (ioscar==3) Write (96, *) iaevt, miss, ihnt2(1), ihnt2(3)
  Do i = 1, ihnt2(1)
     Do j = 1, npj(i)
        ityp = kfpj(i, j)
        gx = yp(1, i) + 0.5*bb*cos(phirp)
        gy = yp(2, i) + 0.5*bb*sin(phirp)
        gz = 0.
        ft = 0.
        px = pjpx(i, j)
        py = pjpy(i, j)
        pz = pjpz(i, j)
        xmass = pjpm(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 1
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 1
           End If
        End If
     End Do
  End Do
  Do i = 1, ihnt2(3)
     Do j = 1, ntj(i)
        ityp = kftj(i, j)
        gx = yt(1, i) - 0.5*bb*cos(phirp)
        gy = yt(2, i) - 0.5*bb*sin(phirp)
        gz = 0.
        ft = 0.
        px = pjtx(i, j)
        py = pjty(i, j)
        pz = pjtz(i, j)
        xmass = pjtm(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 2
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 2
           End If
        End If
     End Do
  End Do
  Do i = 1, nsg
     Do j = 1, njsg(i)
        ityp = k2sg(i, j)
        gx = 0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
        gy = 0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
        gz = 0.
        ft = 0.
        px = pxsg(i, j)
        py = pysg(i, j)
        pz = pzsg(i, j)
        xmass = pmsg(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 3
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 3
           End If
        End If
     End Do
  End Do
  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 2(1X,F8.2), 2(2X,F2.0), 2X, I2)
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 2(1X,E8.2), 2(2X,F2.0), 2X, I2)
End Subroutine minijet_out
