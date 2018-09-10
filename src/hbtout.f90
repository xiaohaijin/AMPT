Subroutine hbtout(nnew, nt, ntmax)
  Parameter (maxstr=150001, maxr=1)
  Parameter (oneminus=0.99999, oneplus=1.00001)
  Dimension lastkp(maxstr), newkp(maxstr), xnew(3)
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /lastt/itimeh, bimp
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  Common /arevt/iaevt, iarun, miss
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Common /hjglbr/nelt, ninthj, nelp, ninp
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  External iarflv, invflv
  Common /para8/idpert, npertd, idxsec
  Common /phihj/iphirp, phirp
  Save
  Do i = 1, max0(nlast, nnew)
     lastkp(i) = 0
  End Do
  Do i = 1, nnew
     newkp(i) = 0
  End Do
  Do ip = 1, nnew
     Do iplast = 1, nlast
        If (p(1,ip)==plast(1,iplast) .And. p(2,ip)==plast(2,iplast) .And. p(3,ip)==plast(3,iplast) .And. e(ip)==plast(4,iplast) .And. lb(ip)==lblast(iplast) .And. dpertp(ip)==dplast(iplast) .And. lastkp(iplast)==0) Then
           deltat = nt*dt - xlast(4, iplast)
           ene = sqrt(plast(1,iplast)**2+plast(2,iplast)**2+plast(3,iplast)**2+plast(4,iplast)**2)
           Do ii = 1, 3
              xnew(ii) = xlast(ii, iplast) + plast(ii, iplast)/ene*deltat
           End Do
           dr = sqrt((r(1,ip)-xnew(1))**2+(r(2,ip)-xnew(2))**2+(r(3,ip)-xnew(3))**2)
           If (dr<=0.01) Then
              lastkp(iplast) = 1
              newkp(ip) = 1
              If (nt==ntmax .And. ftsv(ip)>((ntmax-1)*dt)) xlast(4, iplast) = ftsv(ip)
              Goto 100
           End If
        End If
     End Do
100 End Do
  Do ip = 1, nnew
     If (newkp(ip)==0) Then
        Do iplast = 1, nnew
           If (lastkp(iplast)==0) Then
              xlast(1, iplast) = r(1, ip)
              xlast(2, iplast) = r(2, ip)
              xlast(3, iplast) = r(3, ip)
              xlast(4, iplast) = nt*dt
              If (nt==ntmax) Then
                 If (tfdcy(ip)>(ntmax*dt+0.001)) Then
                    xlast(4, iplast) = tfdcy(ip)
                 Else If (ftsv(ip)>((ntmax-1)*dt)) Then
                    xlast(4, iplast) = ftsv(ip)
                 End If
              End If
              plast(1, iplast) = p(1, ip)
              plast(2, iplast) = p(2, ip)
              plast(3, iplast) = p(3, ip)
              plast(4, iplast) = e(ip)
              lblast(iplast) = lb(ip)
              lastkp(iplast) = 1
              dplast(iplast) = dpertp(ip)
              Goto 150
           End If
        End Do
     End If
150 End Do
  If (nnew<nlast) Then
     Do iplast = 1, nlast
        If (lastkp(iplast)==0) Then
           Do ip2 = iplast + 1, nlast
              If (lastkp(ip2)==1) Then
                 xlast(1, iplast) = xlast(1, ip2)
                 xlast(2, iplast) = xlast(2, ip2)
                 xlast(3, iplast) = xlast(3, ip2)
                 xlast(4, iplast) = xlast(4, ip2)
                 plast(1, iplast) = plast(1, ip2)
                 plast(2, iplast) = plast(2, ip2)
                 plast(3, iplast) = plast(3, ip2)
                 plast(4, iplast) = plast(4, ip2)
                 lblast(iplast) = lblast(ip2)
                 lastkp(iplast) = 1
                 dplast(iplast) = dplast(ip2)
                 Goto 170
              End If
           End Do
        End If
170  End Do
  End If
  nlast = nnew
  If (nt==ntmax) Then
     ndpert = 0
     Do ip = 1, nlast
        If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
        Else
           ndpert = ndpert + 1
        End If
     End Do
     Write (16, 191) iaevt, iarun, nlast - ndpert, bimp, npart1, npart2, nelp, ninp, nelt, ninthj, phirp
     If (idpert==1 .Or. idpert==2) Write (90, 190) iaevt, iarun, ndpert, bimp, npart1, npart2, nelp, ninp, nelt, ninthj
     Do ip = 1, nlast
        If (abs(plast(1,ip))<=epsipt .And. abs(plast(2,ip))<=epsipt .And. (plast(3,ip)>amax1(0.,pzproj-epsipz) .Or. plast(3,ip)<(-pztarg+epsipz)) .And. (lblast(ip)==1 .Or. lblast(ip)==2)) Then
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 200) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 250) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        Else If (amax1(abs(xlast(1,ip)),abs(xlast(2,ip)),abs(xlast(3,ip)),abs(xlast(4,ip)))<9999) Then
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 200) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 250) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip), dplast(ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        Else
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 201) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 251) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip), dplast(ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        End If
     End Do
     If (ioscar==1) Call hoscar
  End If
  Return
190 Format (3(I7), F10.4, 5X, 6(I4))
191 Format (3(I7), F10.4, 5X, 6(I4), 5X, F7.4)
200 Format (I6, 2(1X,F8.3), 1X, F11.4, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F11.4, 1X, F6.3, 4(1X,E8.2))
250 Format (I5, 2(1X,F8.3), 1X, F10.3, 2(1X,F7.1), 1X, F8.2, 1X, F7.2, 1X, E10.4)
251 Format (I5, 2(1X,F8.3), 1X, F10.3, 4(1X,E8.2), 1X, E10.4)
End Subroutine hbtout
