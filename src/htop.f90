Subroutine htop
  Parameter (maxstr=150001)
  Parameter (maxptn=400001)
  Parameter (maxidl=4001)
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs, ptwo, xmdq, ptwox, ptwoy, ptwoz
  Dimension it(4)
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /para1/mul
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /decom/ptwo(2, 5)
  Common /rndf77/nseed
  Common /noprec/nnozpc, itypn(maxidl), gxn(maxidl), gyn(maxidl), gzn(maxidl), ftn(maxidl), pxn(maxidl), pyn(maxidl), pzn(maxidl), een(maxidl), xmn(maxidl)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Common /anim/nevent, isoft, isflag, izpc
  Double Precision vxp0, vyp0, vzp0, xstrg0, ystrg0, xstrg, ystrg
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /arevt/iaevt, iarun, miss
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Save
  npar = 0
  nnozpc = 0
  If ((isoft==4 .Or. isoft==5) .And. (ioscar==2 .Or. ioscar==3)) Then
     nsmbbbar = 0
     nsmmeson = 0
     Do i = 1, natt
        id = itypar(i)
        idabs = iabs(id)
        i2 = mod(idabs/10, 10)
        If (abs(pxar(i))<=epsipt .And. abs(pyar(i))<=epsipt .And. (pzar(i)>amax1(0.,pzproj-epsipz) .Or. pzar(i)<(-pztarg+epsipz)) .And. (id==2112 .Or. id==2212)) Then
        Else If (idabs>1000 .And. i2/=0) Then
           nsmbbbar = nsmbbbar + 1
        Else If ((idabs>100 .And. idabs<1000) .Or. idabs>10000) Then
           nsmmeson = nsmmeson + 1
        End If
     End Do
     If (ioscar==2 .Or. ioscar==3) Then
        Write (92, *) iaevt, miss, 3*nsmbbbar + 2*nsmmeson, nsmbbbar, nsmmeson, natt, natt - nsmbbbar - nsmmeson
     End If
  End If
  Do i = 1, natt
     id = itypar(i)
     idabs = iabs(id)
     i4 = mod(idabs/1000, 10)
     i3 = mod(idabs/100, 10)
     i2 = mod(idabs/10, 10)
     i1 = mod(idabs, 10)
     rnum = ranart(nseed)
     ftime = 0.197*pear(i)/(pxar(i)**2+pyar(i)**2+xmar(i)**2)
     inozpc = 0
     it(1) = 0
     it(2) = 0
     it(3) = 0
     it(4) = 0
     If (abs(pxar(i))<=epsipt .And. abs(pyar(i))<=epsipt .And. (pzar(i)>amax1(0.,pzproj-epsipz) .Or. pzar(i)<(-pztarg+epsipz)) .And. (id==2112 .Or. id==2212)) Then
        inozpc = 1
     Else If (idabs>1000 .And. i2/=0) Then
        If (((i4==1 .Or. i4==2) .And. i4==i3) .Or. (i4==3 .And. i3==3)) Then
           If (i1==2) Then
              If (rnum<=(1./2.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 1
              Else If (rnum<=(2./3.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              Else
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              End If
           Else If (i1==4) Then
              If (rnum<=(2./3.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              Else
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              End If
           End If
        Else If (i4==1 .Or. i4==2) Then
           If (i1==2) Then
              If (rnum<=(1./2.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 1
              Else If (rnum<=(2./3.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              Else
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              End If
           Else If (i1==4) Then
              If (rnum<=(2./3.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              Else
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              End If
           End If
        Else If (i4>=3) Then
           it(1) = i4
           If (i3<i2) Then
              it(2) = i2*1000 + i3*100 + 1
           Else
              it(2) = i3*1000 + i2*100 + 3
           End If
        End If
        If (id<0) Then
           it(1) = -it(1)
           it(2) = -it(2)
        End If
        If (isoft==4 .Or. isoft==5) Then
           it(3) = mod(it(2)/1000, 10)
           it(4) = mod(it(2)/100, 10)
        End If
     Else If ((idabs>100 .And. idabs<1000) .Or. idabs>10000) Then
        If (i3==i2) Then
           If (i3==1 .Or. i3==2) Then
              If (rnum<=0.5) Then
                 it(1) = 1
                 it(2) = -1
              Else
                 it(1) = 2
                 it(2) = -2
              End If
           Else
              it(1) = i3
              it(2) = -i3
           End If
        Else
           If ((isign(1,id)*(-1)**i3)==1) Then
              it(1) = i3
              it(2) = -i2
           Else
              it(1) = i2
              it(2) = -i3
           End If
        End If
     Else
        inozpc = 1
     End If
     If (inozpc==1) Then
        njsgs(i) = 0
        nnozpc = nnozpc + 1
        itypn(nnozpc) = itypar(i)
        pxn(nnozpc) = pxar(i)
        pyn(nnozpc) = pyar(i)
        pzn(nnozpc) = pzar(i)
        een(nnozpc) = pear(i)
        xmn(nnozpc) = xmar(i)
        gxn(nnozpc) = gxar(i)
        gyn(nnozpc) = gyar(i)
        gzn(nnozpc) = gzar(i)
        ftn(nnozpc) = ftar(i)
     Else
        njsgs(i) = 2
        ptwo(1, 5) = dble(ulmass(it(1)))
        ptwo(2, 5) = dble(ulmass(it(2)))
        Call decomp(dble(patt(i,1)), dble(patt(i,2)), dble(patt(i,3)), dble(xmar(i)), i, it(1))
        ipamax = 2
        If ((isoft==4 .Or. isoft==5) .And. iabs(it(2))>1000) ipamax = 1
        Do ipar = 1, ipamax
           npar = npar + 1
           ityp0(npar) = it(ipar)
           px0(npar) = ptwo(ipar, 1)
           py0(npar) = ptwo(ipar, 2)
           pz0(npar) = ptwo(ipar, 3)
           e0(npar) = ptwo(ipar, 4)
           xmass0(npar) = ptwo(ipar, 5)
           gx0(npar) = dble(gxar(i))
           gy0(npar) = dble(gyar(i))
           gz0(npar) = dble(gzar(i))
           ft0(npar) = dble(ftime)
           lstrg0(npar) = i
           lpart0(npar) = ipar
           vxp0(npar) = dble(patt(i,1)/patt(i,4))
           vyp0(npar) = dble(patt(i,2)/patt(i,4))
           vzp0(npar) = dble(patt(i,3)/patt(i,4))
           xstrg(npar) = xstrg0(i)
           ystrg(npar) = ystrg0(i)
           istrg(npar) = istrg0(i)
        End Do
        If ((isoft==4 .Or. isoft==5) .And. iabs(it(2))>1000) Then
           njsgs(i) = 3
           xmdq = ptwo(2, 5)
           ptwo(1, 5) = dble(ulmass(it(3)))
           ptwo(2, 5) = dble(ulmass(it(4)))
           ptwox = ptwo(2, 1)
           ptwoy = ptwo(2, 2)
           ptwoz = ptwo(2, 3)
           Call decomp(ptwox, ptwoy, ptwoz, xmdq, i, it(1))
           Do ipar = 1, 2
              npar = npar + 1
              ityp0(npar) = it(ipar+2)
              px0(npar) = ptwo(ipar, 1)
              py0(npar) = ptwo(ipar, 2)
              pz0(npar) = ptwo(ipar, 3)
              e0(npar) = ptwo(ipar, 4)
              xmass0(npar) = ptwo(ipar, 5)
              gx0(npar) = dble(gxar(i))
              gy0(npar) = dble(gyar(i))
              gz0(npar) = dble(gzar(i))
              ft0(npar) = dble(ftime)
              lstrg0(npar) = i
              lpart0(npar) = ipar + 1
              vxp0(npar) = dble(patt(i,1)/patt(i,4))
              vyp0(npar) = dble(patt(i,2)/patt(i,4))
              vzp0(npar) = dble(patt(i,3)/patt(i,4))
              xstrg(npar) = xstrg0(i)
              ystrg(npar) = ystrg0(i)
              istrg(npar) = istrg0(i)
           End Do
        End If
     End If
  End Do
  mul = npar
  If ((isoft==4 .Or. isoft==5) .And. (ioscar==2 .Or. ioscar==3)) Then
     If ((natt-nsmbbbar-nsmmeson)/=nnozpc) Write (92, *) 'Problem with the total # of initial particles (gamma,e,muon,...) not entering ZPC'
     If ((3*nsmbbbar+2*nsmmeson)/=npar) Write (92, *) 'Problem with the total # of initial partons   after string melting'
  End If
  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
End Subroutine htop
