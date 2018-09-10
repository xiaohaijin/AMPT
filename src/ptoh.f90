Subroutine ptoh
  Parameter (maxstr=150001)
  Double Precision gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
  Double Precision gxp0, gyp0, gzp0, ft0fom, drlocl
  Double Precision enenew, pxnew, pynew, pznew, beta2, gam
  Double Precision ftavg0, gxavg0, gyavg0, gzavg0, bex, bey, bez
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Double Precision xmdiag, px1, py1, pz1, e1, px2, py2, pz2, e2, px3, py3, pz3, e3, xmpair, etot
  Double Precision p1, p2, p3
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Common /rndf77/nseed
  Common /anim/nevent, isoft, isflag, izpc
  Common /prtn23/gxp0(3), gyp0(3), gzp0(3), ft0fom
  Common /nzpc/nattzp
  Common /lor/enenew, pxnew, pynew, pznew
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /lastt/itimeh, bimp
  Common /hjglbr/nelt, ninthj, nelp, ninp
  Common /arevt/iaevt, iarun, miss
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /input1/masspr, massta, iseed, iavoid, dt
  Dimension xmdiag(maxstr), indx(maxstr), ndiag(maxstr)
  Save
  Call coales
  mstj24 = mstj(24)
  mstj(24) = 0
  nuudd = 0
  npich = 0
  nrhoch = 0
  ppi0 = 1.
  prho0 = 0.
  Do isg = 1, nsg
     If (njsgs(isg)/=0) Then
        natt = natt + 1
        k1 = k2sgs(isg, 1)
        k1abs = iabs(k1)
        px1 = pxsgs(isg, 1)
        py1 = pysgs(isg, 1)
        pz1 = pzsgs(isg, 1)
        k2 = k2sgs(isg, 2)
        k2abs = iabs(k2)
        px2 = pxsgs(isg, 2)
        py2 = pysgs(isg, 2)
        pz2 = pzsgs(isg, 2)
        e1 = pesgs(isg, 1)
        e2 = pesgs(isg, 2)
        xmpair = dsqrt((e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2)
        ibs = 2
        imspin = 0
        If (k1==-k2 .And. iabs(k1)<=2 .And. njsgs(isg)==2) Then
           nuudd = nuudd + 1
           xmdiag(nuudd) = xmpair
           ndiag(nuudd) = natt
        End If
        k3 = 0
        If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           k3 = k2sgs(isg, 3)
           k3abs = iabs(k3)
           px3 = pxsgs(isg, 3)
           py3 = pysgs(isg, 3)
           pz3 = pzsgs(isg, 3)
           e3 = pesgs(isg, 3)
           xmpair = dsqrt((e1+e2+e3)**2-(px1+px2+px3)**2-(py1+py2+py3)**2-(pz1+pz2+pz3)**2)
        End If
        If (isoft==3 .And. (k1abs>1000 .Or. k2abs>1000)) Then
           If (k1abs>1000) Then
              kdq = k1abs
              kk = k2abs
           Else
              kdq = k2abs
              kk = k1abs
           End If
           ki = mod(kdq/1000, 10)
           kj = mod(kdq/100, 10)
           If (mod(kdq,10)==1) Then
              idqspn = 0
           Else
              idqspn = 1
           End If
           If (kk>ki) Then
              ktemp = kk
              kk = kj
              kj = ki
              ki = ktemp
           Else If (kk>kj) Then
              ktemp = kk
              kk = kj
              kj = ktemp
           End If
           If (ki/=kj .And. ki/=kk .And. kj/=kk) Then
              If (idqspn==0) Then
                 kf = 1000*ki + 100*kk + 10*kj + ibs
              Else
                 kf = 1000*ki + 100*kj + 10*kk + ibs
              End If
           Else If (ki==kj .And. ki==kk) Then
              kf = 1000*ki + 100*kj + 10*kk + 4
           Else
              kf = 1000*ki + 100*kj + 10*kk + ibs
           End If
           If (kf==2112 .Or. kf==2212) Then
              If (abs(sngl(xmpair)-ulmass(kf))>abs(sngl(xmpair)-ulmass(kf+2))) kf = kf + 2
           End If
           If (k1<0) kf = -kf
        Else If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           If (k1abs>k2abs) Then
              ki = k1abs
              kk = k2abs
           Else
              ki = k2abs
              kk = k1abs
           End If
           If (k3abs>ki) Then
              kj = ki
              ki = k3abs
           Else If (k3abs<kk) Then
              kj = kk
              kk = k3abs
           Else
              kj = k3abs
           End If
           If (ki==kj .And. ki==kk) Then
              ibs = 4
              kf = 1000*ki + 100*kj + 10*kk + ibs
           Else If (ki/=kj .And. ki/=kk .And. kj/=kk) Then
              ibs = 2
              kf1 = 1000*ki + 100*kj + 10*kk + ibs
              kf2 = 1000*ki + 100*kk + 10*kj + ibs
              kf = kf1
              If (abs(sngl(xmpair)-ulmass(kf1))>abs(sngl(xmpair)-ulmass(kf2))) kf = kf2
           Else
              ibs = 2
              kf = 1000*ki + 100*kj + 10*kk + ibs
              If (kf==2112 .Or. kf==2212) Then
                 If (abs(sngl(xmpair)-ulmass(kf))>abs(sngl(xmpair)-ulmass(kf+2))) kf = kf + 2
              End If
           End If
           If (k1<0) kf = -kf
        Else
           If (k1abs==k2abs) Then
              If (k1abs<=2) Then
                 kf = 0
              Else If (k1abs<=3) Then
                 kf = 333
              Else
                 kf = 100*k1abs + 10*k1abs + 2*imspin + 1
              End If
           Else
              If (k1abs>k2abs) Then
                 kmax = k1abs
                 kmin = k2abs
              Else If (k1abs<k2abs) Then
                 kmax = k2abs
                 kmin = k1abs
              End If
              kf = (100*kmax+10*kmin+2*imspin+1)*isign(1, k1+k2)*(-1)**kmax
              If (mod(iabs(kf),10)==1) Then
                 If (abs(sngl(xmpair)-ulmass(iabs(kf)))>abs(sngl(xmpair)-ulmass(iabs(kf)+2))) kf = (iabs(kf)+2)*isign(1, kf)
              End If
           End If
        End If
        itypar(natt) = kf
        katt(natt, 1) = kf
        If (iabs(kf)==211) Then
           npich = npich + 1
        Else If (iabs(kf)==213) Then
           nrhoch = nrhoch + 1
        End If
     End If
  End Do
  If (nuudd/=0) Then
     ppi0 = float(npich/2)/float(nuudd)
     prho0 = float(nrhoch/2)/float(nuudd)
  End If
  npi0 = 0
  Do isg = 1, nsg
     If (k2sgs(isg,1)==-k2sgs(isg,2) .And. iabs(k2sgs(isg,1))<=2 .And. njsgs(isg)==2) Then
        If (ranart(nseed)<=ppi0) npi0 = npi0 + 1
     End If
  End Do
  If (nuudd>1) Then
     Call index1(maxstr, nuudd, xmdiag, indx)
  Else
     indx(1) = 1
  End If
  Do ix = 1, nuudd
     iuudd = indx(ix)
     inatt = ndiag(iuudd)
     If (ix<=npi0) Then
        kf = 111
     Else If (ranart(nseed)<=(prho0/(1-ppi0+0.00001))) Then
        kf = 113
     Else
        If (ranart(nseed)<=0.5) Then
           kf = 221
        Else
           kf = 223
        End If
     End If
     itypar(inatt) = kf
     katt(inatt, 1) = kf
  End Do
  inatt = 0
  If (ioscar==3) Then
     Write (85, 395) iaevt, 3*nsmbbbar + 2*nsmmeson, nsmbbbar, nsmmeson, bimp, nelp, ninp, nelt, ninthj, miss
  End If
  Do isg = 1, nsg
     If (njsgs(isg)/=0) Then
        inatt = inatt + 1
        k1 = k2sgs(isg, 1)
        k1abs = iabs(k1)
        px1 = pxsgs(isg, 1)
        py1 = pysgs(isg, 1)
        pz1 = pzsgs(isg, 1)
        k2 = k2sgs(isg, 2)
        k2abs = iabs(k2)
        px2 = pxsgs(isg, 2)
        py2 = pysgs(isg, 2)
        pz2 = pzsgs(isg, 2)
        e1 = pesgs(isg, 1)
        e2 = pesgs(isg, 2)
        If (njsgs(isg)==2) Then
           pxar(inatt) = sngl(px1+px2)
           pyar(inatt) = sngl(py1+py2)
           pzar(inatt) = sngl(pz1+pz2)
           patt(inatt, 1) = pxar(inatt)
           patt(inatt, 2) = pyar(inatt)
           patt(inatt, 3) = pzar(inatt)
           etot = e1 + e2
           p1 = px1 + px2
           p2 = py1 + py2
           p3 = pz1 + pz2
        Else If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           px3 = pxsgs(isg, 3)
           py3 = pysgs(isg, 3)
           pz3 = pzsgs(isg, 3)
           e3 = pesgs(isg, 3)
           pxar(inatt) = sngl(px1+px2+px3)
           pyar(inatt) = sngl(py1+py2+py3)
           pzar(inatt) = sngl(pz1+pz2+pz3)
           patt(inatt, 1) = pxar(inatt)
           patt(inatt, 2) = pyar(inatt)
           patt(inatt, 3) = pzar(inatt)
           etot = e1 + e2 + e3
           p1 = px1 + px2 + px3
           p2 = py1 + py2 + py3
           p3 = pz1 + pz2 + pz3
        End If
        xmar(inatt) = ulmass(itypar(inatt))
        kf = katt(inatt, 1)
        If (kf==113 .Or. abs(kf)==213 .Or. kf==221 .Or. kf==223 .Or. abs(kf)==313 .Or. abs(kf)==323 .Or. kf==333 .Or. abs(kf)==1114 .Or. abs(kf)==2114 .Or. abs(kf)==2214 .Or. abs(kf)==2224) Then
           xmar(inatt) = resmass(kf)
        End If
        pear(inatt) = sqrt(pxar(inatt)**2+pyar(inatt)**2+pzar(inatt)**2+xmar(inatt)**2)
        patt(inatt, 4) = pear(inatt)
        eatt = eatt + pear(inatt)
        ipartn = njsgs(isg)
        Do i = 1, ipartn
           ftp(i) = ftsgs(isg, i)
           gxp(i) = gxsgs(isg, i)
           gyp(i) = gysgs(isg, i)
           gzp(i) = gzsgs(isg, i)
           pxp(i) = pxsgs(isg, i)
           pyp(i) = pysgs(isg, i)
           pzp(i) = pzsgs(isg, i)
           pmp(i) = pmsgs(isg, i)
           pep(i) = pesgs(isg, i)
        End Do
        Call locldr(ipartn, drlocl)
        tau0 = arpar1(1)
        ftavg0 = ft0fom + dble(tau0)
        gxavg0 = 0D0
        gyavg0 = 0D0
        gzavg0 = 0D0
        Do i = 1, ipartn
           gxavg0 = gxavg0 + gxp0(i)/ipartn
           gyavg0 = gyavg0 + gyp0(i)/ipartn
           gzavg0 = gzavg0 + gzp0(i)/ipartn
        End Do
        bex = p1/etot
        bey = p2/etot
        bez = p3/etot
        beta2 = bex**2 + bey**2 + bez**2
        gam = 1.D0/dsqrt(1.D0-beta2)
        If (beta2>=0.9999999999999D0) Then
           Write (6, *) '2', bex, bey, bez, beta2, gam
        End If
        Call lorenz(ftavg0, gxavg0, gyavg0, gzavg0, -bex, -bey, -bez)
        gxar(inatt) = sngl(pxnew)
        gyar(inatt) = sngl(pynew)
        gzar(inatt) = sngl(pznew)
        ftar(inatt) = sngl(enenew)
        If (ioscar==3) Then
           Write (85, 313) k2sgs(isg, 1), px1, py1, pz1, pmsgs(isg, 1), inatt, katt(inatt, 1), xmar(inatt)
           Write (85, 312) k2sgs(isg, 2), px2, py2, pz2, pmsgs(isg, 2), inatt, katt(inatt, 1)
           If (njsgs(isg)==3) Write (85, 312) k2sgs(isg, 3), px3, py3, pz3, pmsgs(isg, 3), inatt, katt(inatt, 1)
        End If
     End If
  End Do
  nattzp = natt
  mstj(24) = mstj24
  Return
395 Format (4I8, F10.4, 5I5)
312 Format (I6, 4(1X,F10.3), 1X, I6, 1X, I6)
313 Format (I6, 4(1X,F10.3), 1X, I6, 1X, I6, 1X, F10.3)
End Subroutine ptoh
