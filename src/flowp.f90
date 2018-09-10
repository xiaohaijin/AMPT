Subroutine flowp(idd)
  Implicit Double Precision (A-H, O-Z)
  Real dt
  Parameter (maxptn=400001)
  Parameter (bmt=0.05D0)
  Dimension nlfile(3), nsfile(3), nmfile(3)
  Dimension v2pp(3), xnpp(3), v2psum(3), v2p2sm(3), nfile(3)
  Dimension tsp(31), v2pevt(3), v2pavg(3), varv2p(3)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /para1/mul
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /pflow/v2p(30, 3), xnpart(30, 3), etp(30, 3), s2p(30, 3), v2p2(30, 3), nevt(30)
  Common /pflowf/v2pf(30, 3), xnpf(30, 3), etpf(30, 3), xncoll(30), s2pf(30, 3), v2pf2(30, 3)
  Common /pfrz/v2pfrz(30, 3), xnpfrz(30, 3), etpfrz(30, 3), s2pfrz(30, 3), v2p2fz(30, 3), tscatt(31), nevtfz(30), iscatt(30)
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  Common /arevt/iaevt, iarun, miss
  Common /anim/nevent, isoft, isflag, izpc
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /precpb/vxp(maxptn), vyp(maxptn), vzp(maxptn)
  Save
  Dimension etpl(30, 3), etps(30, 3), etplf(30, 3), etpsf(30, 3), etlfrz(30, 3), etsfrz(30, 3), xnpl(30, 3), xnps(30, 3), xnplf(30, 3), xnpsf(30, 3), xnlfrz(30, 3), xnsfrz(30, 3), v2pl(30, 3), v2ps(30, 3), v2plf(30, 3), v2psf(30, 3), s2pl(30, 3), s2ps(30, 3), s2plf(30, 3), s2psf(30, 3), dmyil(50, 3), dmyfl(50, 3), dmyis(50, 3), dmyfs(50, 3)
  Data tsp/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30/
  If (idd==0) Then
     nfile(1) = 60
     nfile(2) = 64
     nfile(3) = 20
     Open (nfile(1), File='ana1/v2p.dat', Status='UNKNOWN')
     Open (nfile(1)+1, File='ana1/v2p-formed.dat', Status='UNKNOWN')
     Open (nfile(1)+2, File='ana1/v2p-active.dat', Status='UNKNOWN')
     Open (nfile(1)+3, File='ana1/v2ph.dat', Status='UNKNOWN')
     Open (nfile(2), File='ana1/v2p-y2.dat', Status='UNKNOWN')
     Open (nfile(2)+1, File='ana1/v2p-formed2.dat', Status='UNKNOWN')
     Open (nfile(2)+2, File='ana1/v2p-active2.dat', Status='UNKNOWN')
     Open (nfile(2)+3, File='ana1/v2ph-y2.dat', Status='UNKNOWN')
     Open (nfile(3), File='ana1/v2p-y1.dat', Status='UNKNOWN')
     Open (nfile(3)+1, File='ana1/v2p-formed1.dat', Status='UNKNOWN')
     Open (nfile(3)+2, File='ana1/v2p-active1.dat', Status='UNKNOWN')
     Open (nfile(3)+3, File='ana1/v2ph-y1.dat', Status='UNKNOWN')
     Open (49, File='ana1/v2p-ebe.dat', Status='UNKNOWN')
     Write (49, *) '    ievt,  v2p,  v2p_y2,   v2p_y1'
     Open (59, File='ana1/v2h.dat', Status='UNKNOWN')
     Open (68, File='ana1/v2h-y2.dat', Status='UNKNOWN')
     Open (69, File='ana1/v2h-y1.dat', Status='UNKNOWN')
     Open (88, File='ana1/v2h-ebe.dat', Status='UNKNOWN')
     Write (88, *) '    ievt,  v2h,  v2h_y2,   v2h_y1'
     nlfile(1) = 70
     nlfile(2) = 72
     nlfile(3) = 74
     Open (nlfile(1), File='ana1/mtl.dat', Status='UNKNOWN')
     Open (nlfile(1)+1, File='ana1/mtl-formed.dat', Status='UNKNOWN')
     Open (nlfile(2), File='ana1/mtl-y2.dat', Status='UNKNOWN')
     Open (nlfile(2)+1, File='ana1/mtl-formed2.dat', Status='UNKNOWN')
     Open (nlfile(3), File='ana1/mtl-y1.dat', Status='UNKNOWN')
     Open (nlfile(3)+1, File='ana1/mtl-formed1.dat', Status='UNKNOWN')
     nsfile(1) = 76
     nsfile(2) = 78
     nsfile(3) = 80
     Open (nsfile(1), File='ana1/mts.dat', Status='UNKNOWN')
     Open (nsfile(1)+1, File='ana1/mts-formed.dat', Status='UNKNOWN')
     Open (nsfile(2), File='ana1/mts-y2.dat', Status='UNKNOWN')
     Open (nsfile(2)+1, File='ana1/mts-formed2.dat', Status='UNKNOWN')
     Open (nsfile(3), File='ana1/mts-y1.dat', Status='UNKNOWN')
     Open (nsfile(3)+1, File='ana1/mts-formed1.dat', Status='UNKNOWN')
     nmfile(1) = 82
     nmfile(2) = 83
     nmfile(3) = 84
     Open (nmfile(1), File='ana1/Nmt.dat', Status='UNKNOWN')
     Open (nmfile(2), File='ana1/Nmt-y2.dat', Status='UNKNOWN')
     Open (nmfile(3), File='ana1/Nmt-y1.dat', Status='UNKNOWN')
     ifanim = 0
     If (ifanim==1) Then
        Open (10, File='ana1/h-animate.dat', Status='UNKNOWN')
        Write (10, *) ntmax, dt
        Open (11, File='ana1/p-animate.dat', Status='UNKNOWN')
        Open (15, File='ana1/p-finalft.dat', Status='UNKNOWN')
     End If
     If (nevent<1) Open (93, File='ana1/parton-t.dat', Status='UNKNOWN')
     itimep = 0
     itanim = 0
     iaevtp = 0
     Do ii = 1, 50
        Do iy = 1, 3
           dmyil(ii, iy) = 0D0
           dmyfl(ii, iy) = 0D0
           dmyis(ii, iy) = 0D0
           dmyfs(ii, iy) = 0D0
        End Do
     End Do
     Do ii = 1, 31
        tscatt(ii) = 0D0
     End Do
     Do ii = 1, 30
        nevt(ii) = 0
        xncoll(ii) = 0D0
        nevtfz(ii) = 0
        iscatt(ii) = 0
        Do iy = 1, 3
           v2p(ii, iy) = 0D0
           v2p2(ii, iy) = 0D0
           s2p(ii, iy) = 0D0
           etp(ii, iy) = 0D0
           xnpart(ii, iy) = 0D0
           v2pf(ii, iy) = 0D0
           v2pf2(ii, iy) = 0D0
           s2pf(ii, iy) = 0D0
           etpf(ii, iy) = 0D0
           xnpf(ii, iy) = 0D0
           v2pfrz(ii, iy) = 0D0
           v2p2fz(ii, iy) = 0D0
           s2pfrz(ii, iy) = 0D0
           etpfrz(ii, iy) = 0D0
           xnpfrz(ii, iy) = 0D0
           etpl(ii, iy) = 0D0
           etps(ii, iy) = 0D0
           etplf(ii, iy) = 0D0
           etpsf(ii, iy) = 0D0
           etlfrz(ii, iy) = 0D0
           etsfrz(ii, iy) = 0D0
           xnpl(ii, iy) = 0D0
           xnps(ii, iy) = 0D0
           xnplf(ii, iy) = 0D0
           xnpsf(ii, iy) = 0D0
           xnlfrz(ii, iy) = 0D0
           xnsfrz(ii, iy) = 0D0
           v2pl(ii, iy) = 0D0
           v2ps(ii, iy) = 0D0
           v2plf(ii, iy) = 0D0
           v2psf(ii, iy) = 0D0
           s2pl(ii, iy) = 0D0
           s2ps(ii, iy) = 0D0
           s2plf(ii, iy) = 0D0
           s2psf(ii, iy) = 0D0
        End Do
     End Do
     Do iy = 1, 3
        v2pevt(iy) = 0D0
        v2pavg(iy) = 0D0
        varv2p(iy) = 0D0
        v2pp(iy) = 0.D0
        xnpp(iy) = 0D0
        v2psum(iy) = 0.D0
        v2p2sm(iy) = 0.D0
     End Do
  Else If (idd==1) Then
     If (iaevt/=iaevtp .And. ianp==31) itanim = 0
     t2time = ft5(iscat)
     Do ianp = 1, 30
        If (t2time<tsp(ianp+1) .And. t2time>=tsp(ianp)) Then
           xncoll(ianp) = xncoll(ianp) + 1D0
           If (ianp<=itimep .And. iaevt==iaevtp) Goto 101
           nevt(ianp) = nevt(ianp) + 1
           tscatt(ianp+1) = t2time
           iscatt(ianp) = 1
           nevtfz(ianp) = nevtfz(ianp) + 1
           Do i = 1, mul
              delta = 1D-8
              If ((e5(i)-dabs(pz5(i))+delta)<=0) Then
                 ypartn = 1000000.D0*sign(1.D0, pz5(i))
                 Write (6, *) 'ypartn error', e5(i) - dabs(pz5(i))
              Else
                 ypartn = 0.5D0*dlog((e5(i)+pz5(i)+delta)/(e5(i)-pz5(i)+delta))
              End If
              pt2 = px5(i)**2 + py5(i)**2
              iloop = 1
              If (dabs(ypartn)<=1D0) Then
                 iloop = 2
                 If (dabs(ypartn)<=0.5D0) Then
                    iloop = 3
                 End If
              End If
              Do iy = 1, iloop
                 If (pt2>0D0) Then
                    v2prtn = (px5(i)**2-py5(i)**2)/pt2
                    If (dabs(v2prtn)>1D0) Write (nfile(iy), *) 'v2prtn>1', v2prtn
                    v2p(ianp, iy) = v2p(ianp, iy) + v2prtn
                    v2p2(ianp, iy) = v2p2(ianp, iy) + v2prtn**2
                 End If
                 xperp2 = gx5(i)**2 + gy5(i)**2
                 If (xperp2>0D0) s2p(ianp, iy) = s2p(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                 xnpart(ianp, iy) = xnpart(ianp, iy) + 1D0
                 etp(ianp, iy) = etp(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                 If (ft5(i)<=t2time) Then
                    v2pf(ianp, iy) = v2pf(ianp, iy) + v2prtn
                    v2pf2(ianp, iy) = v2pf2(ianp, iy) + v2prtn**2
                    If (xperp2>0D0) s2pf(ianp, iy) = s2pf(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                    xnpf(ianp, iy) = xnpf(ianp, iy) + 1D0
                    etpf(ianp, iy) = etpf(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                 End If
              End Do
           End Do
           itimep = ianp
           iaevtp = iaevt
           If (ianp==30) Then
              Do iy = 1, 3
                 npartn = idint(xnpart(ianp,iy)-xnpp(iy))
                 If (npartn/=0) Then
                    v2pevt(iy) = (v2p(ianp,iy)-v2pp(iy))/npartn
                    v2psum(iy) = v2psum(iy) + v2pevt(iy)
                    v2p2sm(iy) = v2p2sm(iy) + v2pevt(iy)**2
                    v2pp(iy) = v2p(ianp, iy)
                    xnpp(iy) = xnpart(ianp, iy)
                 End If
              End Do
              Write (49, 160) iaevt, v2pevt
           End If
           Goto 101
        End If
     End Do
101  If (nevent<1) Then
        Do nt = 1, ntmax
           time1 = dble(nt*dt)
           time2 = dble((nt+1)*dt)
           If (t2time<time2 .And. t2time>=time1) Then
              If (nt<=itanim) Return
              If (ifanim==1) Write (11, *) t2time
              iform = 0
              ne1all = 0
              ne1form = 0
              Do i = 1, mul
                 gz_now = gz5(i) + (t2time-ft5(i))*pz5(i)/e5(i)
                 If (dabs(gz_now)<t2time) Then
                    etap = 0.5D0*dlog((t2time+gz_now)/(t2time-gz_now))
                 Else
                    etap = 1000000.D0*sign(1.D0, gz_now)
                 End If
                 ne1all = ne1all + 1
                 If (ft5(i)<=t2time) ne1form = ne1form + 1
                 If (ft5(i)<=t2time) iform = iform + 1
              End Do
              If (ifanim==1) Write (11, *) iform
              Write (93, 184) 'evt#,t,np,npformed=', iaevt, t2time, ne1all, ne1form
              Do i = 1, mul
                 If (ft5(i)<=t2time) Then
                    gz_now = gz5(i) + (t2time-ft5(i))*pz5(i)/e5(i)
                 Else
                    gz_now = gz5(i) + (t2time-ft5(i))*vzp(i)
                 End If
                 If (dabs(gz_now)<t2time) Then
                    etap = 0.5D0*dlog((t2time+gz_now)/(t2time-gz_now))
                 Else
                    etap = 1000000.D0*sign(1.D0, gz_now)
                 End If
                 If (ft5(i)<=t2time) Then
                    gx_now = gx5(i) + (t2time-ft5(i))*px5(i)/e5(i)
                    gy_now = gy5(i) + (t2time-ft5(i))*py5(i)/e5(i)
                 Else
                    gx_now = gx5(i) + (t2time-ft5(i))*vxp(i)
                    gy_now = gy5(i) + (t2time-ft5(i))*vyp(i)
                 End If
                 Write (93, 185) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx_now, gy_now, ft5(i), etap
                 If (ifanim==1 .And. ft5(i)<=t2time) Then
                    Write (11, 180) ityp5(i), gx5(i), gy5(i), gz5(i), ft5(i), px5(i), py5(i), pz5(i), e5(i)
                 End If
              End Do
              itanim = nt
           End If
        End Do
     End If
  Else If (idd==3) Then
     Do ianp = 1, 30
        If (iscatt(ianp)==0) tscatt(ianp+1) = tscatt(ianp)
     End Do
     Do i = 1, mul
        delta = 1D-8
        If ((e5(i)-dabs(pz5(i))+delta)<=0) Then
           Write (6, *) 'ypartn error', e5(i) - dabs(pz5(i))
           ypartn = 1000000.D0*sign(1.D0, pz5(i))
        Else
           ypartn = 0.5D0*dlog((e5(i)+pz5(i)+delta)/(e5(i)-pz5(i)+delta))
        End If
        pt2 = px5(i)**2 + py5(i)**2
        iloop = 1
        If (dabs(ypartn)<=1D0) Then
           iloop = 2
           If (dabs(ypartn)<=0.5D0) Then
              iloop = 3
           End If
        End If
        Do ianp = 1, 30
           If (iscatt(ianp)/=0) Then
              If (ft5(i)<tscatt(ianp+1) .And. ft5(i)>=tscatt(ianp)) Then
                 Do iy = 1, iloop
                    If (pt2>0D0) Then
                       v2prtn = (px5(i)**2-py5(i)**2)/pt2
                       v2pfrz(ianp, iy) = v2pfrz(ianp, iy) + v2prtn
                       v2p2fz(ianp, iy) = v2p2fz(ianp, iy) + v2prtn**2
                    End If
                    xperp2 = gx5(i)**2 + gy5(i)**2
                    If (xperp2>0D0) s2pfrz(ianp, iy) = s2pfrz(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                    etpfrz(ianp, iy) = etpfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                    xnpfrz(ianp, iy) = xnpfrz(ianp, iy) + 1D0
                    If (ityp5(i)==1 .Or. ityp5(i)==2) Then
                       etlfrz(ianp, iy) = etlfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                       xnlfrz(ianp, iy) = xnlfrz(ianp, iy) + 1D0
                    Else If (ityp5(i)==3) Then
                       etsfrz(ianp, iy) = etsfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                       xnsfrz(ianp, iy) = xnsfrz(ianp, iy) + 1D0
                    End If
                 End Do
                 Goto 350
              End If
           End If
        End Do
350  End Do
  Else If (idd==2) Then
     Do i = 1, 3
        Write (nfile(i), *) '   tsp,   v2p,     v2p2, ' // '   s2p,  etp,   xmult,    nevt,  xnctot'
        Write ((nfile(i)+1), *) '  tsp,   v2pf,   v2pf2, ' // '   s2pf, etpf,  xnform,  xcrate'
        Write ((nfile(i)+2), *) '  tsp,   v2pa,   v2pa2, ' // '   s2pa, etpa,  xmult_ap,  xnform,   nevt'
        Write ((nfile(i)+3), *) '  tsph,  v2ph,   v2ph2, ' // '   s2ph, etph,  xmult_(ap/2+h),xmult_ap/2,nevt'
        Write (nlfile(i), *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nsfile(i), *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nlfile(i)+1, *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nsfile(i)+1, *) '   tsp,    v2,     s2,    etp,    xmul'
     End Do
     Do ii = 1, 30
        If (nevt(ii)==0) Goto 150
        Do iy = 1, 3
           If (xnpart(ii,iy)>1D0) Then
              v2p(ii, iy) = v2p(ii, iy)/xnpart(ii, iy)
              v2p2(ii, iy) = dsqrt((v2p2(ii,iy)/xnpart(ii,iy)-v2p(ii,iy)**2)/(xnpart(ii,iy)-1))
              s2p(ii, iy) = s2p(ii, iy)/xnpart(ii, iy)
              xmult = dble(xnpart(ii,iy)/dble(nevt(ii)))
              etp(ii, iy) = etp(ii, iy)/dble(nevt(ii))
              etpl(ii, iy) = etpl(ii, iy)/dble(nevt(ii))
              etps(ii, iy) = etps(ii, iy)/dble(nevt(ii))
              xnctot = 0D0
              Do inum = 1, ii
                 xnctot = xnctot + xncoll(inum)
              End Do
              If (nevt(1)/=0) xnctot = xnctot/nevt(1)
              Write (nfile(iy), 200) tsp(ii), v2p(ii, iy), v2p2(ii, iy), s2p(ii, iy), etp(ii, iy), xmult, nevt(ii), xnctot
           End If
           If (nevt(ii)/=0) xcrate = xncoll(ii)/(tsp(ii+1)-tsp(ii))/nevt(ii)
           If (xnpf(ii,iy)>1D0) Then
              v2pf(ii, iy) = v2pf(ii, iy)/xnpf(ii, iy)
              v2pf2(ii, iy) = dsqrt((v2pf2(ii,iy)/xnpf(ii,iy)-v2pf(ii,iy)**2)/(xnpf(ii,iy)-1))
              s2pf(ii, iy) = s2pf(ii, iy)/xnpf(ii, iy)
              xnform = dble(xnpf(ii,iy)/dble(nevt(ii)))
              etpf(ii, iy) = etpf(ii, iy)/dble(nevt(ii))
              etplf(ii, iy) = etplf(ii, iy)/dble(nevt(ii))
              etpsf(ii, iy) = etpsf(ii, iy)/dble(nevt(ii))
              Write (nfile(iy)+1, 210) tsp(ii), v2pf(ii, iy), v2pf2(ii, iy), s2pf(ii, iy), etpf(ii, iy), xnform, xcrate
           End If
           If (xnpl(ii,iy)>1D0) Then
              v2pl(ii, iy) = v2pl(ii, iy)/xnpl(ii, iy)
              s2pl(ii, iy) = s2pl(ii, iy)/xnpl(ii, iy)
              xmult = dble(xnpl(ii,iy)/dble(nevt(ii)))
              etpl(ii, iy) = etpl(ii, iy)/dble(nevt(ii))
              Write (nlfile(iy), 201) tsp(ii), v2pl(ii, iy), s2pl(ii, iy), etpl(ii, iy), xmult
           End If
           If (xnps(ii,iy)>1D0) Then
              v2ps(ii, iy) = v2ps(ii, iy)/xnps(ii, iy)
              s2ps(ii, iy) = s2ps(ii, iy)/xnps(ii, iy)
              xmult = dble(xnps(ii,iy)/dble(nevt(ii)))
              etps(ii, iy) = etps(ii, iy)/dble(nevt(ii))
              Write (nsfile(iy), 201) tsp(ii), v2ps(ii, iy), s2ps(ii, iy), etps(ii, iy), xmult
           End If
           If (xnplf(ii,iy)>1D0) Then
              v2plf(ii, iy) = v2plf(ii, iy)/xnplf(ii, iy)
              s2plf(ii, iy) = s2plf(ii, iy)/xnplf(ii, iy)
              xmult = dble(xnplf(ii,iy)/dble(nevt(ii)))
              etplf(ii, iy) = etplf(ii, iy)/dble(nevt(ii))
              Write (nlfile(iy)+1, 201) tsp(ii), v2plf(ii, iy), s2plf(ii, iy), etplf(ii, iy), xmult
           End If
           If (xnpsf(ii,iy)>1D0) Then
              v2psf(ii, iy) = v2psf(ii, iy)/xnpsf(ii, iy)
              s2psf(ii, iy) = s2psf(ii, iy)/xnpsf(ii, iy)
              xmult = dble(xnpsf(ii,iy)/dble(nevt(ii)))
              etpsf(ii, iy) = etpsf(ii, iy)/dble(nevt(ii))
              Write (nsfile(iy)+1, 201) tsp(ii), v2psf(ii, iy), s2psf(ii, iy), etpsf(ii, iy), xmult
           End If
        End Do
150  End Do
     scalei = 0D0
     scalef = 0D0
     If (nevt(1)/=0) scalei = 1D0/dble(nevt(1))/bmt
     If (nevt(30)/=0) scalef = 1D0/dble(nevt(30))/bmt
     Do iy = 2, 3
        yra = 1D0
        If (iy==2) yra = 2D0
        Do i = 1, 50
           Write (nmfile(iy), 251) bmt*dble(i-0.5), scalei*dmyil(i, iy)/yra, scalef*dmyfl(i, iy)/yra, scalei*dmyis(i, iy)/yra, scalef*dmyfs(i, iy)/yra
        End Do
     End Do
     If (nevt(30)>=1) Then
        Do iy = 1, 3
           v2pavg(iy) = v2psum(iy)/nevt(30)
           v2var0 = v2p2sm(iy)/nevt(30) - v2pavg(iy)**2
           If (v2var0>0D0) varv2p(iy) = dsqrt(v2var0)
        End Do
        Write (49, 240) 'EBE v2p,v2p(y2),v2p(y1): avg=', v2pavg
        Write (49, 240) 'EBE v2p,v2p(y2),v2p(y1): var=', varv2p
     End If
     If (ifanim==1) Then
        Do i = 1, mul
           If (ft5(i)<=t2time) Then
              Write (15, 140) ityp5(i), gx5(i), gy5(i), gz5(i), ft5(i)
           End If
        End Do
        Write (10, *) - 10.
        Write (10, *) 0
        Write (11, *) - 10.
        Write (11, *) 0
        Close (10)
        Close (11)
        Close (15)
     End If
     Do ianp = 1, 30
        Do iy = 1, 3
           v2pact = 0D0
           v2p2ac = 0D0
           s2pact = 0D0
           etpact = 0D0
           xnacti = 0D0
           If (xnpf(ianp,iy)>1D0) Then
              v2pact = v2pf(ianp, iy)*xnpf(ianp, iy)
              v2p2ac = (v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)+v2pf(ianp,iy)**2)*xnpf(ianp, iy)
              s2pact = s2pf(ianp, iy)*xnpf(ianp, iy)
              etpact = etpf(ianp, iy)*dble(nevt(ianp))
              xnpact = xnpf(ianp, iy)
              Do kanp = 1, ianp
                 v2pact = v2pact - v2pfrz(kanp, iy)
                 v2p2ac = v2p2ac - v2p2fz(kanp, iy)
                 s2pact = s2pact - s2pfrz(kanp, iy)
                 etpact = etpact - etpfrz(kanp, iy)
                 xnpact = xnpact - xnpfrz(kanp, iy)
              End Do
              v2ph = v2pact
              v2ph2 = v2p2ac
              s2ph = s2pact
              etph = etpact
              xnp2 = xnpact/2D0
              If (xnpact>1D0 .And. nevt(ianp)/=0) Then
                 v2pact = v2pact/xnpact
                 v2p2ac = dsqrt((v2p2ac/xnpact-v2pact**2)/(xnpact-1))
                 s2pact = s2pact/xnpact
                 xnacti = dble(xnpact/dble(nevt(ianp)))
                 etpact = etpact/dble(nevt(ianp))
                 Write (nfile(iy)+2, 250) tsp(ianp), v2pact, v2p2ac, s2pact, etpact, xnacti, xnpf(ianp, iy)/dble(nevt(ianp)), nevt(ianp)
              End If
           End If
           shadr = dble(nevt(ianp))/dble(nevent)
           ianh = ianp
           v2ph = v2ph + v2h(ianh, iy)*xnhadr(ianh, iy)*shadr
           v2ph2 = v2ph2 + (v2h2(ianh,iy)**2*(xnhadr(ianh,iy)-1)+v2h(ianh,iy)**2)*xnhadr(ianh, iy)*shadr
           s2ph = s2ph + s2h(ianh, iy)*xnhadr(ianh, iy)*shadr
           etph = etph + eth(ianh, iy)*dble(nevent)*shadr
           xnph = xnpact + xnhadr(ianh, iy)*shadr
           xnp2h = xnp2 + xnhadr(ianh, iy)*shadr
           If (xnph>1D0 .And. nevt(ianp)/=0) Then
              v2ph = v2ph/xnph
              v2ph2 = dsqrt((v2ph2/xnph-v2ph**2)/(xnph-1))
              s2ph = s2ph/xnph
              etph = etph/dble(nevt(ianp))
              xnp2 = xnp2/dble(nevt(ianp))
              xnp2h = xnp2h/dble(nevent)
              If (tsp(ianp)<=dble(ntmax*dt)) Write (nfile(iy)+3, 250) tsp(ianp), v2ph, v2ph2, s2ph, etph, xnp2h, xnp2, nevt(ianp)
           End If
        End Do
     End Do
     Do ianp = 1, 30
        Do iy = 1, 3
           v2pact = 0D0
           v2p2ac = 0D0
           s2pact = 0D0
           etpact = 0D0
           xnacti = 0D0
           v2pact = v2pf(ianp, iy)*xnpf(ianp, iy)
           v2p2ac = (v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)+v2pf(ianp,iy)**2)*xnpf(ianp, iy)
           s2pact = s2pf(ianp, iy)*xnpf(ianp, iy)
           etpact = etpf(ianp, iy)*dble(nevt(ianp))
           xnpact = xnpf(ianp, iy)
        End Do
     End Do
     Close (620)
     Close (630)
     Do nf = 1, 3
        Do ifile = 0, 3
           Close (nfile(nf)+ifile)
        End Do
     End Do
     Do nf = 1, 3
        Close (740+nf)
     End Do
  End If
  Return
184 Format (A20, I7, F8.4, 2(1X,I6))
185 Format (I3, 3(1X,F8.3), 1X, F8.4, 1X, 2(F8.3,1X), F11.4, 1X, F8.3)
160 Format (I10, 3(2X,F9.5))
180 Format (I6, 8(1X,F7.2))
140 Format (I10, 4(2X,F7.2))
200 Format (2X, F5.2, 3(2X,F7.4), 2(2X,F9.2), I6, 2X, F9.2)
210 Format (2X, F5.2, 3(2X,F7.4), 3(2X,F9.2))
240 Format (A30, 3(2X,F9.5))
250 Format (2X, F5.2, 3(2X,F7.4), 3(2X,F9.2), I6)
201 Format (2X, F5.2, 4(2X,F9.2))
251 Format (5E15.5)
End Subroutine flowp
