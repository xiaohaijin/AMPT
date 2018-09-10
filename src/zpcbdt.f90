Block Data zpcbdt
   Implicit Double Precision (A-H, O-Z)
   Parameter (maxptn=400001)
   Parameter (maxstr=150001)
   Common /para1/mul
   Common /para2/xmp, xmu, alpha, rscut2, cutof2
   Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
   Common /para4/iftflg, ireflg, igeflg, ibstfg
   Common /para5/iconfg, iordsc
   Common /para6/centy
   Common /para7/ioscar, nsmbbbar, nsmmeson
   Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
   Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
   Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
   Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
   Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
   Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
   Common /aurec1/jxa, jya, jza
   Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
   Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
   Common /ilist2/icell, icel(10, 10, 10)
   Common /ilist3/size1, size2, size3, v1, v2, v3, size
   Common /ilist4/ifmpt, ichkpt, indx(maxptn)
   Common /ilist6/t, iopern, icolln
   Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
   Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
   Common /rndm1/number
   Common /rndm2/iff
   Common /rndm3/iseedp
   Common /ana1/ts(12)
   Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
   Common /ana3/em(4, 4, 12)
   Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
   Save
   Data centy/0D0/
   Data number/0/
   Data ts/0.11D0, 0.12D0, 0.15D0, 0.2D0, 0.3D0, 0.4D0, 0.6D0, 0.8D0, 1D0, 2D0, 4D0, 6D0/
End Block Data
