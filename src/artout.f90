Subroutine artout(nevnt)
  Parameter (maxstr=150001, maxr=1)
  Parameter (ymt1=-1.0, ymt2=1.0)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  Common /arana1/dy1ntb(50), dy1ntp(50), dy1hm(50), dy1kp(50), dy1km(50), dy1k0s(50), dy1la(50), dy1lb(50), dy1phi(50), dm1pip(50), dm1pim(50), dmt1pr(50), dmt1pb(50), dmt1kp(50), dm1km(50), dm1k0s(50), dmt1la(50), dmt1lb(50), dy1msn(50), dy1pip(50), dy1pim(50), dy1pi0(50), dy1pr(50), dy1pb(50), dy1neg(50), dy1ch(50), de1neg(50), de1ch(50)
  Common /arana2/dy2ntb(50), dy2ntp(50), dy2hm(50), dy2kp(50), dy2km(50), dy2k0s(50), dy2la(50), dy2lb(50), dy2phi(50), dm2pip(50), dm2pim(50), dmt2pr(50), dmt2pb(50), dmt2kp(50), dm2km(50), dm2k0s(50), dmt2la(50), dmt2lb(50), dy2msn(50), dy2pip(50), dy2pim(50), dy2pi0(50), dy2pr(50), dy2pb(50), dy2neg(50), dy2ch(50), de2neg(50), de2ch(50)
  Save
  Open (30, File='../data/dndy_netb.dat', Status='UNKNOWN')
  Open (31, File='../data/dndy_netp.dat', Status='UNKNOWN')
  Open (32, File='../data/dndy_nb.dat', Status='UNKNOWN')
  Open (33, File='../data/dndy_neg.dat', Status='UNKNOWN')
  Open (34, File='../data/dndy_ch.dat', Status='UNKNOWN')
  Open (35, File='../data/dnde_neg.dat', Status='UNKNOWN')
  Open (36, File='../data/dnde_ch.dat', Status='UNKNOWN')
  Open (37, File='../data/dndy_kp.dat', Status='UNKNOWN')
  Open (38, File='../data/dndy_km.dat', Status='UNKNOWN')
  Open (39, File='../data/dndy_k0l.dat', Status='UNKNOWN')
  Open (40, File='../data/dndy_la.dat', Status='UNKNOWN')
  Open (41, File='../data/dndy_lb.dat', Status='UNKNOWN')
  Open (42, File='../data/dndy_phi.dat', Status='UNKNOWN')
  Open (43, File='../data/dndy_meson.dat', Status='UNKNOWN')
  Open (44, File='../data/dndy_pip.dat', Status='UNKNOWN')
  Open (45, File='../data/dndy_pim.dat', Status='UNKNOWN')
  Open (46, File='../data/dndy_pi0.dat', Status='UNKNOWN')
  Open (47, File='../data/dndy_pr.dat', Status='UNKNOWN')
  Open (48, File='../data/dndy_pb.dat', Status='UNKNOWN')
  Open (50, File='../data/dndmtdy_pip.dat', Status='UNKNOWN')
  Open (51, File='../data/dndmtdy_0_1_pim.dat', Status='UNKNOWN')
  Open (52, File='../data/dndmtdy_pr.dat', Status='UNKNOWN')
  Open (53, File='../data/dndmtdy_pb.dat', Status='UNKNOWN')
  Open (54, File='../data/dndmtdy_kp.dat', Status='UNKNOWN')
  Open (55, File='../data/dndmtdy_0_5_km.dat', Status='UNKNOWN')
  Open (56, File='../data/dndmtdy_k0s.dat', Status='UNKNOWN')
  Open (57, File='../data/dndmtdy_la.dat', Status='UNKNOWN')
  Open (58, File='../data/dndmtdy_lb.dat', Status='UNKNOWN')
  scale1 = 1./real(nevnt*num)/by
  scale2 = 1./real(nevnt*num)/bmt/(ymt2-ymt1)
  Do i = 1, 50
    ymid = -10. + by*(i-0.5)
    Write (30, 333) ymid, scale1*dy1ntb(i)
    Write (31, 333) ymid, scale1*dy1ntp(i)
    Write (32, 333) ymid, scale1*dy1hm(i)
    Write (37, 333) ymid, scale1*dy1kp(i)
    Write (38, 333) ymid, scale1*dy1km(i)
    Write (39, 333) ymid, scale1*dy1k0s(i)
    Write (40, 333) ymid, scale1*dy1la(i)
    Write (41, 333) ymid, scale1*dy1lb(i)
    Write (42, 333) ymid, scale1*dy1phi(i)
    Write (33, 333) ymid, scale1*dy1neg(i)
    Write (34, 333) ymid, scale1*dy1ch(i)
    Write (35, 333) ymid, scale1*de1neg(i)
    Write (36, 333) ymid, scale1*de1ch(i)
    Write (43, 333) ymid, scale1*dy1msn(i)
    Write (44, 333) ymid, scale1*dy1pip(i)
    Write (45, 333) ymid, scale1*dy1pim(i)
    Write (46, 333) ymid, scale1*dy1pi0(i)
    Write (47, 333) ymid, scale1*dy1pr(i)
    Write (48, 333) ymid, scale1*dy1pb(i)
    If (dm1pip(i)/=0.0) Then
      Write (50, 333) bmt*(i-0.5), scale2*dm1pip(i)
    End If
    If (dm1pim(i)/=0.0) Then
      Write (51, 333) bmt*(i-0.5), scale2*0.1*dm1pim(i)
    End If
    If (dmt1pr(i)/=0.0) Then
      Write (52, 333) bmt*(i-0.5), scale2*dmt1pr(i)
    End If
    If (dmt1pb(i)/=0.0) Then
      Write (53, 333) bmt*(i-0.5), scale2*dmt1pb(i)
    End If
    If (dmt1kp(i)/=0.0) Then
      Write (54, 333) bmt*(i-0.5), scale2*dmt1kp(i)
    End If
    If (dm1km(i)/=0.0) Then
      Write (55, 333) bmt*(i-0.5), scale2*0.5*dm1km(i)
    End If
    If (dm1k0s(i)/=0.0) Then
      Write (56, 333) bmt*(i-0.5), scale2*dm1k0s(i)
    End If
    If (dmt1la(i)/=0.0) Then
      Write (57, 333) bmt*(i-0.5), scale2*dmt1la(i)
    End If
    If (dmt1lb(i)/=0.0) Then
      Write (58, 333) bmt*(i-0.5), scale2*dmt1lb(i)
    End If
  End Do
  Do i = 30, 48
    Write (i, *) 'after hadron evolution'
  End Do
  Do i = 50, 58
    Write (i, *) 'after hadron evolution'
  End Do
  Do i = 1, 50
    ymid = -10. + by*(i-0.5)
    Write (30, 333) ymid, scale1*dy2ntb(i)
    Write (31, 333) ymid, scale1*dy2ntp(i)
    Write (32, 333) ymid, scale1*dy2hm(i)
    Write (37, 333) ymid, scale1*dy2kp(i)
    Write (38, 333) ymid, scale1*dy2km(i)
    Write (39, 333) ymid, scale1*dy2k0s(i)
    Write (40, 333) ymid, scale1*dy2la(i)
    Write (41, 333) ymid, scale1*dy2lb(i)
    Write (42, 333) ymid, scale1*dy2phi(i)
    Write (33, 333) ymid, scale1*dy2neg(i)
    Write (34, 333) ymid, scale1*dy2ch(i)
    Write (35, 333) ymid, scale1*de2neg(i)
    Write (36, 333) ymid, scale1*de2ch(i)
    Write (43, 333) ymid, scale1*dy2msn(i)
    Write (44, 333) ymid, scale1*dy2pip(i)
    Write (45, 333) ymid, scale1*dy2pim(i)
    Write (46, 333) ymid, scale1*dy2pi0(i)
    Write (47, 333) ymid, scale1*dy2pr(i)
    Write (48, 333) ymid, scale1*dy2pb(i)
    If (dm2pip(i)/=0.0) Then
      Write (50, 333) bmt*(i-0.5), scale2*dm2pip(i)
    End If
    If (dm2pim(i)/=0.0) Then
      Write (51, 333) bmt*(i-0.5), scale2*0.1*dm2pim(i)
    End If
    If (dmt2pr(i)/=0.0) Then
      Write (52, 333) bmt*(i-0.5), scale2*dmt2pr(i)
    End If
    If (dmt2pb(i)/=0.0) Then
      Write (53, 333) bmt*(i-0.5), scale2*dmt2pb(i)
    End If
    If (dmt2kp(i)/=0.0) Then
      Write (54, 333) bmt*(i-0.5), scale2*dmt2kp(i)
    End If
    If (dm2km(i)/=0.0) Then
      Write (55, 333) bmt*(i-0.5), scale2*0.5*dm2km(i)
    End If
    If (dm2k0s(i)/=0.0) Then
      Write (56, 333) bmt*(i-0.5), scale2*dm2k0s(i)
    End If
    If (dmt2la(i)/=0.0) Then
      Write (57, 333) bmt*(i-0.5), scale2*dmt2la(i)
    End If
    If (dmt2lb(i)/=0.0) Then
      Write (58, 333) bmt*(i-0.5), scale2*dmt2lb(i)
    End If
  End Do
  Return
  333 Format (2(F12.5,1X))
End Subroutine artout
