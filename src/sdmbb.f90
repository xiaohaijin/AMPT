Subroutine sdmbb(srt, sdm, ianti)
  Parameter (amn=0.939457, amp=0.93828, am0=1.232, am1440=1.44, am1535=1.535, srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /dpifsl/lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1, lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2, lbpp1, lbpp2
  Common /dpifsm/xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2, xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1, xmpp2
  Common /dpisig/sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp, sdmpp
  Common /para8/idpert, npertd, idxsec
  Common /rndf77/nseed
  Save
  sdm = 0.
  sdmel = 0.
  sdmnn = 0.
  sdmnd = 0.
  sdmns = 0.
  sdmnp = 0.
  sdmdd = 0.
  sdmds = 0.
  sdmdp = 0.
  sdmss = 0.
  sdmsp = 0.
  sdmpp = 0.
  If (srt<=(em1+em2)) Return
  s = srt**2
  pinitial = sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt
  fs = fnndpi(s)
  If (idxsec==1 .Or. idxsec==2) Then
    If ((lb1>=3 .And. lb1<=5) .Or. (lb2>=3 .And. lb2<=5)) Then
      xnnfactor = 8./9.
    Else If ((lb1>=25 .And. lb1<=27) .Or. (lb2>=25 .And. lb2<=27)) Then
      xnnfactor = 8./27.
    Else If (lb1==28 .Or. lb2==28) Then
      xnnfactor = 8./9.
    Else If (lb1==0 .Or. lb2==0) Then
      xnnfactor = 8./3.
    End If
  Else
  End If
  If (idxsec==1 .Or. idxsec==3) Then
    sdmel = fdpiel(s)
  Else If (idxsec==2 .Or. idxsec==4) Then
    threshold = em1 + em2
    snew = (srt-threshold+srt0)**2
    sdmel = fdpiel(snew)
  End If
  If (((lb1==5 .Or. lb2==5 .Or. lb1==27 .Or. lb2==27) .And. ianti==0) .Or. ((lb1==3 .Or. lb2==3 .Or. lb1==25 .Or. lb2==25) .And. ianti==1)) Then
    lbnn1 = 1
    lbnn2 = 1
    xmnn1 = amp
    xmnn2 = amp
  Else If (lb1==3 .Or. lb2==3 .Or. lb1==26 .Or. lb2==26 .Or. lb1==28 .Or. lb2==28 .Or. lb1==0 .Or. lb2==0) Then
    lbnn1 = 2
    lbnn2 = 1
    xmnn1 = amn
    xmnn2 = amp
  Else
    lbnn1 = 2
    lbnn2 = 2
    xmnn1 = amn
    xmnn2 = amn
  End If
  If (srt>(xmnn1+xmnn2)) Then
    pfinal = sqrt((s-(xmnn1+xmnn2)**2)*(s-(xmnn1-xmnn2)**2))/2./srt
    If (idxsec==1) Then
      sdmnn = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnn1+xmnn2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmnn = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
        sdmnn = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmnn = fs*pfinal/pinitial/6.
    End If
  End If
  lbnd1 = 1 + int(2*ranart(nseed))
  lbnd2 = 6 + int(4*ranart(nseed))
  If (lbnd1==1) Then
    xmnd1 = amp
  Else If (lbnd1==2) Then
    xmnd1 = amn
  End If
  xmnd2 = am0
  If (srt>(xmnd1+xmnd2)) Then
    pfinal = sqrt((s-(xmnd1+xmnd2)**2)*(s-(xmnd1-xmnd2)**2))/2./srt
    If (idxsec==1) Then
      sdmnd = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnd1+xmnd2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmnd = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmnd = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmnd = fs*pfinal/pinitial/6.
    End If
  End If
  lbns1 = 1 + int(2*ranart(nseed))
  lbns2 = 10 + int(2*ranart(nseed))
  If (lbns1==1) Then
    xmns1 = amp
  Else If (lbns1==2) Then
    xmns1 = amn
  End If
  xmns2 = am1440
  If (srt>(xmns1+xmns2)) Then
    pfinal = sqrt((s-(xmns1+xmns2)**2)*(s-(xmns1-xmns2)**2))/2./srt
    If (idxsec==1) Then
      sdmns = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmns1+xmns2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmns = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmns = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmns = fs*pfinal/pinitial/6.
    End If
  End If
  lbnp1 = 1 + int(2*ranart(nseed))
  lbnp2 = 12 + int(2*ranart(nseed))
  If (lbnp1==1) Then
    xmnp1 = amp
  Else If (lbnp1==2) Then
    xmnp1 = amn
  End If
  xmnp2 = am1535
  If (srt>(xmnp1+xmnp2)) Then
    pfinal = sqrt((s-(xmnp1+xmnp2)**2)*(s-(xmnp1-xmnp2)**2))/2./srt
    If (idxsec==1) Then
      sdmnp = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnp1+xmnp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmnp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmnp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmnp = fs*pfinal/pinitial/6.
    End If
  End If
  lbdd1 = 6 + int(4*ranart(nseed))
  lbdd2 = 6 + int(4*ranart(nseed))
  xmdd1 = am0
  xmdd2 = am0
  If (srt>(xmdd1+xmdd2)) Then
    pfinal = sqrt((s-(xmdd1+xmdd2)**2)*(s-(xmdd1-xmdd2)**2))/2./srt
    If (idxsec==1) Then
      sdmdd = fs*pfinal/pinitial*3./16.*(xnnfactor*16.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmdd1+xmdd2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmdd = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*16.)
      Else If (idxsec==4) Then
        sdmdd = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmdd = fs*pfinal/pinitial/6.
    End If
  End If
  lbds1 = 6 + int(4*ranart(nseed))
  lbds2 = 10 + int(2*ranart(nseed))
  xmds1 = am0
  xmds2 = am1440
  If (srt>(xmds1+xmds2)) Then
    pfinal = sqrt((s-(xmds1+xmds2)**2)*(s-(xmds1-xmds2)**2))/2./srt
    If (idxsec==1) Then
      sdmds = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmds1+xmds2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmds = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmds = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmds = fs*pfinal/pinitial/6.
    End If
  End If
  lbdp1 = 6 + int(4*ranart(nseed))
  lbdp2 = 12 + int(2*ranart(nseed))
  xmdp1 = am0
  xmdp2 = am1535
  If (srt>(xmdp1+xmdp2)) Then
    pfinal = sqrt((s-(xmdp1+xmdp2)**2)*(s-(xmdp1-xmdp2)**2))/2./srt
    If (idxsec==1) Then
      sdmdp = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmdp1+xmdp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmdp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmdp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmdp = fs*pfinal/pinitial/6.
    End If
  End If
  lbss1 = 10 + int(2*ranart(nseed))
  lbss2 = 10 + int(2*ranart(nseed))
  xmss1 = am1440
  xmss2 = am1440
  If (srt>(xmss1+xmss2)) Then
    pfinal = sqrt((s-(xmss1+xmss2)**2)*(s-(xmss1-xmss2)**2))/2./srt
    If (idxsec==1) Then
      sdmss = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmss1+xmss2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmss = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
        sdmss = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmns = fs*pfinal/pinitial/6.
    End If
  End If
  lbsp1 = 10 + int(2*ranart(nseed))
  lbsp2 = 12 + int(2*ranart(nseed))
  xmsp1 = am1440
  xmsp2 = am1535
  If (srt>(xmsp1+xmsp2)) Then
    pfinal = sqrt((s-(xmsp1+xmsp2)**2)*(s-(xmsp1-xmsp2)**2))/2./srt
    If (idxsec==1) Then
      sdmsp = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmsp1+xmsp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmsp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmsp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmsp = fs*pfinal/pinitial/6.
    End If
  End If
  lbpp1 = 12 + int(2*ranart(nseed))
  lbpp2 = 12 + int(2*ranart(nseed))
  xmpp1 = am1535
  xmpp2 = am1535
  If (srt>(xmpp1+xmpp2)) Then
    pfinal = sqrt((s-(xmpp1+xmpp2)**2)*(s-(xmpp1-xmpp2)**2))/2./srt
    If (idxsec==1) Then
      sdmpp = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmpp1+xmpp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmpp = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
        sdmpp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmpp = fs*pfinal/pinitial/6.
    End If
  End If
  sdm = sdmel + sdmnn + sdmnd + sdmns + sdmnp + sdmdd + sdmds + sdmdp + sdmss + sdmsp + sdmpp
  If (ianti==1) Then
    lbnn1 = -lbnn1
    lbnn2 = -lbnn2
    lbnd1 = -lbnd1
    lbnd2 = -lbnd2
    lbns1 = -lbns1
    lbns2 = -lbns2
    lbnp1 = -lbnp1
    lbnp2 = -lbnp2
    lbdd1 = -lbdd1
    lbdd2 = -lbdd2
    lbds1 = -lbds1
    lbds2 = -lbds2
    lbdp1 = -lbdp1
    lbdp2 = -lbdp2
    lbss1 = -lbss1
    lbss2 = -lbss2
    lbsp1 = -lbsp1
    lbsp2 = -lbsp2
    lbpp1 = -lbpp1
    lbpp2 = -lbpp2
  End If
  Return
End Subroutine sdmbb
