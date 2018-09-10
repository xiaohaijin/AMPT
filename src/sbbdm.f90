Subroutine sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  Parameter (xmd=1.8756, ap1=0.13496, ap2=0.13957, xmrho=0.770, xmomega=0.782, xmeta=0.548, srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Common /rndf77/nseed
  Save
  sdprod = 0.
  sbbdpi = 0.
  sbbdrho = 0.
  sbbdomega = 0.
  sbbdeta = 0.
  If (srt<=(em1+em2)) Return
  ilb1 = iabs(lb1)
  ilb2 = iabs(lb2)
  s = srt**2
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<=0) Then
    Write (99, *) 'scheck50: ', scheck
    Stop
  End If
  pinitial = sqrt(scheck)/2./srt
  fs = fnndpi(s)
  If (idxsec==1 .Or. idxsec==2) Then
  Else
    If (ilb1>=1 .And. ilb1<=2 .And. ilb2>=1 .And. ilb2<=2) Then
      pifactor = 9./8.
    Else If ((ilb1>=1 .And. ilb1<=2 .And. ilb2>=6 .And. ilb2<=9) .Or. (ilb2>=1 .And. ilb2<=2 .And. ilb1>=6 .And. ilb1<=9)) Then
      pifactor = 9./64.
    Else If ((ilb1>=1 .And. ilb1<=2 .And. ilb2>=10 .And. ilb2<=13) .Or. (ilb2>=1 .And. ilb2<=2 .And. ilb1>=10 .And. ilb1<=13)) Then
      pifactor = 9./16.
    Else If (ilb1>=6 .And. ilb1<=9 .And. ilb2>=6 .And. ilb2<=9) Then
      pifactor = 9./128.
    Else If ((ilb1>=6 .And. ilb1<=9 .And. ilb2>=10 .And. ilb2<=13) .Or. (ilb2>=6 .And. ilb2<=9 .And. ilb1>=10 .And. ilb1<=13)) Then
      pifactor = 9./64.
    Else If ((ilb1>=10 .And. ilb1<=11 .And. ilb2>=10 .And. ilb2<=11) .Or. (ilb2>=12 .And. ilb2<=13 .And. ilb1>=12 .And. ilb1<=13)) Then
      pifactor = 9./8.
    Else If ((ilb1>=10 .And. ilb1<=11 .And. ilb2>=12 .And. ilb2<=13) .Or. (ilb2>=10 .And. ilb2<=11 .And. ilb1>=12 .And. ilb1<=13)) Then
      pifactor = 9./16.
    End If
  End If
  If ((ilb1*ilb2)==1) Then
    lbm = 5
    If (ianti==1) lbm = 3
    xmm = ap2
  Else If (ilb1==2 .And. ilb2==2) Then
    lbm = 3
    If (ianti==1) lbm = 5
    xmm = ap2
  Else If ((ilb1*ilb2)==2) Then
    lbm = 4
    xmm = ap1
  Else
    lbm = 3 + int(3*ranart(nseed))
    If (lbm==4) Then
      xmm = ap1
    Else
      xmm = ap2
    End If
  End If
  If (srt>=(xmd+xmm)) Then
    pfinal = sqrt((s-(xmd+xmm)**2)*(s-(xmd-xmm)**2))/2./srt
    If ((ilb1==1 .And. ilb2==1) .Or. (ilb1==2 .And. ilb2==2)) Then
      sbbdpi = fs*pfinal/pinitial/4.
    Else If ((ilb1==1 .And. ilb2==2) .Or. (ilb1==2 .And. ilb2==1)) Then
      sbbdpi = fs*pfinal/pinitial/4./2.
    Else
      If (idxsec==1) Then
        sbbdpi = fs*pfinal/pinitial*3./16.
      Else If (idxsec==2 .Or. idxsec==4) Then
        threshold = amax1(xmd+xmm, em1+em2)
        snew = (srt-threshold+srt0)**2
        If (idxsec==2) Then
          sbbdpi = fnndpi(snew)*pfinal/pinitial*3./16.
        Else If (idxsec==4) Then
          sbbdpi = fnndpi(snew)*pfinal/pinitial/6.*pifactor
        End If
      Else If (idxsec==3) Then
        sbbdpi = fs*pfinal/pinitial/6.*pifactor
      End If
    End If
  End If
  If (srt>(xmd+xmrho)) Then
    pfinal = sqrt((s-(xmd+xmrho)**2)*(s-(xmd-xmrho)**2))/2./srt
    If (idxsec==1) Then
      sbbdrho = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmrho, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdrho = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
        sbbdrho = fnndpi(snew)*pfinal/pinitial/6.*(pifactor*3.)
      End If
    Else If (idxsec==3) Then
      sbbdrho = fs*pfinal/pinitial/6.*(pifactor*3.)
    End If
  End If
  If (srt>(xmd+xmomega)) Then
    pfinal = sqrt((s-(xmd+xmomega)**2)*(s-(xmd-xmomega)**2))/2./srt
    If (idxsec==1) Then
      sbbdomega = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmomega, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdomega = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
        sbbdomega = fnndpi(snew)*pfinal/pinitial/6.*pifactor
      End If
    Else If (idxsec==3) Then
      sbbdomega = fs*pfinal/pinitial/6.*pifactor
    End If
  End If
  If (srt>(xmd+xmeta)) Then
    pfinal = sqrt((s-(xmd+xmeta)**2)*(s-(xmd-xmeta)**2))/2./srt
    If (idxsec==1) Then
      sbbdeta = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmeta, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdeta = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
        sbbdeta = fnndpi(snew)*pfinal/pinitial/6.*(pifactor/3.)
      End If
    Else If (idxsec==3) Then
      sbbdeta = fs*pfinal/pinitial/6.*(pifactor/3.)
    End If
  End If
  sdprod = sbbdpi + sbbdrho + sbbdomega + sbbdeta
  If (sdprod<=0) Return
  x1 = ranart(nseed)
  If (x1<=sbbdpi/sdprod) Then
  Else If (x1<=(sbbdpi+sbbdrho)/sdprod) Then
    lbm = 25 + int(3*ranart(nseed))
    xmm = xmrho
  Else If (x1<=(sbbdpi+sbbdrho+sbbdomega)/sdprod) Then
    lbm = 28
    xmm = xmomega
  Else
    lbm = 0
    xmm = xmeta
  End If
  Return
End Subroutine sbbdm
