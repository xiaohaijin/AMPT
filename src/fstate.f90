Subroutine fstate(iseed, srt, dm3, dm4, px, py, pz, iflag)
  Dimension px(4), py(4), pz(4), pe(4)
  Common /rndf77/nseed
  Save
  iflag = -1
  pio = 3.1415926
  aka = 0.498
  icount = 0
  ekmax = (srt-dm3-dm4)/2.
  If (ekmax<=aka) Return
  pkmax = sqrt(ekmax**2-aka**2)
  If (dm3<=0.0 .Or. dm4<=0.0) Then
    Write (1, *) 'error: minus mass!!!'
    Return
  End If
  50 Continue
  icount = icount + 1
  If (icount>10) Return
  ptkmi2 = -1./4.145*alog(ranart(nseed))
  ptkm = sqrt(ptkmi2)
  3 v1 = ranart(nseed)
  v2 = ranart(nseed)
  rsq = v1**2 + v2**2
  If (rsq>=1.0 .Or. rsq<=0.) Goto 3
  fac = sqrt(-2.*alog(rsq)/rsq)
  guass = v1*fac
  If (guass>=5.) Goto 3
  xstar = guass/5.
  pzkm = pkmax*xstar
  ekm = sqrt(aka**2+pzkm**2+ptkm**2)
  If (ranart(nseed)>aka/ekm) Goto 50
  bbb = ranart(nseed)
  px(3) = ptkm*cos(2.*pio*bbb)
  py(3) = ptkm*sin(2.*pio*bbb)
  If (ranart(nseed)>0.5) pzkm = -1.*pzkm
  pz(3) = pzkm
  pe(3) = ekm
  150 ptkpl2 = -1./3.68*alog(ranart(nseed))
  ptkp = sqrt(ptkpl2)
  13 v1 = ranart(nseed)
  v2 = ranart(nseed)
  rsq = v1**2 + v2**2
  If (rsq>=1.0 .Or. rsq<=0.) Goto 13
  fac = sqrt(-2.*alog(rsq)/rsq)
  guass = v1*fac
  If (guass>=3.25) Goto 13
  xstar = guass/3.25
  pzkp = pkmax*xstar
  ekp = sqrt(aka**2+pzkp**2+ptkp**2)
  If (ranart(nseed)>aka/ekp) Goto 150
  bbb = ranart(nseed)
  px(4) = ptkp*cos(2.*pio*bbb)
  py(4) = ptkp*sin(2.*pio*bbb)
  If (ranart(nseed)>0.5) pzkp = -1.*pzkp
  pz(4) = pzkp
  pe(4) = ekp
  resten = srt - pe(3) - pe(4)
  restpz = -pz(3) - pz(4)
  If (resten<=abs(restpz)) Goto 50
  restms = sqrt(resten**2-restpz**2)
  If (restms<(dm3+dm4)) Goto 50
  ptp2 = -1./2.76*alog(ranart(nseed))
  ptp = sqrt(ptp2)
  bbb = ranart(nseed)
  px(2) = ptp*cos(2.*pio*bbb)
  py(2) = ptp*sin(2.*pio*bbb)
  px(1) = -1.*(px(4)+px(3)+px(2))
  py(1) = -1.*(py(4)+py(3)+py(2))
  rmt3 = sqrt(dm3**2+px(1)**2+py(1)**2)
  rmt4 = sqrt(dm4**2+px(2)**2+py(2)**2)
  If (restms<(rmt3+rmt4)) Goto 50
  pzcms = sqrt((restms**2-(rmt3+rmt4)**2)*(restms**2-(rmt3-rmt4)**2))/2./restms
  If (ranart(nseed)>0.5) Then
    pz(1) = pzcms
    pz(2) = -pzcms
  Else
    pz(1) = -pzcms
    pz(2) = pzcms
  End If
  beta = restpz/resten
  gama = 1./sqrt(1.-beta**2)
  pz(1) = pz(1)*gama + beta*gama*sqrt(rmt3**2+pz(1)**2)
  pz(2) = pz(2)*gama + beta*gama*sqrt(rmt4**2+pz(2)**2)
  pe(1) = sqrt(rmt3**2+pz(1)**2)
  pe(2) = sqrt(rmt4**2+pz(2)**2)
  iflag = 1
  Return
End Subroutine fstate
