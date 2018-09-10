Subroutine pyinre
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint6/proc(0:200)
  Character proc*28
  Save /pyint6/
  Dimension wdtp(0:40), wdte(0:40, 0:5)
  aem = paru(101)
  xw = paru(102)
  Do i = 21, 40
    Do j = 0, 40
      widp(i, j) = 0.
      wide(i, j) = 0.
    End Do
  End Do
  wmas = pmas(24, 1)
  wfac = aem/(24.*xw)*wmas
  Call pywidt(24, wmas, wdtp, wdte)
  wids(24, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(24, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(24, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(24, i) = wfac*wdtp(i)
    wide(24, i) = wfac*wdte(i, 0)
  End Do
  hcmas = pmas(37, 1)
  hcfac = aem/(8.*xw)*(hcmas/wmas)**2*hcmas
  Call pywidt(37, hcmas, wdtp, wdte)
  wids(37, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(37, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(37, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(37, i) = hcfac*wdtp(i)
    wide(37, i) = hcfac*wdte(i, 0)
  End Do
  zmas = pmas(23, 1)
  zfac = aem/(48.*xw*(1.-xw))*zmas
  Call pywidt(23, zmas, wdtp, wdte)
  wids(23, 1) = ((wdte(0,1)+wdte(0,2))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(23, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(23, 3) = 0.
  Do i = 0, 40
    widp(23, i) = zfac*wdtp(i)
    wide(23, i) = zfac*wdte(i, 0)
  End Do
  hmas = pmas(25, 1)
  hfac = aem/(8.*xw)*(hmas/wmas)**2*hmas
  Call pywidt(25, hmas, wdtp, wdte)
  wids(25, 1) = ((wdte(0,1)+wdte(0,2))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(25, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(25, 3) = 0.
  Do i = 0, 40
    widp(25, i) = hfac*wdtp(i)
    wide(25, i) = hfac*wdte(i, 0)
  End Do
  zpmas = pmas(32, 1)
  zpfac = aem/(48.*xw*(1.-xw))*zpmas
  Call pywidt(32, zpmas, wdtp, wdte)
  wids(32, 1) = ((wdte(0,1)+wdte(0,2)+wdte(0,3))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(32, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(32, 3) = 0.
  Do i = 0, 40
    widp(32, i) = zpfac*wdtp(i)
    wide(32, i) = zpfac*wdte(i, 0)
  End Do
  rmas = pmas(40, 1)
  rfac = 0.08*rmas/((mstp(1)-1)*(1.+6.*(1.+ulalps(rmas**2)/paru(1))))
  Call pywidt(40, rmas, wdtp, wdte)
  wids(40, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(40, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(40, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(40, i) = wfac*wdtp(i)
    wide(40, i) = wfac*wdte(i, 0)
  End Do
  kflqm = 1
  Do i = 1, min(8, mdcy(21,3))
    idc = i + mdcy(21, 2) - 1
    If (mdme(idc,1)<=0) Goto 170
    kflqm = i
  170 End Do
  mint(46) = kflqm
  kfpr(81, 1) = kflqm
  kfpr(81, 2) = kflqm
  kfpr(82, 1) = kflqm
  kfpr(82, 2) = kflqm
  Do i = 1, 6
    If (i<=3) kc = i + 22
    If (i==4) kc = 32
    If (i==5) kc = 37
    If (i==6) kc = 40
    pmas(kc, 2) = widp(kc, 0)
    pmas(kc, 3) = min(0.9*pmas(kc,1), 10.*pmas(kc,2))
    Do j = 1, mdcy(kc, 3)
      idc = j + mdcy(kc, 2) - 1
      brat(idc) = wide(kc, j)/wide(kc, 0)
    End Do
  End Do
  If (mstp(43)==1) Then
    proc(1) = 'f + fb -> gamma*'
  Else If (mstp(43)==2) Then
    proc(1) = 'f + fb -> Z0'
  Else If (mstp(43)==3) Then
    proc(1) = 'f + fb -> gamma*/Z0'
  End If
  If (mstp(44)==1) Then
    proc(141) = 'f + fb -> gamma*'
  Else If (mstp(44)==2) Then
    proc(141) = 'f + fb -> Z0'
  Else If (mstp(44)==3) Then
    proc(141) = 'f + fb -> Z''0'
  Else If (mstp(44)==4) Then
    proc(141) = 'f + fb -> gamma*/Z0'
  Else If (mstp(44)==5) Then
    proc(141) = 'f + fb -> gamma*/Z''0'
  Else If (mstp(44)==6) Then
    proc(141) = 'f + fb -> Z0/Z''0'
  Else If (mstp(44)==7) Then
    proc(141) = 'f + fb -> gamma*/Z0/Z''0'
  End If
  Return
End Subroutine pyinre
