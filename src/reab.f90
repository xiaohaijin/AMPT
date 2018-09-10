Real Function reab(i1, i2, srt, ictrl)
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (amn=0.938, ap1=0.14, arho=0.77, aomega=0.782)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /ee/id(maxstr), lb(maxstr)
  Save
  lb1 = lb(i1)
  lb2 = lb(i2)
  reab = 0
  If (ictrl==1 .And. srt<=(amn+2.*ap1+0.02)) Return
  If (ictrl==3 .And. srt<=(amn+ap1+aomega+0.02)) Return
  pin2 = ((srt**2+ap1**2-amn**2)/(2.*srt))**2 - ap1**2
  If (pin2<=0) Return
  If (ictrl==1) Then
    If (e(i1)>1) Then
      ed = e(i1)
    Else
      ed = e(i2)
    End If
    pout2 = ((srt**2+ap1**2-ed**2)/(2.*srt))**2 - ap1**2
    If (pout2<=0) Return
    xpro = twopi(srt)/10.
    factor = 1/3.
    If (((lb1==8 .And. lb2==5) .Or. (lb1==5 .And. lb2==8)) .Or. ((lb1==-8 .And. lb2==3) .Or. (lb1==3 .And. lb2==-8))) factor = 1/4.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1.
    reab = factor*pin2/pout2*xpro
    Return
  End If
  If (ictrl==2) Then
    If (lb(i2)>=25) Then
      ed = e(i1)
      arho1 = e(i2)
    Else
      ed = e(i2)
      arho1 = e(i1)
    End If
    If (srt<=(amn+ap1+arho1+0.02)) Return
    pout2 = ((srt**2+arho1**2-ed**2)/(2.*srt))**2 - arho1**2
    If (pout2<=0) Return
    xpro = threpi(srt)/10.
    factor = 1/3.
    If (((lb1==8 .And. lb2==27) .Or. (lb1==27 .And. lb2==8)) .Or. ((lb1==-8 .And. lb2==25) .Or. (lb1==25 .And. lb2==-8))) factor = 1/4.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1.
    reab = factor*pin2/pout2*xpro
    Return
  End If
  If (ictrl==3) Then
    If (e(i1)>1) ed = e(i1)
    If (e(i2)>1) ed = e(i2)
    pout2 = ((srt**2+aomega**2-ed**2)/(2.*srt))**2 - aomega**2
    If (pout2<=0) Return
    xpro = fourpi(srt)/10.
    factor = 1/6.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1./3.
    reab = factor*pin2/pout2*xpro
  End If
  Return
End Function reab
