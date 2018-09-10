Real Function reab2d(i1, i2, srt)
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
  reab2d = 0
  lb1 = iabs(lb(i1))
  lb2 = iabs(lb(i2))
  ed1 = e(i1)
  ed2 = e(i2)
  pin2 = (srt/2.)**2 - amn**2
  pout2 = ((srt**2+ed1**2-ed2**2)/(2.*srt))**2 - ed1**2
  If (pout2<=0) Return
  xpro = x2pi(srt)
  factor = 1/4.
  If ((lb1>=10 .And. lb1<=13) .And. (lb2>=10 .And. lb2<=13)) factor = 1.
  If ((lb1>=6 .And. lb1<=9) .And. (lb2>10 .And. lb2<=13)) factor = 1/2.
  If ((lb2>=6 .And. lb2<=9) .And. (lb1>10 .And. lb1<=13)) factor = 1/2.
  reab2d = factor*pin2/pout2*xpro
  Return
End Function reab2d
