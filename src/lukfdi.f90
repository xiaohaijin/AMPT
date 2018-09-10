Subroutine lukfdi(kfl1, kfl2, kfl3, kf)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  kf1a = iabs(kfl1)
  kf2a = iabs(kfl2)
  kfl3 = 0
  kf = 0
  If (kf1a==0) Return
  If (kf2a/=0) Then
    If (kf1a<=10 .And. kf2a<=10 .And. kfl1*kfl2>0) Return
    If (kf1a>10 .And. kf2a>10) Return
    If ((kf1a>10 .Or. kf2a>10) .And. kfl1*kfl2<0) Return
  End If
  If (mstj(15)==1) Then
    ktab1 = -1
    If (kf1a>=1 .And. kf1a<=6) ktab1 = kf1a
    kfl1a = mod(kf1a/1000, 10)
    kfl1b = mod(kf1a/100, 10)
    kfl1s = mod(kf1a, 10)
    If (kfl1a>=1 .And. kfl1a<=4 .And. kfl1b>=1 .And. kfl1b<=4) ktab1 = 6 + kfl1a*(kfl1a-2) + 2*kfl1b + (kfl1s-1)/2
    If (kfl1a>=1 .And. kfl1a<=4 .And. kfl1a==kfl1b) ktab1 = ktab1 - 1
    If (kf1a>=1 .And. kf1a<=6) kfl1a = kf1a
    ktab2 = 0
    If (kf2a/=0) Then
      ktab2 = -1
      If (kf2a>=1 .And. kf2a<=6) ktab2 = kf2a
      kfl2a = mod(kf2a/1000, 10)
      kfl2b = mod(kf2a/100, 10)
      kfl2s = mod(kf2a, 10)
      If (kfl2a>=1 .And. kfl2a<=4 .And. kfl2b>=1 .And. kfl2b<=4) ktab2 = 6 + kfl2a*(kfl2a-2) + 2*kfl2b + (kfl2s-1)/2
      If (kfl2a>=1 .And. kfl2a<=4 .And. kfl2a==kfl2b) ktab2 = ktab2 - 1
    End If
    If (ktab1>=0 .And. ktab2>=0) Goto 140
  End If
  100 par2 = parj(2)
  par3 = parj(3)
  par4 = 3.*parj(4)
  If (mstj(12)>=2) Then
    par3m = sqrt(parj(3))
    par4m = 1./(3.*sqrt(parj(4)))
    pardm = parj(7)/(parj(7)+par3m*parj(6))
    pars0 = parj(5)*(2.+(1.+par2*par3m*parj(7))*(1.+par4m))
    pars1 = parj(7)*pars0/(2.*par3m) + parj(5)*(parj(6)*(1.+par4m)+par2*par3m*parj(6)*parj(7))
    pars2 = parj(5)*2.*parj(6)*parj(7)*(par2*parj(7)+(1.+par4m)/par3m)
    parsm = max(pars0, pars1, pars2)
    par4 = par4*(1.+parsm)/(1.+parsm/(3.*par4m))
  End If
  mbary = 0
  kfda = 0
  If (kf1a<=10) Then
    If (kf2a==0 .And. mstj(12)>=1 .And. (1.+parj(1))*rlu(0)>1.) mbary = 1
    If (kf2a>10) mbary = 2
    If (kf2a>10 .And. kf2a<=10000) kfda = kf2a
  Else
    mbary = 2
    If (kf1a<=10000) kfda = kf1a
  End If
  If (kfda/=0 .And. mstj(12)>=2) Then
    kflda = mod(kfda/1000, 10)
    kfldb = mod(kfda/100, 10)
    kflds = mod(kfda, 10)
    wtdq = pars0
    If (max(kflda,kfldb)==3) wtdq = pars1
    If (min(kflda,kfldb)==3) wtdq = pars2
    If (kflds==1) wtdq = wtdq/(3.*par4m)
    If ((1.+wtdq)*rlu(0)>1.) mbary = -1
    If (mbary==-1 .And. kf2a/=0) Return
  End If
  If (mbary<=0) Then
    kfs = isign(1, kfl1)
    If (mbary==0) Then
      If (kf2a==0) kfl3 = isign(1+int((2.+par2)*rlu(0)), -kfl1)
      kfla = max(kf1a, kf2a+iabs(kfl3))
      kflb = min(kf1a, kf2a+iabs(kfl3))
      If (kfla/=kf1a) kfs = -kfs
    Else
      kfl1a = mod(kf1a/1000, 10)
      kfl1b = mod(kf1a/100, 10)
      110 kfl1d = kfl1a + int(rlu(0)+0.5)*(kfl1b-kfl1a)
      kfl1e = kfl1a + kfl1b - kfl1d
      If ((kfl1d==3 .And. rlu(0)>pardm) .Or. (kfl1e==3 .And. rlu(0)<pardm)) Then
        kfl1d = kfl1a + kfl1b - kfl1d
        kfl1e = kfl1a + kfl1b - kfl1e
      End If
      kfl3a = 1 + int((2.+par2*par3m*parj(7))*rlu(0))
      If ((kfl1e/=kfl3a .And. rlu(0)>(1.+par4m)/max(2.,1.+par4m)) .Or. (kfl1e==kfl3a .And. rlu(0)>2./max(2.,1.+par4m))) Goto 110
      kflds = 3
      If (kfl1e/=kfl3a) kflds = 2*int(rlu(0)+1./(1.+par4m)) + 1
      kfl3 = isign(10000+1000*max(kfl1e,kfl3a)+100*min(kfl1e,kfl3a)+kflds, -kfl1)
      kfla = max(kfl1d, kfl3a)
      kflb = min(kfl1d, kfl3a)
      If (kfla/=kfl1d) kfs = -kfs
    End If
    If (kfla<=2) kmul = int(parj(11)+rlu(0))
    If (kfla==3) kmul = int(parj(12)+rlu(0))
    If (kfla>=4) kmul = int(parj(13)+rlu(0))
    If (kmul==0 .And. parj(14)>0.) Then
      If (rlu(0)<parj(14)) kmul = 2
    Else If (kmul==1 .And. parj(15)+parj(16)+parj(17)>0.) Then
      rmul = rlu(0)
      If (rmul<parj(15)) kmul = 3
      If (kmul==1 .And. rmul<parj(15)+parj(16)) kmul = 4
      If (kmul==1 .And. rmul<parj(15)+parj(16)+parj(17)) kmul = 5
    End If
    kfls = 3
    If (kmul==0 .Or. kmul==3) kfls = 1
    If (kmul==5) kfls = 5
    If (kfla/=kflb) Then
      kf = (100*kfla+10*kflb+kfls)*kfs*(-1)**kfla
    Else
      rmix = rlu(0)
      imix = 2*kfla + 10*kmul
      If (kfla<=3) kf = 110*(1+int(rmix+parf(imix-1))+int(rmix+parf(imix))) + kfls
      If (kfla>=4) kf = 110*kfla + kfls
    End If
    If (kmul==2 .Or. kmul==3) kf = kf + isign(10000, kf)
    If (kmul==4) kf = kf + isign(20000, kf)
  Else
    120 If (kf1a<=10 .And. kf2a==0) Then
      kfla = kf1a
      130 kflb = 1 + int((2.+par2*par3)*rlu(0))
      kflc = 1 + int((2.+par2*par3)*rlu(0))
      kflds = 1
      If (kflb>=kflc) kflds = 3
      If (kflds==1 .And. par4*rlu(0)>1.) Goto 130
      If (kflds==3 .And. par4<rlu(0)) Goto 130
      kfl3 = isign(1000*max(kflb,kflc)+100*min(kflb,kflc)+kflds, kfl1)
    Else If (kf1a<=10) Then
      kfla = kf1a
      kflb = mod(kf2a/1000, 10)
      kflc = mod(kf2a/100, 10)
      kflds = mod(kf2a, 10)
    Else
      If (kf2a==0) kfl3 = isign(1+int((2.+par2)*rlu(0)), kfl1)
      kfla = kf2a + iabs(kfl3)
      kflb = mod(kf1a/1000, 10)
      kflc = mod(kf1a/100, 10)
      kflds = mod(kf1a, 10)
    End If
    kbary = kflds
    If (kflds==3 .And. kflb/=kflc) kbary = 5
    If (kfla/=kflb .And. kfla/=kflc) kbary = kbary + 1
    wt = parf(60+kbary) + parj(18)*parf(70+kbary)
    If (mbary==1 .And. mstj(12)>=2) Then
      wtdq = pars0
      If (max(kflb,kflc)==3) wtdq = pars1
      If (min(kflb,kflc)==3) wtdq = pars2
      If (kflds==1) wtdq = wtdq/(3.*par4m)
      If (kflds==1) wt = wt*(1.+wtdq)/(1.+parsm/(3.*par4m))
      If (kflds==3) wt = wt*(1.+wtdq)/(1.+parsm)
    End If
    If (kf2a==0 .And. wt<rlu(0)) Goto 120
    kfld = max(kfla, kflb, kflc)
    kflf = min(kfla, kflb, kflc)
    kfle = kfla + kflb + kflc - kfld - kflf
    kfls = 2
    If ((parf(60+kbary)+parj(18)*parf(70+kbary))*rlu(0)>parf(60+kbary)) kfls = 4
    kfll = 0
    If (kfls==2 .And. kfld>kfle .And. kfle>kflf) Then
      If (kflds==1 .And. kfla==kfld) kfll = 1
      If (kflds==1 .And. kfla/=kfld) kfll = int(0.25+rlu(0))
      If (kflds==3 .And. kfla/=kfld) kfll = int(0.75+rlu(0))
    End If
    If (kfll==0) kf = isign(1000*kfld+100*kfle+10*kflf+kfls, kfl1)
    If (kfll==1) kf = isign(1000*kfld+100*kflf+10*kfle+kfls, kfl1)
  End If
  Return
  140 If (ktab2==0 .And. mstj(12)<=0) Then
    kt3l = 1
    kt3u = 6
  Else If (ktab2==0 .And. ktab1>=7 .And. mstj(12)<=1) Then
    kt3l = 1
    kt3u = 6
  Else If (ktab2==0) Then
    kt3l = 1
    kt3u = 22
  Else
    kt3l = ktab2
    kt3u = ktab2
  End If
  rfl = 0.
  Do kts = 0, 2
    Do kt3 = kt3l, kt3u
      rfl = rfl + parf(120+80*ktab1+25*kts+kt3)
    End Do
  End Do
  rfl = rlu(0)*rfl
  Do kts = 0, 2
    ktabs = kts
    Do kt3 = kt3l, kt3u
      ktab3 = kt3
      rfl = rfl - parf(120+80*ktab1+25*kts+kt3)
      If (rfl<=0.) Goto 170
    End Do
  End Do
  170 Continue
  If (ktab3<=6) Then
    kfl3a = ktab3
    kfl3b = 0
    kfl3 = isign(kfl3a, kfl1*(2*ktab1-13))
  Else
    kfl3a = 1
    If (ktab3>=8) kfl3a = 2
    If (ktab3>=11) kfl3a = 3
    If (ktab3>=16) kfl3a = 4
    kfl3b = (ktab3-6-kfl3a*(kfl3a-2))/2
    kfl3 = 1000*kfl3a + 100*kfl3b + 1
    If (kfl3a==kfl3b .Or. ktab3/=6+kfl3a*(kfl3a-2)+2*kfl3b) kfl3 = kfl3 + 2
    kfl3 = isign(kfl3, kfl1*(13-2*ktab1))
  End If
  If (kfl3a==kfl1a .And. kfl3b==kfl1b .And. (kfl3a<=3 .Or. kfl3b/=0)) Then
    rfl = rlu(0)*(parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+25*ktabs)+parf(145+80*ktab1+25*ktabs))
    kf = 110 + 2*ktabs + 1
    If (rfl>parf(143+80*ktab1+25*ktabs)) kf = 220 + 2*ktabs + 1
    If (rfl>parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+25*ktabs)) kf = 330 + 2*ktabs + 1
  Else If (ktab1<=6 .And. ktab3<=6) Then
    kfla = max(ktab1, ktab3)
    kflb = min(ktab1, ktab3)
    kfs = isign(1, kfl1)
    If (kfla/=kf1a) kfs = -kfs
    kf = (100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla
  Else If (ktab1>=7 .And. ktab3>=7) Then
    kfs = isign(1, kfl1)
    If (kfl1a==kfl3a) Then
      kfla = max(kfl1b, kfl3b)
      kflb = min(kfl1b, kfl3b)
      If (kfla/=kfl1b) kfs = -kfs
    Else If (kfl1a==kfl3b) Then
      kfla = kfl3a
      kflb = kfl1b
      kfs = -kfs
    Else If (kfl1b==kfl3a) Then
      kfla = kfl1a
      kflb = kfl3b
    Else If (kfl1b==kfl3b) Then
      kfla = max(kfl1a, kfl3a)
      kflb = min(kfl1a, kfl3a)
      If (kfla/=kfl1a) kfs = -kfs
    Else
      Call luerrm(2, '(LUKFDI:) no matching flavours for qq -> qq')
      Goto 100
    End If
    kf = (100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla
  Else
    If (ktab1>=7) Then
      kfla = kfl3a
      kflb = kfl1a
      kflc = kfl1b
    Else
      kfla = kfl1a
      kflb = kfl3a
      kflc = kfl3b
    End If
    kfld = max(kfla, kflb, kflc)
    kflf = min(kfla, kflb, kflc)
    kfle = kfla + kflb + kflc - kfld - kflf
    If (ktabs==0) kf = isign(1000*kfld+100*kflf+10*kfle+2, kfl1)
    If (ktabs>=1) kf = isign(1000*kfld+100*kfle+10*kflf+2*ktabs, kfl1)
  End If
  If (kfl2/=0) kfl3 = 0
  kc = lucomp(kf)
  If (kc==0) Then
    Call luerrm(2, '(LUKFDI:) user-defined flavour probabilities '//'failed')
    Goto 100
  End If
  Return
End Subroutine lukfdi
