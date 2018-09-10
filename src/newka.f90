Subroutine newka(icase, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, iblock)
  Parameter (maxstr=150001, maxr=1)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /nn/nnn
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /rndf77/nseed
  Save
  Logical lb1bn, lb2bn, lb1mn, lb2mn
  Logical lb1bn1, lb2bn1, lb1bn0, lb2bn0
  Logical lb1mn0, lb2mn0, lb1mn1, lb2mn1
  Logical lb1mn2, lb2mn2
  icase = -1
  nchrg = -100
  ictrl = 1
  lb1 = lb(i1)
  lb2 = lb(i2)
  em1 = e(i1)
  em2 = e(i2)
  lb1bn = lb1 == 1 .Or. lb1 == 2 .Or. (lb1>5 .And. lb1<=13)
  lb2bn = lb2 == 1 .Or. lb2 == 2 .Or. (lb2>5 .And. lb2<=13)
  lb1bn0 = lb1 == 2 .Or. lb1 == 7 .Or. lb1 == 10 .Or. lb1 == 12
  lb2bn0 = lb2 == 2 .Or. lb2 == 7 .Or. lb2 == 10 .Or. lb2 == 12
  lb1bn1 = lb1 == 1 .Or. lb1 == 8 .Or. lb1 == 11 .Or. lb1 == 13
  lb2bn1 = lb2 == 1 .Or. lb2 == 8 .Or. lb2 == 11 .Or. lb2 == 13
  lb1mn = em1 < 0.2 .Or. lb1 == 0 .Or. (lb1>=25 .And. lb1<=29)
  lb2mn = em2 < 0.2 .Or. lb2 == 0 .Or. (lb2>=25 .And. lb2<=29)
  lb1mn0 = lb1 == 0 .Or. lb1 == 4 .Or. lb1 == 26 .Or. lb1 == 28 .Or. lb1 == 29
  lb2mn0 = lb2 == 0 .Or. lb2 == 4 .Or. lb2 == 26 .Or. lb2 == 28 .Or. lb2 == 29
  lb1mn1 = lb1 == 5 .Or. lb1 == 27
  lb2mn1 = lb2 == 5 .Or. lb2 == 27
  lb1mn2 = lb1 == 3 .Or. lb1 == 25
  lb2mn2 = lb2 == 3 .Or. lb2 == 25
  If (lb1bn .And. lb2bn) Then
    icase = 1
    sig = 40.
    If (lb1==9 .And. lb2==9) Then
      nchrg = 4
    End If
    If ((lb1bn1 .And. lb2==9) .Or. (lb2bn1 .And. lb1==9)) Then
      nchrg = 3
    End If
    If ((lb1bn0 .And. lb2==9) .Or. (lb2bn0 .And. lb1==9) .Or. (lb1bn1 .And. lb2bn1)) Then
      nchrg = 2
    End If
    If ((lb1bn1 .And. lb2bn0) .Or. (lb1==6 .And. lb2==9) .Or. (lb2bn1 .And. lb1bn0) .Or. (lb2==6 .And. lb1==9)) Then
      nchrg = 1
    End If
    If ((lb1bn0 .And. lb2bn0) .Or. (lb1bn1 .And. lb2==6) .Or. (lb2bn1 .And. lb1==6)) Then
      nchrg = 0
    End If
    If ((lb1bn0 .And. lb2==6) .Or. (lb2bn0 .And. lb1==6)) Then
      nchrg = -1
    End If
    If (lb1==6 .And. lb2==6) Then
      nchrg = -2
    End If
    If (nchrg>=-1 .And. nchrg<=2) Then
      brsig = x2kaon(srt)
    Else
      brsig = 0.0
    End If
    brsig = 2.0*brsig
  End If
  If ((lb1bn .And. lb2mn) .Or. (lb2bn .And. lb1mn)) Then
    icase = 2
    sig = 20.
    sigma0 = pinsg0(srt)
    brsig = 0.0
    If ((lb1bn1 .And. lb2mn0) .Or. (lb2bn1 .And. lb1mn0) .Or. (lb1bn0 .And. lb2mn1) .Or. (lb2bn0 .And. lb1mn1) .Or. (lb1==9 .And. lb2mn2) .Or. (lb2==9 .And. lb1mn2)) Then
      nchrg = 1
      If (lb1bn1 .Or. lb2bn1) brsig = 0.5*sigma0
      If (lb1bn0 .Or. lb2bn0) brsig = 2.0*sigma0
    End If
    If ((lb1bn0 .And. lb2mn0) .Or. (lb2bn0 .And. lb1mn0) .Or. (lb1bn1 .And. lb2mn2) .Or. (lb2bn1 .And. lb1mn2) .Or. (lb1==6 .And. lb2mn1) .Or. (lb2==6 .And. lb1mn1)) Then
      nchrg = 0
      If (lb1bn1 .Or. lb2bn1) Then
        brsig = 3.0*sigma0
        ratiok = 2./3.
      End If
      If (lb1bn0 .Or. lb2bn0) Then
        brsig = 2.5*sigma0
        ratiok = 0.2
      End If
    End If
    If ((lb1bn0 .And. lb2mn2) .Or. (lb2bn0 .And. lb1mn2) .Or. (lb1==6 .And. lb2mn0) .Or. (lb2==6 .And. lb1mn0)) Then
      nchrg = -1
      If (lb1bn0 .Or. lb2bn0) brsig = sigma0
    End If
    If ((lb1==6 .And. lb2mn2) .Or. (lb2==6 .And. lb1mn2)) Then
      nchrg = -2
    End If
    If ((lb1bn1 .And. lb2mn1) .Or. (lb2bn1 .And. lb1mn1) .Or. (lb1==9 .And. lb2mn0) .Or. (lb2==9 .And. lb1mn0)) Then
      nchrg = 2
    End If
    If (nchrg>=-2 .And. nchrg<=2) Then
      brsig = 3.0*sigma0
    End If
  End If
  If ((lb1bn .And. (lb2==21 .Or. lb2==-30)) .Or. (lb2bn .And. (lb1==21 .Or. lb1==-30))) Then
    bmass = 0.938
    If (srt<=(bmass+aka)) Then
      pkaon = 0.
    Else
      pkaon = sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
    End If
    sig = 0.
    If (lb1==1 .Or. lb2==1 .Or. lb1==8 .Or. lb2==8 .Or. lb1==11 .Or. lb2==11 .Or. lb1==13 .Or. lb2==13) Then
      nchrg = 0
      sigela = akpel(pkaon)
      sigsgm = 3.*akpsgm(pkaon)
      sig = sigela + sigsgm + akplam(pkaon)
    End If
    If (lb1==2 .Or. lb2==2 .Or. lb1==7 .Or. lb2==7 .Or. lb1==10 .Or. lb2==10 .Or. lb1==12 .Or. lb2==12) Then
      nchrg = -1
      sigela = aknel(pkaon)
      sigsgm = 2.*aknsgm(pkaon)
      sig = sigela + sigsgm + aknlam(pkaon)
    End If
    If (lb1==6 .Or. lb2==6) Then
      nchrg = -2
      sigela = aknel(pkaon)
      sigsgm = aknsgm(pkaon)
      sig = sigela + sigsgm
    End If
    If (lb1==9 .Or. lb2==9) Then
      nchrg = 1
      sigela = akpel(pkaon)
      sigsgm = 2.*akpsgm(pkaon)
      sig = sigela + sigsgm + akplam(pkaon)
    End If
    sigela = 0.5*(akpel(pkaon)+aknel(pkaon))
    sigsgm = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
    sig = sigela + sigsgm + akplam(pkaon)
    If (sig>1.E-7) Then
      icase = 3
      brel = sigela/sig
      brsgm = sigsgm/sig
      brsig = sig
    End If
  End If
  If (((lb1>=14 .And. lb1<=17) .And. (lb2>=3 .And. lb2<=5)) .Or. ((lb2>=14 .And. lb2<=17) .And. (lb1>=3 .And. lb1<=5))) Then
    nchrg = -100
    If ((lb1==15 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==15 .And. (lb1==3 .Or. lb1==25))) Then
      nchrg = -2
      bmass = 1.232
    End If
    If ((lb1==15 .And. lb2mn0) .Or. (lb2==15 .And. lb1mn0) .Or. ((lb1==14 .Or. lb1==16) .And. (lb2==3 .Or. lb2==25)) .Or. ((lb2==14 .Or. lb2==16) .And. (lb1==3 .Or. lb1==25))) Then
      nchrg = -1
      bmass = 0.938
    End If
    If ((lb1==15 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==15 .And. (lb1==5 .Or. lb1==27)) .Or. (lb1==17 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==17 .And. (lb1==3 .Or. lb1==25)) .Or. ((lb1==14 .Or. lb1==16) .And. lb2mn0) .Or. ((lb2==14 .Or. lb2==16) .And. lb1mn0)) Then
      nchrg = 0
      bmass = 0.938
    End If
    If ((lb1==17 .And. lb2mn0) .Or. (lb2==17 .And. lb1mn0) .Or. ((lb1==14 .Or. lb1==16) .And. (lb2==5 .Or. lb2==27)) .Or. ((lb2==14 .Or. lb2==16) .And. (lb1==5 .Or. lb1==27))) Then
      nchrg = 1
      bmass = 1.232
    End If
    sig = 0.
    If (nchrg/=-100 .And. srt>(aka+bmass)) Then
      icase = 4
      pkaon = sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
      If (lb1==14 .Or. lb2==14) Then
        If (nchrg>=0) sigma0 = akplam(pkaon)
        If (nchrg<0) sigma0 = aknlam(pkaon)
      Else
        If (nchrg>=0) sigma0 = akpsgm(pkaon)
        If (nchrg<0) sigma0 = aknsgm(pkaon)
        sigma0 = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
      End If
      sig = (srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
      If (nchrg==-2 .Or. nchrg==2) sig = 2.*sig
      If (lb1==14 .Or. lb2==14) Then
        sig = 4.0/3.0*sig
      Else If (nchrg==-2 .Or. nchrg==2) Then
        sig = 8.0/9.0*sig
      Else
        sig = 4.0/9.0*sig
      End If
      brsig = sig
      If (sig<1.E-7) sig = 1.E-7
    End If
    icase = 4
    brsig = sig
    sigela = 10.
    sig = sig + sigela
    brel = sigela/sig
  End If
  If (icase==-1) Then
    ictrl = -1
    Return
  End If
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
  If (ic==-1) Then
    ictrl = -1
    Return
  End If
  ik = 0
  ik0 = 0
  ik1 = 0
  ik2 = 0
  ik3 = 0
  il = 0
  im = 0
  im3 = 0
  im4 = 0
  in = 0
  inpion = 0
  ipipi = 0
  sgsum = 0.
  sgsum1 = 0.
  sgsum3 = 0.
  If (icase==1) Then
    ik = ik + 1
    If (srt>2.8639) Then
      ik0 = ik0 + 1
      If (em1<1.0 .And. em2<1.0) Then
        ik1 = ik1 + 1
        sgsum1 = sgsum1 + brsig
      End If
      If (em1>1.0 .And. em2>1.0) Then
        ik3 = ik3 + 1
        sgsum3 = sgsum3 + brsig
      End If
      If (em1>1.0 .And. em2<1.0) ik2 = ik2 + 1
      If (em1<1.0 .And. em2>1.0) ik2 = ik2 + 1
      sgsum = sgsum + brsig
    End If
  End If
  If (icase==2) inpion = inpion + 1
  If (icase==5) ipipi = ipipi + 1
  If (ranart(nseed)>(brsig/sig)) Then
    ictrl = -1
    Return
  End If
  il = il + 1
  If (icase==1) Then
    in = in + 1
    Call nnkaon(irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
  End If
  If (icase==2) Then
    im = im + 1
    Call npik(irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, ratiok, iblock)
  End If
  If (icase==3) Then
    im3 = im3 + 1
    Call kaonn(brel, brsgm, irun, iseed, dt, nt, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
  End If
  If (icase==4) Then
    im4 = im4 + 1
    If (ranart(nseed)<brel) Then
      ielstc = 1
    Else
      ielstc = 0
    End If
    Call pihypn(ielstc, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, iblock)
  End If
  Return
End Subroutine newka
