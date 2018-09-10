Subroutine xnd(px, py, pz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /input/nstar, ndirct, dir
  Common /nn/nnn
  Common /bg/betax, betay, betaz, gamma
  Common /run/num
  Common /pa/rpion(3, maxstr, maxr)
  Common /pb/ppion(3, maxstr, maxr)
  Common /pc/epion(maxstr, maxr)
  Common /pd/lpion(maxstr, maxr)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Save
  xinel = 0.
  sigk = 0
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  xsk5 = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  If (srt<2.04) Return
  prf = sqrt(0.25*srt**2-avmass**2)
  If (em1>1.) Then
    deltam = em1
  Else
    deltam = em2
  End If
  renom = deltam*prf**2/denom(srt, 1.)/pr
  renomn = deltam*prf**2/denom(srt, 2.)/pr
  renom1 = deltam*prf**2/denom(srt, -1.)/pr
  If ((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) renom = 0.
  If ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) renom = 0.
  If ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) renom = 0.
  If ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9)) renom = 0.
  Call m1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
  x1440 = (3./4.)*sigma(srt, 2, 0, 1)
  akp = 0.498
  ak0 = 0.498
  ana = 0.94
  ada = 1.232
  al = 1.1157
  as = 1.1197
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  xsk5 = 0
  t1nlk = ana + al + akp
  If (srt<=t1nlk) Goto 222
  xsk1 = 1.5*pplpk(srt)
  t1dlk = ada + al + akp
  t2dlk = ada + al - akp
  If (srt<=t1dlk) Goto 222
  es = srt
  pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
  pmdlk = sqrt(pmdlk2)
  xsk3 = 1.5*pplpk(srt)
  t1nsk = ana + as + akp
  t2nsk = ana + as - akp
  If (srt<=t1nsk) Goto 222
  pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
  pmnsk = sqrt(pmnsk2)
  xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
  t1dsk = ada + as + akp
  t2dsk = ada + as - akp
  If (srt<=t1dsk) Goto 222
  pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
  pmdsk = sqrt(pmdsk2)
  xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
  If (srt<=(2.*amn+aphi)) Goto 222
  xsk5 = 0.0001
  222 sigk = xsk1 + xsk2 + xsk3 + xsk4
  xsk1 = 2.0*xsk1
  xsk2 = 2.0*xsk2
  xsk3 = 2.0*xsk3
  xsk4 = 2.0*xsk4
  sigk = 2.0*sigk + xsk5
  If (((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) .Or. ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) .Or. ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) .Or. ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9))) Then
    xinel = sigk
    Return
  End If
  If (lb(i1)*lb(i2)==18 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==6 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==8 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = 1.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==14 .And. (iabs(lb(i1))==2 .And. iabs(lb(i2))==2)) Then
    signd = 1.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==16 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
    sigdn = 0.5*signd*renom
    xinel = sigdn + 2.*x1440 + 2.*x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==7) Then
    signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
    sigdn = 0.5*signd*renom
    xinel = sigdn + 2.*x1440 + 2.*x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==10 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = (3./4.)*sigma(srt, 2, 0, 1)
    sigdn = signd*renomn
    xinel = sigdn + x1535 + sigk
    Return
  End If
  If (lb(i1)*lb(i2)==22 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = (3./4.)*sigma(srt, 2, 0, 1)
    sigdn = signd*renomn
    xinel = sigdn + x1535 + sigk
    Return
  End If
  If ((iabs(lb(i1))==12) .Or. (iabs(lb(i1))==13) .Or. (iabs(lb(i2))==12) .Or. (iabs(lb(i2))==13)) Then
    signd = x1535
    sigdn = signd*renom1
    xinel = sigdn + sigk
    Return
  End If
  Return
End Subroutine xnd
