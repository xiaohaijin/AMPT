Subroutine xddin(px, py, pz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
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
  xinel = 0
  sigk = 0
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  xsk5 = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
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
  idd = iabs(lb(i1)*lb(i2))
  s2d = reab2d(i1, i2, srt)
  s2d = 0.
  If (((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=12)) .Or. ((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=6)) .Or. ((iabs(lb(i2))>=12) .And. (iabs(lb(i1))>=6))) Then
    xinel = sigk + s2d
    Return
  End If
  If ((idd==63) .Or. (idd==64) .Or. (idd==48) .Or. (idd==49) .Or. (idd==11*11) .Or. (idd==10*10) .Or. (idd==88) .Or. (idd==66) .Or. (idd==90) .Or. (idd==70)) Then
    xinel = x1535 + sigk + s2d
    Return
  End If
  If ((idd==110) .Or. (idd==77) .Or. (idd==80)) Then
    xinel = x1535 + sigk + s2d
    Return
  End If
  If ((idd==54) .Or. (idd==56)) Then
    sig2 = (3./4.)*sigma(srt, 2, 0, 1)
    xinel = 2.*(sig2+x1535) + sigk + s2d
    Return
  End If
  Return
End Subroutine xddin
