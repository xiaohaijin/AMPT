Real Function x4pi(srt)
  Save
  akp = 0.498
  ak0 = 0.498
  ana = 0.94
  ada = 1.232
  al = 1.1157
  as = 1.1197
  pmass = 0.9383
  es = srt
  If (es<=4) Then
     x4pi = 0.
  Else
     xpp2pi = 4.*x2pi(es)
     xpp3pi = 3.*(x3pi(es)+x33pi(es))
     pps1 = sigma(es, 1, 1, 0) + 0.5*sigma(es, 1, 1, 1)
     pps2 = 1.5*sigma(es, 1, 1, 1)
     ppsngl = pps1 + pps2 + s1535(es)
     xk1 = 0
     xk2 = 0
     xk3 = 0
     xk4 = 0
     t1nlk = ana + al + akp
     t2nlk = ana + al - akp
     If (es<=t1nlk) Goto 333
     pmnlk2 = (es**2-t1nlk**2)*(es**2-t2nlk**2)/(4.*es**2)
     pmnlk = sqrt(pmnlk2)
     xk1 = pplpk(es)
     t1dlk = ada + al + akp
     t2dlk = ada + al - akp
     If (es<=t1dlk) Goto 333
     pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
     pmdlk = sqrt(pmdlk2)
     xk3 = pplpk(es)
     t1nsk = ana + as + akp
     t2nsk = ana + as - akp
     If (es<=t1nsk) Goto 333
     pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
     pmnsk = sqrt(pmnsk2)
     xk2 = ppk1(es) + ppk0(es)
     t1dsk = ada + as + akp
     t2dsk = ada + as - akp
     If (es<=t1dsk) Goto 333
     pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
     pmdsk = sqrt(pmdsk2)
     xk4 = ppk1(es) + ppk0(es)
333  xkaon = 3.*(xk1+xk2+xk3+xk4)
     x4pi = pp1(es) - ppsngl - xpp2pi - xpp3pi - xkaon
     If (x4pi<=0) x4pi = 1.E-06
  End If
  Return
End Function x4pi
