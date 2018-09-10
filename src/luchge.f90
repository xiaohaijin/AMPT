Function luchge(kf)
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  luchge = 0
  kfa = iabs(kf)
  kc = lucomp(kfa)
  If (kc==0) Then
  Else If (kfa<=100 .Or. kc<=80 .Or. kc>100) Then
    luchge = kchg(kc, 1)
  Else If (mod(kfa/1000,10)==0) Then
    luchge = (kchg(mod(kfa/100,10),1)-kchg(mod(kfa/10,10),1))*(-1)**mod(kfa/100, 10)
  Else If (mod(kfa/10,10)==0) Then
    luchge = kchg(mod(kfa/1000,10), 1) + kchg(mod(kfa/100,10), 1)
  Else
    luchge = kchg(mod(kfa/1000,10), 1) + kchg(mod(kfa/100,10), 1) + kchg(mod(kfa/10,10), 1)
  End If
  luchge = luchge*isign(1, kf)
  Return
End Function luchge
