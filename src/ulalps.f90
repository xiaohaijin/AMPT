Function ulalps(q2)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  If (mstu(111)<=0) Then
    ulalps = paru(111)
    mstu(118) = mstu(112)
    paru(117) = 0.
    paru(118) = paru(111)
    Return
  End If
  q2eff = q2
  If (mstu(115)>=2) q2eff = max(q2, paru(114))
  nf = mstu(112)
  alam2 = paru(112)**2
  100 If (nf>max(2,mstu(113))) Then
    q2thr = paru(113)*pmas(nf, 1)**2
    If (q2eff<q2thr) Then
      nf = nf - 1
      alam2 = alam2*(q2thr/alam2)**(2./(33.-2.*nf))
      Goto 100
    End If
  End If
  110 If (nf<min(8,mstu(114))) Then
    q2thr = paru(113)*pmas(nf+1, 1)**2
    If (q2eff>q2thr) Then
      nf = nf + 1
      alam2 = alam2*(alam2/q2thr)**(2./(33.-2.*nf))
      Goto 110
    End If
  End If
  If (mstu(115)==1) q2eff = q2eff + alam2
  paru(117) = sqrt(alam2)
  b0 = (33.-2.*nf)/6.
  algq = log(q2eff/alam2)
  If (mstu(111)==1) Then
    ulalps = paru(2)/(b0*algq)
  Else
    b1 = (153.-19.*nf)/6.
    ulalps = paru(2)/(b0*algq)*(1.-b1*log(algq)/(b0**2*algq))
  End If
  mstu(118) = nf
  paru(118) = ulalps
  Return
End Function ulalps
