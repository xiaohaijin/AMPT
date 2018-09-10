Subroutine lambar(i1, i2, srt, siglab)
  Parameter (maxstr=150001)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Save
  siglab = 1.E-06
  If (iabs(lb(i1))>=14 .And. iabs(lb(i1))<=17) Then
    eml = e(i1)
    emb = e(i2)
  Else
    eml = e(i2)
    emb = e(i1)
  End If
  pthr = srt**2 - eml**2 - emb**2
  If (pthr>0.) Then
    plab2 = (pthr/2./emb)**2 - eml**2
    If (plab2>0) Then
      plab = sqrt(plab2)
      siglab = 12. + 0.43/(plab**3.3)
      If (siglab>200.) siglab = 200.
    End If
  End If
  Return
End Subroutine lambar
