Subroutine pyovly(movly)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Dimension wti(0:100)
  Save imax, wti, wts
  If (movly==1) Then
    vint(131) = vint(106)
    If (mstp(132)>=2) vint(131) = vint(131) + vint(104)
    If (mstp(132)>=3) vint(131) = vint(131) + vint(103)
    If (mstp(132)>=4) vint(131) = vint(131) + vint(102)
    If (mstp(133)==1) Then
      xnave = vint(131)*parp(131)
      If (xnave>40.) Write (mstu(11), 1000) xnave
      wti(0) = exp(-min(50.,xnave))
      wts = 0.
      wtn = 0.
      Do i = 1, 100
        wti(i) = wti(i-1)*xnave/i
        If (i-2.5>xnave .And. wti(i)<1E-6) Goto 110
        wts = wts + wti(i)
        wtn = wtn + wti(i)*i
        imax = i
      End Do
      110 vint(132) = xnave
      vint(133) = wtn/wts
      vint(134) = wts
    Else If (mstp(133)==2) Then
      xnave = vint(131)*parp(131)
      If (xnave>40.) Write (mstu(11), 1000) xnave
      wti(1) = exp(-min(50.,xnave))*xnave
      wts = wti(1)
      wtn = wti(1)
      Do i = 2, 100
        wti(i) = wti(i-1)*xnave/(i-1)
        If (i-2.5>xnave .And. wti(i)<1E-6) Goto 130
        wts = wts + wti(i)
        wtn = wtn + wti(i)*i
        imax = i
      End Do
      130 vint(132) = xnave
      vint(133) = wtn/wts
      vint(134) = wts
    End If
  Else
    If (mstp(133)==0) Then
      mint(81) = max(1, mstp(134))
    Else
      wtr = wts*rlu(0)
      Do i = 1, imax
        mint(81) = i
        wtr = wtr - wti(i)
        If (wtr<=0.) Goto 150
      End Do
      150 Continue
    End If
  End If
  Return
  1000 Format (1X, 'Warning: requested average number of events per bunch', 'crossing too large, ', 1P, E12.4)
End Subroutine pyovly
