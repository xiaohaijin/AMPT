Subroutine pyfram(iframe)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  If (iframe<1 .Or. iframe>2) Then
    Write (mstu(11), 1000) iframe, mint(6)
    Return
  End If
  If (iframe==mint(6)) Return
  If (mint(6)==1) Then
    Call lurobo(0., 0., -vint(8), -vint(9), -vint(10))
    Call lurobo(0., -vint(7), 0., 0., 0.)
    Call lurobo(-vint(6), 0., 0., 0., 0.)
    mint(6) = 2
  Else
    Call lurobo(vint(6), vint(7), vint(8), vint(9), vint(10))
    mint(6) = 1
  End If
  msti(6) = mint(6)
  Return
  1000 Format (1X, 'Error: illegal values in subroutine PYFRAM.', 1X, 'No transformation performed.'/1X, 'IFRAME =', 1X, I5, '; MINT(6) =', 1X, I5)
End Subroutine pyfram
