Subroutine luerrm(merr, chmess)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Character chmess*(*)
  Write (6, *) 'merr,chmess=', merr, chmess
  If (merr<=10) Then
    mstu(27) = mstu(27) + 1
    mstu(28) = merr
    If (mstu(25)==1 .And. mstu(27)<=mstu(26)) Write (mstu(11), 1000) merr, mstu(31), chmess
  Else If (merr<=20) Then
    mstu(23) = mstu(23) + 1
    mstu(24) = merr - 10
    If (mstu(21)>=1 .And. mstu(23)<=mstu(22)) Write (mstu(11), 1100) merr - 10, mstu(31), chmess
    If (mstu(21)>=2 .And. mstu(23)>mstu(22)) Then
      Write (mstu(11), 1100) merr - 10, mstu(31), chmess
      Write (mstu(11), 1200)
      If (merr/=17) Call lulist(2)
      Stop
    End If
  Else
    Write (mstu(11), 1300) merr - 20, mstu(31), chmess
    Stop
  End If
  Return
  1000 Format (/5X, 'Advisory warning type', I2, ' given after', I6, ' LUEXEC calls:'/5X, A)
  1100 Format (/5X, 'Error type', I2, ' has occured after', I6, ' LUEXEC calls:'/5X, A)
  1200 Format (5X, 'Execution will be stopped after listing of last ', 'event!')
  1300 Format (/5X, 'Fatal error type', I2, ' has occured after', I6, ' LUEXEC calls:'/5X, A/5X, 'Execution will now be stopped!')
End Subroutine luerrm
