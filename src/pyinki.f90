Subroutine pyinki(chfram, chbeam, chtarg, win)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Character chfram*8, chbeam*8, chtarg*8, chcom(3)*8, chalp(2)*26, chidnt(3)*8, chtemp*8, chcde(18)*8, chinit*76
  Dimension len(3), kcde(18)
  Data chalp/'abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
  Data chcde/'e-      ', 'e+      ', 'nue     ', 'nue~    ', 'mu-     ', 'mu+     ', 'numu    ', 'numu~   ', 'tau-    ', 'tau+    ', 'nutau   ', 'nutau~  ', 'pi+     ', 'pi-     ', 'n       ', 'n~      ', 'p       ', 'p~      '/
  Data kcde/11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16, 211, -211, 2112, -2112, 2212, -2212/
  chcom(1) = chfram
  chcom(2) = chbeam
  chcom(3) = chtarg
  Do i = 1, 3
    len(i) = 8
    Do ll = 8, 1, -1
      If (len(i)==ll .And. chcom(i)(ll:ll)==' ') len(i) = ll - 1
      Do la = 1, 26
        If (chcom(i)(ll:ll)==chalp(2)(la:la)) chcom(i)(ll:ll) = chalp(1)(la:la)
      End Do
    End Do
    chidnt(i) = chcom(i)
    Do ll = 1, 6
      If (chidnt(i)(ll:ll+2)=='bar') Then
        chtemp = chidnt(i)
        chidnt(i) = chtemp(1:ll-1) // '~' // chtemp(ll+3:8) // '  '
      End If
    End Do
    Do ll = 1, 8
      If (chidnt(i)(ll:ll)=='_') Then
        chtemp = chidnt(i)
        chidnt(i) = chtemp(1:ll-1) // chtemp(ll+1:8) // ' '
      End If
    End Do
  End Do
  n = 2
  Do i = 1, 2
    k(i, 2) = 0
    Do j = 1, 18
      If (chidnt(i+1)==chcde(j)) k(i, 2) = kcde(j)
    End Do
    p(i, 5) = ulmass(k(i,2))
    mint(40+i) = 1
    If (iabs(k(i,2))>100) mint(40+i) = 2
    Do j = 1, 5
      v(i, j) = 0.
    End Do
  End Do
  If (k(1,2)==0) Write (mstu(11), 1000) chbeam(1:len(2))
  If (k(2,2)==0) Write (mstu(11), 1100) chtarg(1:len(3))
  If (k(1,2)==0 .Or. k(2,2)==0) Stop
  Do j = 6, 10
    vint(j) = 0.
  End Do
  chinit = ' '
  If (chcom(1)(1:2)=='cm') Then
    If (chcom(2)(1:1)/='e') Then
      loffs = (34-(len(2)+len(3)))/2
      chinit(loffs+1:76) = 'PYTHIA will be initialized for a ' // chcom(2)(1:len(2)) // '-' // chcom(3)(1:len(3)) // ' collider' // ' '
    Else
      loffs = (33-(len(2)+len(3)))/2
      chinit(loffs+1:76) = 'PYTHIA will be initialized for an ' // chcom(2)(1:len(2)) // '-' // chcom(3)(1:len(3)) // ' collider' // ' '
    End If
    s = win**2
    p(1, 1) = 0.
    p(1, 2) = 0.
    p(2, 1) = 0.
    p(2, 2) = 0.
    p(1, 3) = sqrt(((s-p(1,5)**2-p(2,5)**2)**2-(2.*p(1,5)*p(2,5))**2)/(4.*s))
    p(2, 3) = -p(1, 3)
    p(1, 4) = sqrt(p(1,3)**2+p(1,5)**2)
    p(2, 4) = sqrt(p(2,3)**2+p(2,5)**2)
  Else If (chcom(1)(1:3)=='fix') Then
    loffs = (29-(len(2)+len(3)))/2
    chinit(loffs+1:76) = 'PYTHIA will be initialized for ' // chcom(2)(1:len(2)) // ' on ' // chcom(3)(1:len(3)) // ' fixed target' // ' '
    p(1, 1) = 0.
    p(1, 2) = 0.
    p(2, 1) = 0.
    p(2, 2) = 0.
    p(1, 3) = win
    p(1, 4) = sqrt(p(1,3)**2+p(1,5)**2)
    p(2, 3) = 0.
    p(2, 4) = p(2, 5)
    s = p(1, 5)**2 + p(2, 5)**2 + 2.*p(2, 4)*p(1, 4)
    vint(10) = p(1, 3)/(p(1,4)+p(2,4))
    Call lurobo(0., 0., 0., 0., -vint(10))
  Else If (chcom(1)(1:3)=='use') Then
    loffs = (13-(len(1)+len(2)))/2
    chinit(loffs+1:76) = 'PYTHIA will be initialized for ' // chcom(2)(1:len(2)) // ' on ' // chcom(3)(1:len(3)) // 'user-specified configuration' // ' '
    p(1, 4) = sqrt(p(1,1)**2+p(1,2)**2+p(1,3)**2+p(1,5)**2)
    p(2, 4) = sqrt(p(2,1)**2+p(2,2)**2+p(2,3)**2+p(2,5)**2)
    Do j = 1, 3
      vint(7+j) = sngl((dble(p(1,j))+dble(p(2,j)))/dble(p(1,4)+p(2,4)))
    End Do
    Call lurobo(0., 0., -vint(8), -vint(9), -vint(10))
    vint(7) = ulangl(p(1,1), p(1,2))
    Call lurobo(0., -vint(7), 0., 0., 0.)
    vint(6) = ulangl(p(1,3), p(1,1))
    Call lurobo(-vint(6), 0., 0., 0., 0.)
    s = p(1, 5)**2 + p(2, 5)**2 + 2.*(p(1,4)*p(2,4)-p(1,3)*p(2,3))
  Else
    Write (mstu(11), 1800) chfram(1:len(1))
    Stop
  End If
  If (s<parp(2)**2) Then
    Write (mstu(11), 1900) sqrt(s)
    Stop
  End If
  mint(11) = k(1, 2)
  mint(12) = k(2, 2)
  mint(43) = 2*mint(41) + mint(42) - 2
  vint(1) = sqrt(s)
  vint(2) = s
  vint(3) = p(1, 5)
  vint(4) = p(2, 5)
  vint(5) = p(1, 3)
  If (mstp(82)<=1) vint(149) = 4.*parp(81)**2/s
  If (mstp(82)>=2) vint(149) = 4.*parp(82)**2/s
  Return
  1000 Format (1X, 'Error: unrecognized beam particle ''', A, '''.'/1X, 'Execution stopped!')
  1100 Format (1X, 'Error: unrecognized target particle ''', A, '''.'/1X, 'Execution stopped!')
  1800 Format (1X, 'Error: unrecognized coordinate frame ''', A, '''.'/1X, 'Execution stopped!')
  1900 Format (1X, 'Error: too low CM energy,', F8.3, ' GeV for event ', 'generation.'/1X, 'Execution stopped!')
End Subroutine pyinki
