Subroutine lugive(chin)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Character chin*(*), chfix*104, chbit*104, chold*8, chnew*8, chnam*4, chvar(17)*4, chalp(2)*26, chind*8, chini*10, chinr*16
  Data chvar/'N', 'K', 'P', 'V', 'MSTU', 'PARU', 'MSTJ', 'PARJ', 'KCHG', 'PMAS', 'PARF', 'VCKM', 'MDCY', 'MDME', 'BRAT', 'KFDP', 'CHAF'/
  Data chalp/'abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
  If (mstu(12)>=1) Call lulist(0)
  chbit = chin // ' '
  lbit = 101
  100 lbit = lbit - 1
  If (chbit(lbit:lbit)==' ') Goto 100
  ltot = 0
  Do lcom = 1, lbit
    If (chbit(lcom:lcom)==' ') Goto 110
    ltot = ltot + 1
    chfix(ltot:ltot) = chbit(lcom:lcom)
  110 End Do
  llow = 0
  120 lhig = llow + 1
  130 lhig = lhig + 1
  If (lhig<=ltot .And. chfix(lhig:lhig)/=';') Goto 130
  lbit = lhig - llow - 1
  chbit(1:lbit) = chfix(llow+1:lhig-1)
  lnam = 1
  140 lnam = lnam + 1
  If (chbit(lnam:lnam)/='(' .And. chbit(lnam:lnam)/='=' .And. lnam<=4) Goto 140
  chnam = chbit(1:lnam-1) // ' '
  Do lcom = 1, lnam - 1
    Do lalp = 1, 26
      If (chnam(lcom:lcom)==chalp(1)(lalp:lalp)) chnam(lcom:lcom) = chalp(2)(lalp:lalp)
    End Do
  End Do
  ivar = 0
  Do iv = 1, 17
    If (chnam==chvar(iv)) ivar = iv
  End Do
  If (ivar==0) Then
    Call luerrm(18, '(LUGIVE:) do not recognize variable '//chnam)
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If
  i = 0
  j = 0
  If (chbit(lnam:lnam)=='(') Then
    lind = lnam
    170 lind = lind + 1
    If (chbit(lind:lind)/=')' .And. chbit(lind:lind)/=',') Goto 170
    chind = ' '
    If ((chbit(lnam+1:lnam+1)=='C' .Or. chbit(lnam+1:lnam+1)=='c') .And. (ivar==9 .Or. ivar==10 .Or. ivar==13 .Or. ivar==17)) Then
      chind(lnam-lind+11:8) = chbit(lnam+2:lind-1)
      Read (chind, '(I8)') i1
      i = lucomp(i1)
    Else
      chind(lnam-lind+10:8) = chbit(lnam+1:lind-1)
      Read (chind, '(I8)') i
    End If
    lnam = lind
    If (chbit(lnam:lnam)==')') lnam = lnam + 1
  End If
  If (chbit(lnam:lnam)==',') Then
    lind = lnam
    180 lind = lind + 1
    If (chbit(lind:lind)/=')' .And. chbit(lind:lind)/=',') Goto 180
    chind = ' '
    chind(lnam-lind+10:8) = chbit(lnam+1:lind-1)
    Read (chind, '(I8)') j
    lnam = lind + 1
  End If
  ierr = 1
  If (chbit(lnam:lnam)/='=') Goto 190
  If (ivar==1) Then
    If (i/=0 .Or. j/=0) Goto 190
    iold = n
  Else If (ivar==2) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    iold = k(i, j)
  Else If (ivar==3) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    rold = p(i, j)
  Else If (ivar==4) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    rold = v(i, j)
  Else If (ivar==5) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    iold = mstu(i)
  Else If (ivar==6) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    rold = paru(i)
  Else If (ivar==7) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    iold = mstj(i)
  Else If (ivar==8) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    rold = parj(i)
  Else If (ivar==9) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>3) Goto 190
    iold = kchg(i, j)
  Else If (ivar==10) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>4) Goto 190
    rold = pmas(i, j)
  Else If (ivar==11) Then
    If (i<1 .Or. i>2000 .Or. j/=0) Goto 190
    rold = parf(i)
  Else If (ivar==12) Then
    If (i<1 .Or. i>4 .Or. j<1 .Or. j>4) Goto 190
    rold = vckm(i, j)
  Else If (ivar==13) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>3) Goto 190
    iold = mdcy(i, j)
  Else If (ivar==14) Then
    If (i<1 .Or. i>mstu(7) .Or. j<1 .Or. j>2) Goto 190
    iold = mdme(i, j)
  Else If (ivar==15) Then
    If (i<1 .Or. i>mstu(7) .Or. j/=0) Goto 190
    rold = brat(i)
  Else If (ivar==16) Then
    If (i<1 .Or. i>mstu(7) .Or. j<1 .Or. j>5) Goto 190
    iold = kfdp(i, j)
  Else If (ivar==17) Then
    If (i<1 .Or. i>mstu(6) .Or. j/=0) Goto 190
    chold = chaf(i)
  End If
  ierr = 0
  190 If (ierr==1) Then
    Call luerrm(18, '(LUGIVE:) unallowed indices for '//chbit(1:lnam-1))
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If
  If (lnam>=lbit) Then
    chbit(lnam:14) = ' '
    chbit(15:60) = ' has the value                                '
    If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
      Write (chbit(51:60), '(I10)') iold
    Else If (ivar/=17) Then
      Write (chbit(47:60), '(F14.5)') rold
    Else
      chbit(53:60) = chold
    End If
    If (mstu(13)>=1) Write (mstu(11), 1000) chbit(1:60)
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If
  If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
    chini = ' '
    chini(lnam-lbit+11:10) = chbit(lnam+1:lbit)
    Read (chini, '(I10)') inew
  Else If (ivar/=17) Then
    chinr = ' '
    chinr(lnam-lbit+17:16) = chbit(lnam+1:lbit)
    Read (chinr, '(F16.2)') rnew
  Else
    chnew = chbit(lnam+1:lbit) // ' '
  End If
  If (ivar==1) Then
    n = inew
  Else If (ivar==2) Then
    k(i, j) = inew
  Else If (ivar==3) Then
    p(i, j) = rnew
  Else If (ivar==4) Then
    v(i, j) = rnew
  Else If (ivar==5) Then
    mstu(i) = inew
  Else If (ivar==6) Then
    paru(i) = rnew
  Else If (ivar==7) Then
    mstj(i) = inew
  Else If (ivar==8) Then
    parj(i) = rnew
  Else If (ivar==9) Then
    kchg(i, j) = inew
  Else If (ivar==10) Then
    pmas(i, j) = rnew
  Else If (ivar==11) Then
    parf(i) = rnew
  Else If (ivar==12) Then
    vckm(i, j) = rnew
  Else If (ivar==13) Then
    mdcy(i, j) = inew
  Else If (ivar==14) Then
    mdme(i, j) = inew
  Else If (ivar==15) Then
    brat(i) = rnew
  Else If (ivar==16) Then
    kfdp(i, j) = inew
  Else If (ivar==17) Then
    chaf(i) = chnew
  End If
  chbit(lnam:14) = ' '
  chbit(15:60) = ' changed from                to               '
  If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
    Write (chbit(33:42), '(I10)') iold
    Write (chbit(51:60), '(I10)') inew
  Else If (ivar/=17) Then
    Write (chbit(29:42), '(F14.5)') rold
    Write (chbit(47:60), '(F14.5)') rnew
  Else
    chbit(35:42) = chold
    chbit(53:60) = chnew
  End If
  If (mstu(13)>=1) Write (mstu(11), 1000) chbit(1:60)
  llow = lhig
  If (llow<ltot) Goto 120
  Return
  1000 Format (5X, A60)
End Subroutine lugive
