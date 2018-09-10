Subroutine lulist(mlist)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Character chap*16, chac*16, chan*16, chad(5)*16, chmo(12)*3, chdl(7)*4
  Dimension ps(6)
  Data chmo/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/, chdl/'(())', ' ', '()', '!!', '<>', '==', '(==)'/
  If (mlist>=1 .And. mlist<=3) Then
    If (mlist==1) Write (mstu(11), 1100)
    If (mlist==2) Write (mstu(11), 1200)
    If (mlist==3) Write (mstu(11), 1300)
    lmx = 12
    If (mlist>=2) lmx = 16
    istr = 0
    imax = n
    If (mstu(2)>0) imax = mstu(2)
    Do i = max(1, mstu(1)), max(imax, n+max(0,mstu(3)))
      If ((i>imax .And. i<=n) .Or. k(i,1)<0) Goto 120
      Call luname(k(i,2), chap)
      len = 0
      Do lem = 1, 16
        If (chap(lem:lem)/=' ') len = lem
      End Do
      mdl = (k(i,1)+19)/10
      ldl = 0
      If (mdl==2 .Or. mdl>=8) Then
        chac = chap
        If (len>lmx) chac(lmx:lmx) = '?'
      Else
        ldl = 1
        If (mdl==1 .Or. mdl==7) ldl = 2
        If (len==0) Then
          chac = chdl(mdl)(1:2*ldl) // ' '
        Else
          chac = chdl(mdl)(1:ldl) // chap(1:min(len,lmx-2*ldl)) // chdl(mdl)(ldl+1:2*ldl) // ' '
          If (len+2*ldl>lmx) chac(lmx:lmx) = '?'
        End If
      End If
      If (k(i,1)==1 .Or. k(i,1)==2 .Or. k(i,1)==11 .Or. k(i,1)==12) Then
        kc = lucomp(k(i,2))
        kcc = 0
        If (kc/=0) kcc = kchg(kc, 2)
        If (kcc/=0 .And. istr==0) Then
          istr = 1
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'A'
        Else If (kcc/=0 .And. (k(i,1)==2 .Or. k(i,1)==12)) Then
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'I'
        Else If (kcc/=0) Then
          istr = 0
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'V'
        End If
      End If
      If (mlist==1 .And. abs(p(i,4))<9999.) Then
        Write (mstu(11), 1400) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mlist==1 .And. abs(p(i,4))<99999.) Then
        Write (mstu(11), 1500) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mlist==1) Then
        Write (mstu(11), 1600) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mstu(5)==10000 .And. (k(i,1)==3 .Or. k(i,1)==13 .Or. k(i,1)==14)) Then
        Write (mstu(11), 1700) i, chac, (k(i,j1), j1=1, 3), k(i, 4)/100000000, mod(k(i,4)/10000, 10000), mod(k(i,4), 10000), k(i, 5)/100000000, mod(k(i,5)/10000, 10000), mod(k(i,5), 10000), (p(i,j2), j2=1, 5)
      Else
        Write (mstu(11), 1800) i, chac, (k(i,j1), j1=1, 5), (p(i,j2), j2=1, 5)
      End If
      If (mlist==3) Write (mstu(11), 1900)(v(i,j), j=1, 5)
      If (mstu(70)>=1) Then
        isep = 0
        Do j = 1, min(10, mstu(70))
          If (i==mstu(70+j)) isep = 1
        End Do
        If (isep==1 .And. mlist==1) Write (mstu(11), 2000)
        If (isep==1 .And. mlist>=2) Write (mstu(11), 2100)
      End If
    120 End Do
    Do j = 1, 6
      ps(j) = plu(0, j)
    End Do
    If (mlist==1 .And. abs(ps(4))<9999.) Then
      Write (mstu(11), 2200) ps(6), (ps(j), j=1, 5)
    Else If (mlist==1 .And. abs(ps(4))<99999.) Then
      Write (mstu(11), 2300) ps(6), (ps(j), j=1, 5)
    Else If (mlist==1) Then
      Write (mstu(11), 2400) ps(6), (ps(j), j=1, 5)
    Else
      Write (mstu(11), 2500) ps(6), (ps(j), j=1, 5)
    End If
  Else If (mlist==11) Then
    Write (mstu(11), 2600)
    Do kf = 1, 40
      Call luname(kf, chap)
      Call luname(-kf, chan)
      If (chap/=' ' .And. chan==' ') Write (mstu(11), 2700) kf, chap
      If (chan/=' ') Write (mstu(11), 2700) kf, chap, -kf, chan
    End Do
    Do kfls = 1, 3, 2
      Do kfla = 1, 8
        Do kflb = 1, kfla - (3-kfls)/2
          kf = 1000*kfla + 100*kflb + kfls
          Call luname(kf, chap)
          Call luname(-kf, chan)
          Write (mstu(11), 2700) kf, chap, -kf, chan
        End Do
      End Do
    End Do
    Do kmul = 0, 5
      kfls = 3
      If (kmul==0 .Or. kmul==3) kfls = 1
      If (kmul==5) kfls = 5
      kflr = 0
      If (kmul==2 .Or. kmul==3) kflr = 1
      If (kmul==4) kflr = 2
      Do kflb = 1, 8
        Do kflc = 1, kflb - 1
          kf = 10000*kflr + 100*kflb + 10*kflc + kfls
          Call luname(kf, chap)
          Call luname(-kf, chan)
          Write (mstu(11), 2700) kf, chap, -kf, chan
        End Do
        kf = 10000*kflr + 110*kflb + kfls
        Call luname(kf, chap)
        Write (mstu(11), 2700) kf, chap
      End Do
    End Do
    kf = 130
    Call luname(kf, chap)
    Write (mstu(11), 2700) kf, chap
    kf = 310
    Call luname(kf, chap)
    Write (mstu(11), 2700) kf, chap
    Do kflsp = 1, 3
      kfls = 2 + 2*(kflsp/3)
      Do kfla = 1, 8
        Do kflb = 1, kfla
          Do kflc = 1, kflb
            If (kflsp==1 .And. (kfla==kflb .Or. kflb==kflc)) Goto 180
            If (kflsp==2 .And. kfla==kflc) Goto 180
            If (kflsp==1) kf = 1000*kfla + 100*kflc + 10*kflb + kfls
            If (kflsp>=2) kf = 1000*kfla + 100*kflb + 10*kflc + kfls
            Call luname(kf, chap)
            Call luname(-kf, chan)
            Write (mstu(11), 2700) kf, chap, -kf, chan
          180 End Do
        End Do
      End Do
    End Do
  Else If (mlist==12) Then
    Write (mstu(11), 2800)
    mstj24 = mstj(24)
    mstj(24) = 0
    kfmax = 20883
    If (mstu(2)/=0) kfmax = mstu(2)
    Do kf = max(1, mstu(1)), kfmax
      kc = lucomp(kf)
      If (kc==0) Goto 220
      If (mstu(14)==0 .And. kf>100 .And. kc<=100) Goto 220
      If (mstu(14)>0 .And. kf>100 .And. max(mod(kf/1000,10),mod(kf/100,10))>mstu(14)) Goto 220
      Call luname(kf, chap)
      If (kf<=100 .And. chap==' ' .And. mdcy(kc,2)==0) Goto 220
      Call luname(-kf, chan)
      pm = ulmass(kf)
      Write (mstu(11), 2900) kf, kc, chap, chan, kchg(kc, 1), kchg(kc, 2), kchg(kc, 3), pm, pmas(kc, 2), pmas(kc, 3), pmas(kc, 4), mdcy(kc, 1)
      If (kf>100 .And. kc<=100) Goto 220
      Do idc = mdcy(kc, 2), mdcy(kc, 2) + mdcy(kc, 3) - 1
        Do j = 1, 5
          Call luname(kfdp(idc,j), chad(j))
        End Do
        Write (mstu(11), 3000) idc, mdme(idc, 1), mdme(idc, 2), brat(idc), (chad(j), j=1, 5)
      End Do
    220 End Do
    mstj(24) = mstj24
  Else If (mlist==13) Then
    Write (mstu(11), 3100)
    Do i = 1, 200
      Write (mstu(11), 3200) i, mstu(i), paru(i), mstj(i), parj(i), parf(i)
    End Do
  End If
  Return
  1100 Format (///28X, 'Event listing (summary)'//4X, 'I  particle/jet KS', 5X, 'KF orig    p_x      p_y      p_z       E        m'/)
  1200 Format (///28X, 'Event listing (standard)'//4X, 'I  particle/jet', '  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', '       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/)
  1300 Format (///28X, 'Event listing (with vertices)'//4X, 'I  particle/j', 'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', '       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X, 'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/)
  1400 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.3)
  1500 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.2)
  1600 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.1)
  1700 Format (1X, I4, 2X, A16, 1X, I3, 1X, I8, 2X, I4, 2(3X,I1,2I4), 5F13.5)
  1800 Format (1X, I4, 2X, A16, 1X, I3, 1X, I8, 2X, I4, 2(3X,I9), 5F13.5)
  1900 Format (66X, 5(1X,F12.3))
  2000 Format (1X, 78('='))
  2100 Format (1X, 130('='))
  2200 Format (19X, 'sum:', F6.2, 5X, 5F9.3)
  2300 Format (19X, 'sum:', F6.2, 5X, 5F9.2)
  2400 Format (19X, 'sum:', F6.2, 5X, 5F9.1)
  2500 Format (19X, 'sum charge:', F6.2, 3X, 'sum momentum and inv. mass:', 5F13.5)
  2600 Format (///20X, 'List of KF codes in program'/)
  2700 Format (4X, I6, 4X, A16, 6X, I6, 4X, A16)
  2800 Format (///30X, 'Particle/parton data table'//5X, 'KF', 5X, 'KC', 4X, 'particle', 8X, 'antiparticle', 6X, 'chg  col  anti', 8X, 'mass', 7X, 'width', 7X, 'w-cut', 5X, 'lifetime', 1X, 'decay'/11X, 'IDC', 1X, 'on/off', 1X, 'ME', 3X, 'Br.rat.', 4X, 'decay products')
  2900 Format (/1X, I6, 3X, I4, 4X, A16, A16, 3I5, 1X, F12.5, 2(1X,F11.5), 2X, F12.5, 3X, I2)
  3000 Format (10X, I4, 2X, I3, 2X, I3, 2X, F8.5, 4X, 5A16)
  3100 Format (///20X, 'Parameter value table'//4X, 'I', 3X, 'MSTU(I)', 8X, 'PARU(I)', 3X, 'MSTJ(I)', 8X, 'PARJ(I)', 8X, 'PARF(I)')
  3200 Format (1X, I4, 1X, I9, 1X, F14.5, 1X, I9, 1X, F14.5, 1X, F14.5)
End Subroutine lulist
