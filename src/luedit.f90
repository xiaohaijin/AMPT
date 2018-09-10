Subroutine luedit(medit)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension ns(2), pts(2), pls(2)
  If ((medit>=0 .And. medit<=3) .Or. medit==5) Then
    imax = n
    If (mstu(2)>0) imax = mstu(2)
    i1 = max(1, mstu(1)) - 1
    Do i = max(1, mstu(1)), imax
      If (k(i,1)==0 .Or. k(i,1)>20) Goto 110
      If (medit==1) Then
        If (k(i,1)>10) Goto 110
      Else If (medit==2) Then
        If (k(i,1)>10) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 110
      Else If (medit==3) Then
        If (k(i,1)>10) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0) Goto 110
        If (kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 110
      Else If (medit==5) Then
        If (k(i,1)==13 .Or. k(i,1)==14) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0) Goto 110
        If (k(i,1)>=11 .And. kchg(kc,2)==0) Goto 110
      End If
      i1 = i1 + 1
      Do j = 1, 5
        k(i1, j) = k(i, j)
        p(i1, j) = p(i, j)
        v(i1, j) = v(i, j)
      End Do
      k(i1, 3) = 0
    110 End Do
    n = i1
  Else If (medit>=11 .And. medit<=15) Then
    i1 = 0
    Do i = 1, n
      k(i, 3) = mod(k(i,3), mstu(5))
      If (medit==11 .And. k(i,1)<0) Goto 120
      If (medit==12 .And. k(i,1)==0) Goto 120
      If (medit==13 .And. (k(i,1)==11 .Or. k(i,1)==12 .Or. k(i,1)==15) .And. k(i,2)/=94) Goto 120
      If (medit==14 .And. (k(i,1)==13 .Or. k(i,1)==14 .Or. k(i,2)==94)) Goto 120
      If (medit==15 .And. k(i,1)>=21) Goto 120
      i1 = i1 + 1
      k(i, 3) = k(i, 3) + mstu(5)*i1
    120 End Do
    Do i = 1, n
      If (k(i,1)<=0 .Or. k(i,1)>20 .Or. k(i,3)/mstu(5)==0) Goto 140
      id = i
      130 im = mod(k(id,3), mstu(5))
      If (medit==13 .And. im>0 .And. im<=n) Then
        If ((k(im,1)==11 .Or. k(im,1)==12 .Or. k(im,1)==15) .And. k(im,2)/=94) Then
          id = im
          Goto 130
        End If
      Else If (medit==14 .And. im>0 .And. im<=n) Then
        If (k(im,1)==13 .Or. k(im,1)==14 .Or. k(im,2)==94) Then
          id = im
          Goto 130
        End If
      End If
      k(i, 3) = mstu(5)*(k(i,3)/mstu(5))
      If (im/=0) k(i, 3) = k(i, 3) + k(im, 3)/mstu(5)
      If (k(i,1)/=3 .And. k(i,1)/=13 .And. k(i,1)/=14) Then
        If (k(i,4)>0 .And. k(i,4)<=mstu(4)) k(i, 4) = k(k(i,4), 3)/mstu(5)
        If (k(i,5)>0 .And. k(i,5)<=mstu(4)) k(i, 5) = k(k(i,5), 3)/mstu(5)
      Else
        kcm = mod(k(i,4)/mstu(5), mstu(5))
        If (kcm>0 .And. kcm<=mstu(4)) kcm = k(kcm, 3)/mstu(5)
        kcd = mod(k(i,4), mstu(5))
        If (kcd>0 .And. kcd<=mstu(4)) kcd = k(kcd, 3)/mstu(5)
        k(i, 4) = mstu(5)**2*(k(i,4)/mstu(5)**2) + mstu(5)*kcm + kcd
        kcm = mod(k(i,5)/mstu(5), mstu(5))
        If (kcm>0 .And. kcm<=mstu(4)) kcm = k(kcm, 3)/mstu(5)
        kcd = mod(k(i,5), mstu(5))
        If (kcd>0 .And. kcd<=mstu(4)) kcd = k(kcd, 3)/mstu(5)
        k(i, 5) = mstu(5)**2*(k(i,5)/mstu(5)**2) + mstu(5)*kcm + kcd
      End If
    140 End Do
    i1 = 0
    Do i = 1, n
      If (k(i,3)/mstu(5)==0) Goto 160
      i1 = i1 + 1
      Do j = 1, 5
        k(i1, j) = k(i, j)
        p(i1, j) = p(i, j)
        v(i1, j) = v(i, j)
      End Do
      k(i1, 3) = mod(k(i1,3), mstu(5))
    160 End Do
    n = i1
  Else If (medit==21) Then
    If (2*n>=mstu(4)) Then
      Call luerrm(11, '(LUEDIT:) no more memory left in LUJETS')
      Return
    End If
    Do i = 1, n
      Do j = 1, 5
        k(mstu(4)-i, j) = k(i, j)
        p(mstu(4)-i, j) = p(i, j)
        v(mstu(4)-i, j) = v(i, j)
      End Do
    End Do
    mstu(32) = n
  Else If (medit==22) Then
    Do i = 1, mstu(32)
      Do j = 1, 5
        k(i, j) = k(mstu(4)-i, j)
        p(i, j) = p(mstu(4)-i, j)
        v(i, j) = v(mstu(4)-i, j)
      End Do
    End Do
    n = mstu(32)
  Else If (medit==23) Then
    i1 = 0
    Do i = 1, n
      kh = k(i, 3)
      If (kh>=1) Then
        If (k(kh,1)>20) kh = 0
      End If
      If (kh/=0) Goto 200
      i1 = i1 + 1
      If (k(i,1)>10 .And. k(i,1)<=20) k(i, 1) = k(i, 1) - 10
    End Do
    200 n = i1
  Else If (medit==31 .Or. medit==32) Then
    Call ludbrb(1, n+mstu(3), 0., -ulangl(p(mstu(61),1),p(mstu(61),2)), 0D0, 0D0, 0D0)
    Call ludbrb(1, n+mstu(3), -ulangl(p(mstu(61),3),p(mstu(61),1)), 0., 0D0, 0D0, 0D0)
    Call ludbrb(1, n+mstu(3), 0., -ulangl(p(mstu(61)+1,1),p(mstu(61)+1,2)), 0D0, 0D0, 0D0)
    If (medit==31) Return
    Do is = 1, 2
      ns(is) = 0
      pts(is) = 0.
      pls(is) = 0.
    End Do
    Do i = 1, n
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 220
      If (mstu(41)>=2) Then
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 220
        If (mstu(41)>=3 .And. kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 220
      End If
      is = int(2.-sign(0.5,p(i,3)))
      ns(is) = ns(is) + 1
      pts(is) = pts(is) + sqrt(p(i,1)**2+p(i,2)**2)
    220 End Do
    If (ns(1)*pts(2)**2<ns(2)*pts(1)**2) Call ludbrb(1, n+mstu(3), paru(1), 0., 0D0, 0D0, 0D0)
    Do i = 1, n
      If (p(i,3)>=0.) Goto 230
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 230
      If (mstu(41)>=2) Then
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 230
        If (mstu(41)>=3 .And. kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 230
      End If
      is = int(2.-sign(0.5,p(i,1)))
      pls(is) = pls(is) - p(i, 3)
    230 End Do
    If (pls(2)>pls(1)) Call ludbrb(1, n+mstu(3), 0., paru(1), 0D0, 0D0, 0D0)
  End If
  Return
End Subroutine luedit
