Subroutine lustrf(ip)
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension dps(5), kfl(3), pmq(3), px(3), py(3), gam(3), ie(2), pr(2), in(9), dhm(4), dhg(4), dp(5, 5), irank(2), mju(4), iju(3), pju(5, 5), tju(5), kfjh(2), njs(2), kfjs(2), pjs(4, 5)
  four(i, j) = p(i, 4)*p(j, 4) - p(i, 1)*p(j, 1) - p(i, 2)*p(j, 2) - p(i, 3)*p(j, 3)
  dfour(i, j) = dp(i, 4)*dp(j, 4) - dp(i, 1)*dp(j, 1) - dp(i, 2)*dp(j, 2) - dp(i, 3)*dp(j, 3)
  mstj(91) = 0
  nsav = n
  np = 0
  kqsum = 0
  Do j = 1, 5
    dps(j) = 0D0
  End Do
  mju(1) = 0
  mju(2) = 0
  i = ip - 1
  110 i = i + 1
  If (i>min(n,mstu(4)-mstu(32))) Then
    Call luerrm(12, '(LUSTRF:) failed to reconstruct jet system')
    If (mstu(21)>=1) Return
  End If
  If (k(i,1)/=1 .And. k(i,1)/=2 .And. k(i,1)/=41) Goto 110
  kc = lucomp(k(i,2))
  If (kc==0) Goto 110
  kq = kchg(kc, 2)*isign(1, k(i,2))
  If (kq==0) Goto 110
  If (n+5*np+11>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  np = np + 1
  Do j = 1, 5
    k(n+np, j) = k(i, j)
    p(n+np, j) = p(i, j)
    dps(j) = dps(j) + dble(p(i,j))
  End Do
  k(n+np, 3) = i
  If (p(n+np,4)**2<p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2) Then
    p(n+np, 4) = sqrt(p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2+p(n+np,5)**2)
    dps(4) = dps(4) + dble(max(0.,p(n+np,4)-p(i,4)))
  End If
  If (kq/=2) kqsum = kqsum + kq
  If (k(i,1)==41) Then
    kqsum = kqsum + 2*kq
    If (kqsum==kq) mju(1) = n + np
    If (kqsum/=kq) mju(2) = n + np
  End If
  If (k(i,1)==2 .Or. k(i,1)==41) Goto 110
  If (kqsum/=0) Then
    Call luerrm(12, '(LUSTRF:) unphysical flavour combination')
    If (mstu(21)>=1) Return
  End If
  Call ludbrb(n+1, n+np, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))
  ntryr = 0
  paru12 = paru(12)
  paru13 = paru(13)
  mju(3) = mju(1)
  mju(4) = mju(2)
  nr = np
  130 If (nr>=3) Then
    pdrmin = 2.*paru12
    Do i = n + 1, n + nr
      If (i==n+nr .And. iabs(k(n+1,2))/=21) Goto 140
      i1 = i + 1
      If (i==n+nr) i1 = n + 1
      If (k(i,1)==41 .Or. k(i1,1)==41) Goto 140
      If (mju(1)/=0 .And. i1<mju(1) .And. iabs(k(i1,2))/=21) Goto 140
      If (mju(2)/=0 .And. i>mju(2) .And. iabs(k(i,2))/=21) Goto 140
      pap = sqrt((p(i,1)**2+p(i,2)**2+p(i,3)**2)*(p(i1,1)**2+p(i1,2)**2+p(i1,3)**2))
      pvp = p(i, 1)*p(i1, 1) + p(i, 2)*p(i1, 2) + p(i, 3)*p(i1, 3)
      pdr = 4.*(pap-pvp)**2/(paru13**2*pap+2.*(pap-pvp))
      If (pdr<pdrmin) Then
        ir = i
        pdrmin = pdr
      End If
    140 End Do
    If (pdrmin<paru12 .And. ir==n+nr) Then
      Do j = 1, 4
        p(n+1, j) = p(n+1, j) + p(n+nr, j)
      End Do
      p(n+1, 5) = sqrt(max(0.,p(n+1,4)**2-p(n+1,1)**2-p(n+1,2)**2-p(n+1,3)**2))
      nr = nr - 1
      Goto 130
    Else If (pdrmin<paru12) Then
      Do j = 1, 4
        p(ir, j) = p(ir, j) + p(ir+1, j)
      End Do
      p(ir, 5) = sqrt(max(0.,p(ir,4)**2-p(ir,1)**2-p(ir,2)**2-p(ir,3)**2))
      Do i = ir + 1, n + nr - 1
        k(i, 2) = k(i+1, 2)
        Do j = 1, 5
          p(i, j) = p(i+1, j)
        End Do
      End Do
      If (ir==n+nr-1) k(ir, 2) = k(n+nr, 2)
      nr = nr - 1
      If (mju(1)>ir) mju(1) = mju(1) - 1
      If (mju(2)>ir) mju(2) = mju(2) - 1
      Goto 130
    End If
  End If
  ntryr = ntryr + 1
  nrs = max(5*nr+11, np)
  ntry = 0
  180 ntry = ntry + 1
  If (ntry>100 .And. ntryr<=4) Then
    paru12 = 4.*paru12
    paru13 = 2.*paru13
    Goto 130
  Else If (ntry>100) Then
    Call luerrm(14, '(LUSTRF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = n + nrs
  If (mju(1)==0 .And. mju(2)==0) Goto 500
  Do jt = 1, 2
    njs(jt) = 0
    If (mju(jt)==0) Goto 490
    js = 3 - 2*jt
    Do iu = 1, 3
      iju(iu) = 0
      Do j = 1, 5
        pju(iu, j) = 0.
      End Do
    End Do
    iu = 0
    Do i1 = n + 1 + (jt-1)*(nr-1), n + nr + (jt-1)*(1-nr), js
      If (k(i1,2)/=21 .And. iu<=2) Then
        iu = iu + 1
        iju(iu) = i1
      End If
      Do j = 1, 4
        pju(iu, j) = pju(iu, j) + p(i1, j)
      End Do
    End Do
    Do iu = 1, 3
      pju(iu, 5) = sqrt(pju(iu,1)**2+pju(iu,2)**2+pju(iu,3)**2)
    End Do
    If (k(iju(3),2)/100/=10*k(iju(1),2)+k(iju(2),2) .And. k(iju(3),2)/100/=10*k(iju(2),2)+k(iju(1),2)) Then
      Call luerrm(12, '(LUSTRF:) unphysical flavour combination')
      If (mstu(21)>=1) Return
    End If
    t12 = (pju(1,1)*pju(2,1)+pju(1,2)*pju(2,2)+pju(1,3)*pju(2,3))/(pju(1,5)*pju(2,5))
    t13 = (pju(1,1)*pju(3,1)+pju(1,2)*pju(3,2)+pju(1,3)*pju(3,3))/(pju(1,5)*pju(3,5))
    t23 = (pju(2,1)*pju(3,1)+pju(2,2)*pju(3,2)+pju(2,3)*pju(3,3))/(pju(2,5)*pju(3,5))
    t11 = sqrt((2./3.)*(1.-t12)*(1.-t13)/(1.-t23))
    t22 = sqrt((2./3.)*(1.-t12)*(1.-t23)/(1.-t13))
    tsq = sqrt((2.*t11*t22+t12-1.)*(1.+t12))
    t1f = (tsq-t22*(1.+t12))/(1.-t12**2)
    t2f = (tsq-t11*(1.+t12))/(1.-t12**2)
    Do j = 1, 3
      tju(j) = -(t1f*pju(1,j)/pju(1,5)+t2f*pju(2,j)/pju(2,5))
    End Do
    tju(4) = sqrt(1.+tju(1)**2+tju(2)**2+tju(3)**2)
    Do iu = 1, 3
      pju(iu, 5) = tju(4)*pju(iu, 4) - tju(1)*pju(iu, 1) - tju(2)*pju(iu, 2) - tju(3)*pju(iu, 3)
    End Do
    If (pju(1,5)+pju(2,5)>pju(1,4)+pju(2,4)) Then
      Do j = 1, 3
        tju(j) = 0.
      End Do
      tju(4) = 1.
      pju(1, 5) = pju(1, 4)
      pju(2, 5) = pju(2, 4)
      pju(3, 5) = pju(3, 4)
    End If
    ista = i
    Do iu = 1, 2
      ns = iju(iu+1) - iju(iu)
      Do is = 1, ns
        is1 = iju(iu) + is - 1
        is2 = iju(iu) + is
        Do j = 1, 5
          dp(1, j) = dble(0.5*p(is1,j))
          If (is==1) dp(1, j) = dble(p(is1,j))
          dp(2, j) = dble(0.5*p(is2,j))
          If (is==ns) dp(2, j) = -dble(pju(iu,j))
        End Do
        If (is==ns) dp(2, 4) = dble(sqrt(pju(iu,1)**2+pju(iu,2)**2+pju(iu,3)**2))
        If (is==ns) dp(2, 5) = 0D0
        dp(3, 5) = dfour(1, 1)
        dp(4, 5) = dfour(2, 2)
        dhkc = dfour(1, 2)
        If (dp(3,5)+2D0*dhkc+dp(4,5)<=0D0) Then
          dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
          dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
          dp(3, 5) = 0D0
          dp(4, 5) = 0D0
          dhkc = dfour(1, 2)
        End If
        dhks = sqrt(dhkc**2-dp(3,5)*dp(4,5))
        dhk1 = 0.5D0*((dp(4,5)+dhkc)/dhks-1D0)
        dhk2 = 0.5D0*((dp(3,5)+dhkc)/dhks-1D0)
        in1 = n + nr + 4*is - 3
        p(in1, 5) = sngl(sqrt(dp(3,5)+2D0*dhkc+dp(4,5)))
        Do j = 1, 4
          p(in1, j) = sngl((1D0+dhk1)*dp(1,j)-dhk2*dp(2,j))
          p(in1+1, j) = sngl((1D0+dhk2)*dp(2,j)-dhk1*dp(1,j))
        End Do
      End Do
      isav = i
      270 ntry = ntry + 1
      If (ntry>100 .And. ntryr<=4) Then
        paru12 = 4.*paru12
        paru13 = 2.*paru13
        Goto 130
      Else If (ntry>100) Then
        Call luerrm(14, '(LUSTRF:) caught in infinite loop')
        If (mstu(21)>=1) Return
      End If
      i = isav
      irankj = 0
      ie(1) = k(n+1+(jt/2)*(np-1), 3)
      in(4) = n + nr + 1
      in(5) = in(4) + 1
      in(6) = n + nr + 4*ns + 1
      Do jq = 1, 2
        Do in1 = n + nr + 2 + jq, n + nr + 4*ns - 2 + jq, 4
          p(in1, 1) = 2 - jq
          p(in1, 2) = jq - 1
          p(in1, 3) = 1.
        End Do
      End Do
      kfl(1) = k(iju(iu), 2)
      px(1) = 0.
      py(1) = 0.
      gam(1) = 0.
      Do j = 1, 5
        pju(iu+3, j) = 0.
      End Do
      Do j = 1, 4
        dp(1, j) = dble(p(in(4),j))
        dp(2, j) = dble(p(in(4)+1,j))
        dp(3, j) = 0D0
        dp(4, j) = 0D0
      End Do
      dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
      dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
      dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
      dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
      dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
      If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1D0
      If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1D0
      If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1D0
      If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1D0
      dhc12 = dfour(1, 2)
      dhcx1 = dfour(3, 1)/dhc12
      dhcx2 = dfour(3, 2)/dhc12
      dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
      dhcy1 = dfour(4, 1)/dhc12
      dhcy2 = dfour(4, 2)/dhc12
      dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
      dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
      Do j = 1, 4
        dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
        p(in(6), j) = sngl(dp(3,j))
        p(in(6)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
      End Do
      320 i = i + 1
      If (2*i-nsav>=mstu(4)-mstu(32)-5) Then
        Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
        If (mstu(21)>=1) Return
      End If
      irankj = irankj + 1
      k(i, 1) = 1
      k(i, 3) = ie(1)
      k(i, 4) = 0
      k(i, 5) = 0
      330 Call lukfdi(kfl(1), 0, kfl(3), k(i,2))
      If (k(i,2)==0) Goto 270
      If (mstj(12)>=3 .And. irankj==1 .And. iabs(kfl(1))<=10 .And. iabs(kfl(3))>10) Then
        If (rlu(0)>parj(19)) Goto 330
      End If
      p(i, 5) = ulmass(k(i,2))
      Call luptdi(kfl(1), px(3), py(3))
      pr(1) = p(i, 5)**2 + (px(1)+px(3))**2 + (py(1)+py(3))**2
      Call luzdis(kfl(1), kfl(3), pr(1), z)
      gam(3) = (1.-z)*(gam(1)+pr(1)/z)
      Do j = 1, 3
        in(j) = in(3+j)
      End Do
      If (in(1)+1==in(2) .And. z*p(in(1)+2,3)*p(in(2)+2,3)*p(in(1),5)**2>=pr(1)) Then
        p(in(1)+2, 4) = z*p(in(1)+2, 3)
        p(in(2)+2, 4) = pr(1)/(p(in(1)+2,4)*p(in(1),5)**2)
        Do j = 1, 4
          p(i, j) = (px(1)+px(3))*p(in(3), j) + (py(1)+py(3))*p(in(3)+1, j)
        End Do
        Goto 420
      Else If (in(1)+1==in(2)) Then
        p(in(2)+2, 4) = p(in(2)+2, 3)
        p(in(2)+2, 1) = 1.
        in(2) = in(2) + 4
        If (in(2)>n+nr+4*ns) Goto 270
        If (four(in(1),in(2))<=1E-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
        End If
      End If
      360 If (in(1)>n+nr+4*ns .Or. in(2)>n+nr+4*ns .Or. in(1)>in(2)) Goto 270
      If (in(1)/=in(4) .Or. in(2)/=in(5)) Then
        Do j = 1, 4
          dp(1, j) = dble(p(in(1),j))
          dp(2, j) = dble(p(in(2),j))
          dp(3, j) = 0D0
          dp(4, j) = 0D0
        End Do
        dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
        dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
        dhc12 = dfour(1, 2)
        If (dhc12<=1D-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
          Goto 360
        End If
        in(3) = n + nr + 4*ns + 5
        dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
        dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
        dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
        If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1D0
        If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1D0
        If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1D0
        If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1D0
        dhcx1 = dfour(3, 1)/dhc12
        dhcx2 = dfour(3, 2)/dhc12
        dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
        dhcy1 = dfour(4, 1)/dhc12
        dhcy2 = dfour(4, 2)/dhc12
        dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
        dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
        Do j = 1, 4
          dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
          p(in(3), j) = sngl(dp(3,j))
          p(in(3)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
        End Do
        pxp = -(px(3)*four(in(6),in(3))+py(3)*four(in(6)+1,in(3)))
        pyp = -(px(3)*four(in(6),in(3)+1)+py(3)*four(in(6)+1,in(3)+1))
        If (abs(pxp**2+pyp**2-px(3)**2-py(3)**2)<0.01) Then
          px(3) = pxp
          py(3) = pyp
        End If
      End If
      Do j = 1, 4
        dhg(j) = 0D0
        p(i, j) = px(1)*p(in(6), j) + py(1)*p(in(6)+1, j) + px(3)*p(in(3), j) + py(3)*p(in(3)+1, j)
        Do in1 = in(4), in(1) - 4, 4
          p(i, j) = p(i, j) + p(in1+2, 3)*p(in1, j)
        End Do
        Do in2 = in(5), in(2) - 4, 4
          p(i, j) = p(i, j) + p(in2+2, 3)*p(in2, j)
        End Do
      End Do
      dhm(1) = dble(four(i,i))
      dhm(2) = dble(2.*four(i,in(1)))
      dhm(3) = dble(2.*four(i,in(2)))
      dhm(4) = dble(2.*four(in(1),in(2)))
      Do in2 = in(1) + 1, in(2), 4
        Do in1 = in(1), in2 - 1, 4
          dhc = dble(2.*four(in1,in2))
          dhg(1) = dhg(1) + dble(p(in1+2,1)*p(in2+2,1))*dhc
          If (in1==in(1)) dhg(2) = dhg(2) - dble(p(in2+2,1))*dhc
          If (in2==in(2)) dhg(3) = dhg(3) + dble(p(in1+2,1))*dhc
          If (in1==in(1) .And. in2==in(2)) dhg(4) = dhg(4) - dhc
        End Do
      End Do
      dhs1 = dhm(3)*dhg(4) - dhm(4)*dhg(3)
      If (dabs(dhs1)<1D-4) Goto 270
      dhs2 = dhm(4)*(dble(gam(3))-dhg(1)) - dhm(2)*dhg(3) - dhg(4)*(dble(p(i,5))**2-dhm(1)) + dhg(2)*dhm(3)
      dhs3 = dhm(2)*(dble(gam(3))-dhg(1)) - dhg(2)*(dble(p(i,5))**2-dhm(1))
      p(in(2)+2, 4) = 0.5*sngl(sqrt(max(0D0,dhs2**2-4D0*dhs1*dhs3))/abs(dhs1)-dhs2/dhs1)
      If (dhm(2)+dhm(4)*dble(p(in(2)+2,4))<=0D0) Goto 270
      p(in(1)+2, 4) = (p(i,5)**2-sngl(dhm(1))-sngl(dhm(3))*p(in(2)+2,4))/(sngl(dhm(2))+sngl(dhm(4))*p(in(2)+2,4))
      If (p(in(2)+2,4)>p(in(2)+2,3)) Then
        p(in(2)+2, 4) = p(in(2)+2, 3)
        p(in(2)+2, 1) = 1.
        in(2) = in(2) + 4
        If (in(2)>n+nr+4*ns) Goto 270
        If (four(in(1),in(2))<=1E-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
        End If
        Goto 360
      Else If (p(in(1)+2,4)>p(in(1)+2,3)) Then
        p(in(1)+2, 4) = p(in(1)+2, 3)
        p(in(1)+2, 1) = 0.
        in(1) = in(1) + js
        Goto 710
      End If
      420 Do j = 1, 4
        p(i, j) = p(i, j) + p(in(1)+2, 4)*p(in(1), j) + p(in(2)+2, 4)*p(in(2), j)
        pju(iu+3, j) = pju(iu+3, j) + p(i, j)
      End Do
      If (p(i,4)<=0.) Goto 270
      pju(iu+3, 5) = tju(4)*pju(iu+3, 4) - tju(1)*pju(iu+3, 1) - tju(2)*pju(iu+3, 2) - tju(3)*pju(iu+3, 3)
      If (pju(iu+3,5)<pju(iu,5)) Then
        kfl(1) = -kfl(3)
        px(1) = -px(3)
        py(1) = -py(3)
        gam(1) = gam(3)
        If (in(3)/=in(6)) Then
          Do j = 1, 4
            p(in(6), j) = p(in(3), j)
            p(in(6)+1, j) = p(in(3)+1, j)
          End Do
        End If
        Do jq = 1, 2
          in(3+jq) = in(jq)
          p(in(jq)+2, 3) = p(in(jq)+2, 3) - p(in(jq)+2, 4)
          p(in(jq)+2, 1) = p(in(jq)+2, 1) - (3-2*jq)*p(in(jq)+2, 4)
        End Do
        Goto 320
      End If
      If (iabs(kfl(1))>10) Goto 270
      i = i - 1
      kfjh(iu) = kfl(1)
      Do j = 1, 4
        pju(iu+3, j) = pju(iu+3, j) - p(i+1, j)
      End Do
    End Do
    njs(jt) = i - ista
    kfjs(jt) = k(k(mju(jt+2),3), 2)
    kfls = 2*int(rlu(0)+3.*parj(4)/(1.+3.*parj(4))) + 1
    If (kfjh(1)==kfjh(2)) kfls = 3
    If (ista/=i) kfjs(jt) = isign(1000*max(iabs(kfjh(1)),iabs(kfjh(2)))+100*min(iabs(kfjh(1)),iabs(kfjh(2)))+kfls, kfjh(1))
    Do j = 1, 4
      pjs(jt, j) = pju(1, j) + pju(2, j) + p(mju(jt), j)
      pjs(jt+2, j) = pju(4, j) + pju(5, j)
    End Do
    pjs(jt, 5) = sqrt(max(0.,pjs(jt,4)**2-pjs(jt,1)**2-pjs(jt,2)**2-pjs(jt,3)**2))
  490 End Do
  500 If (mju(1)/=0 .And. mju(2)/=0) Then
    ns = mju(2) - mju(1)
    nb = mju(1) - n
  Else If (mju(1)/=0) Then
    ns = n + nr - mju(1)
    nb = mju(1) - n
  Else If (mju(2)/=0) Then
    ns = mju(2) - n
    nb = 1
  Else If (iabs(k(n+1,2))/=21) Then
    ns = nr - 1
    nb = 1
  Else
    ns = nr + 1
    w2sum = 0.
    Do is = 1, nr
      p(n+nr+is, 1) = 0.5*four(n+is, n+is+1-nr*(is/nr))
      w2sum = w2sum + p(n+nr+is, 1)
    End Do
    w2ran = rlu(0)*w2sum
    nb = 0
    520 nb = nb + 1
    w2sum = w2sum - p(n+nr+nb, 1)
    If (w2sum>w2ran .And. nb<nr) Goto 520
  End If
  Do is = 1, ns
    is1 = n + is + nb - 1 - nr*((is+nb-2)/nr)
    is2 = n + is + nb - nr*((is+nb-1)/nr)
    Do j = 1, 5
      dp(1, j) = dble(p(is1,j))
      If (iabs(k(is1,2))==21) dp(1, j) = 0.5D0*dp(1, j)
      If (is1==mju(1)) dp(1, j) = dble(pjs(1,j)-pjs(3,j))
      dp(2, j) = dble(p(is2,j))
      If (iabs(k(is2,2))==21) dp(2, j) = 0.5D0*dp(2, j)
      If (is2==mju(2)) dp(2, j) = dble(pjs(2,j)-pjs(4,j))
    End Do
    dp(3, 5) = dfour(1, 1)
    dp(4, 5) = dfour(2, 2)
    dhkc = dfour(1, 2)
    If (dp(3,5)+2.D0*dhkc+dp(4,5)<=0.D0) Then
      dp(3, 5) = dp(1, 5)**2
      dp(4, 5) = dp(2, 5)**2
      dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2+dp(1,5)**2)
      dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2+dp(2,5)**2)
      dhkc = dfour(1, 2)
    End If
    dhks = sqrt(dhkc**2-dp(3,5)*dp(4,5))
    dhk1 = 0.5D0*((dp(4,5)+dhkc)/dhks-1.D0)
    dhk2 = 0.5D0*((dp(3,5)+dhkc)/dhks-1.D0)
    in1 = n + nr + 4*is - 3
    p(in1, 5) = sqrt(sngl(dp(3,5)+2.D0*dhkc+dp(4,5)))
    Do j = 1, 4
      p(in1, j) = sngl((1.D0+dhk1)*dp(1,j)-dhk2*dp(2,j))
      p(in1+1, j) = sngl((1.D0+dhk2)*dp(2,j)-dhk1*dp(1,j))
    End Do
  End Do
  isav = i
  550 ntry = ntry + 1
  If (ntry>100 .And. ntryr<=4) Then
    paru12 = 4.*paru12
    paru13 = 2.*paru13
    Goto 130
  Else If (ntry>100) Then
    Call luerrm(14, '(LUSTRF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = isav
  Do j = 1, 4
    p(n+nrs, j) = 0.
    Do is = 1, nr
      p(n+nrs, j) = p(n+nrs, j) + p(n+is, j)
    End Do
  End Do
  Do jt = 1, 2
    irank(jt) = 0
    If (mju(jt)/=0) irank(jt) = njs(jt)
    If (ns>nr) irank(jt) = 1
    ie(jt) = k(n+1+(jt/2)*(np-1), 3)
    in(3*jt+1) = n + nr + 1 + 4*(jt/2)*(ns-1)
    in(3*jt+2) = in(3*jt+1) + 1
    in(3*jt+3) = n + nr + 4*ns + 2*jt - 1
    Do in1 = n + nr + 2 + jt, n + nr + 4*ns - 2 + jt, 4
      p(in1, 1) = 2 - jt
      p(in1, 2) = jt - 1
      p(in1, 3) = 1.
    End Do
  End Do
  If (ns<nr) Then
    px(1) = 0.
    py(1) = 0.
    If (ns==1 .And. mju(1)+mju(2)==0) Call luptdi(0, px(1), py(1))
    px(2) = -px(1)
    py(2) = -py(1)
    Do jt = 1, 2
      kfl(jt) = k(ie(jt), 2)
      If (mju(jt)/=0) kfl(jt) = kfjs(jt)
      mstj(93) = 1
      pmq(jt) = ulmass(kfl(jt))
      gam(jt) = 0.
    End Do
  Else
    kfl(3) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
    Call lukfdi(kfl(3), 0, kfl(1), kdump)
    kfl(2) = -kfl(1)
    If (iabs(kfl(1))>10 .And. rlu(0)>0.5) Then
      kfl(2) = -(kfl(1)+isign(10000,kfl(1)))
    Else If (iabs(kfl(1))>10) Then
      kfl(1) = -(kfl(2)+isign(10000,kfl(2)))
    End If
    Call luptdi(kfl(1), px(1), py(1))
    px(2) = -px(1)
    py(2) = -py(1)
    pr3 = min(25., 0.1*p(n+nr+1,5)**2)
    590 Call luzdis(kfl(1), kfl(2), pr3, z)
    zr = pr3/(z*p(n+nr+1,5)**2)
    If (zr>=1.) Goto 590
    Do jt = 1, 2
      mstj(93) = 1
      pmq(jt) = ulmass(kfl(jt))
      gam(jt) = pr3*(1.-z)/z
      in1 = n + nr + 3 + 4*(jt/2)*(ns-1)
      p(in1, jt) = 1. - z
      p(in1, 3-jt) = jt - 1
      p(in1, 3) = (2-jt)*(1.-z) + (jt-1)*z
      p(in1+1, jt) = zr
      p(in1+1, 3-jt) = 2 - jt
      p(in1+1, 3) = (2-jt)*(1.-zr) + (jt-1)*zr
    End Do
  End If
  Do jt = 1, 2
    If (jt==1 .Or. ns==nr-1) Then
      in1 = in(3*jt+1)
      in3 = in(3*jt+3)
      Do j = 1, 4
        dp(1, j) = dble(p(in1,j))
        dp(2, j) = dble(p(in1+1,j))
        dp(3, j) = 0.D0
        dp(4, j) = 0.D0
      End Do
      dp(1, 4) = dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
      dp(2, 4) = dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
      dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
      dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
      dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
      If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1.D0
      If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1.D0
      If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1.D0
      If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1.D0
      dhc12 = dfour(1, 2)
      dhcx1 = dfour(3, 1)/dhc12
      dhcx2 = dfour(3, 2)/dhc12
      dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
      dhcy1 = dfour(4, 1)/dhc12
      dhcy2 = dfour(4, 2)/dhc12
      dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
      dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
      Do j = 1, 4
        dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
        p(in3, j) = sngl(dp(3,j))
        p(in3+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
      End Do
    Else
      Do j = 1, 4
        p(in3+2, j) = p(in3, j)
        p(in3+3, j) = p(in3+1, j)
      End Do
    End If
  End Do
  If (mju(1)+mju(2)>0) Then
    Do jt = 1, 2
      If (njs(jt)==0) Goto 660
      Do j = 1, 4
        p(n+nrs, j) = p(n+nrs, j) - pjs(jt+2, j)
      End Do
    660 End Do
  End If
  670 i = i + 1
  If (2*i-nsav>=mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  jt = int(1.5+rlu(0))
  If (iabs(kfl(3-jt))>10) jt = 3 - jt
  jr = 3 - jt
  js = 3 - 2*jt
  irank(jt) = irank(jt) + 1
  k(i, 1) = 1
  k(i, 3) = ie(jt)
  k(i, 4) = 0
  k(i, 5) = 0
  680 Call lukfdi(kfl(jt), 0, kfl(3), k(i,2))
  If (k(i,2)==0) Goto 550
  If (mstj(12)>=3 .And. irank(jt)==1 .And. iabs(kfl(jt))<=10 .And. iabs(kfl(3))>10) Then
    If (rlu(0)>parj(19)) Goto 680
  End If
  p(i, 5) = ulmass(k(i,2))
  Call luptdi(kfl(jt), px(3), py(3))
  pr(jt) = p(i, 5)**2 + (px(jt)+px(3))**2 + (py(jt)+py(3))**2
  mstj(93) = 1
  pmq(3) = ulmass(kfl(3))
  wmin = parj(32+mstj(11)) + pmq(1) + pmq(2) + parj(36)*pmq(3)
  If (iabs(kfl(jt))>10 .And. iabs(kfl(3))>10) wmin = wmin - 0.5*parj(36)*pmq(3)
  wrem2 = four(n+nrs, n+nrs)
  If (wrem2<0.10) Goto 550
  If (wrem2<max(wmin*(1.+(2.*rlu(0)-1.)*parj(37)),parj(32)+pmq(1)+pmq(2))**2) Goto 810
  Call luzdis(kfl(jt), kfl(3), pr(jt), z)
  kfl1a = iabs(kfl(1))
  kfl2a = iabs(kfl(2))
  If (max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),mod(kfl2a/1000,10))>=4) Then
    pr(jr) = (pmq(jr)+pmq(3))**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2
    pw12 = sqrt(max(0.,(wrem2-pr(1)-pr(2))**2-4.*pr(1)*pr(2)))
    z = (wrem2+pr(jt)-pr(jr)+pw12*(2.*z-1.))/(2.*wrem2)
    pr(jr) = (pmq(jr)+parj(32+mstj(11)))**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2
    If ((1.-z)*(wrem2-pr(jt)/z)<pr(jr)) Goto 810
  End If
  gam(3) = (1.-z)*(gam(jt)+pr(jt)/z)
  Do j = 1, 3
    in(j) = in(3*jt+j)
  End Do
  If (in(1)+1==in(2) .And. z*p(in(1)+2,3)*p(in(2)+2,3)*p(in(1),5)**2>=pr(jt)) Then
    p(in(jt)+2, 4) = z*p(in(jt)+2, 3)
    p(in(jr)+2, 4) = pr(jt)/(p(in(jt)+2,4)*p(in(1),5)**2)
    Do j = 1, 4
      p(i, j) = (px(jt)+px(3))*p(in(3), j) + (py(jt)+py(3))*p(in(3)+1, j)
    End Do
    Goto 770
  Else If (in(1)+1==in(2)) Then
    p(in(jr)+2, 4) = p(in(jr)+2, 3)
    p(in(jr)+2, jt) = 1.
    in(jr) = in(jr) + 4*js
    If (js*in(jr)>js*in(4*jr)) Goto 550
    If (four(in(1),in(2))<=1E-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
    End If
  End If
  710 If (js*in(1)>js*in(3*jr+1) .Or. js*in(2)>js*in(3*jr+2) .Or. in(1)>in(2)) Goto 550
  If (in(1)/=in(3*jt+1) .Or. in(2)/=in(3*jt+2)) Then
    Do j = 1, 4
      dp(1, j) = dble(p(in(1),j))
      dp(2, j) = dble(p(in(2),j))
      dp(3, j) = 0.D0
      dp(4, j) = 0.D0
    End Do
    dp(1, 4) = dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
    dp(2, 4) = dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
    dhc12 = dfour(1, 2)
    If (dhc12<=1D-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
      Goto 710
    End If
    in(3) = n + nr + 4*ns + 5
    dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
    dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
    dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
    If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1.D0
    If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1.D0
    If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1.D0
    If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1.D0
    dhcx1 = dfour(3, 1)/dhc12
    dhcx2 = dfour(3, 2)/dhc12
    dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
    dhcy1 = dfour(4, 1)/dhc12
    dhcy2 = dfour(4, 2)/dhc12
    dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
    dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
    Do j = 1, 4
      dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
      p(in(3), j) = sngl(dp(3,j))
      p(in(3)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
    End Do
    pxp = -(px(3)*four(in(3*jt+3),in(3))+py(3)*four(in(3*jt+3)+1,in(3)))
    pyp = -(px(3)*four(in(3*jt+3),in(3)+1)+py(3)*four(in(3*jt+3)+1,in(3)+1))
    If (abs(pxp**2+pyp**2-px(3)**2-py(3)**2)<0.01) Then
      px(3) = pxp
      py(3) = pyp
    End If
  End If
  Do j = 1, 4
    dhg(j) = 0.D0
    p(i, j) = px(jt)*p(in(3*jt+3), j) + py(jt)*p(in(3*jt+3)+1, j) + px(3)*p(in(3), j) + py(3)*p(in(3)+1, j)
    Do in1 = in(3*jt+1), in(1) - 4*js, 4*js
      p(i, j) = p(i, j) + p(in1+2, 3)*p(in1, j)
    End Do
    Do in2 = in(3*jt+2), in(2) - 4*js, 4*js
      p(i, j) = p(i, j) + p(in2+2, 3)*p(in2, j)
    End Do
  End Do
  dhm(1) = dble(four(i,i))
  dhm(2) = dble(2.*four(i,in(1)))
  dhm(3) = dble(2.*four(i,in(2)))
  dhm(4) = dble(2.*four(in(1),in(2)))
  Do in2 = in(1) + 1, in(2), 4
    Do in1 = in(1), in2 - 1, 4
      dhc = dble(2.*four(in1,in2))
      dhg(1) = dhg(1) + dble(p(in1+2,jt)*p(in2+2,jt))*dhc
      If (in1==in(1)) dhg(2) = dhg(2) - dble(float(js)*p(in2+2,jt))*dhc
      If (in2==in(2)) dhg(3) = dhg(3) + dble(float(js)*p(in1+2,jt))*dhc
      If (in1==in(1) .And. in2==in(2)) dhg(4) = dhg(4) - dhc
    End Do
  End Do
  dhs1 = dhm(jr+1)*dhg(4) - dhm(4)*dhg(jr+1)
  If (dabs(dhs1)<1D-4) Goto 550
  dhs2 = dhm(4)*(dble(gam(3))-dhg(1)) - dhm(jt+1)*dhg(jr+1) - dhg(4)*(dble(p(i,5))**2-dhm(1)) + dhg(jt+1)*dhm(jr+1)
  dhs3 = dhm(jt+1)*(dble(gam(3))-dhg(1)) - dhg(jt+1)*(dble(p(i,5))**2-dhm(1))
  p(in(jr)+2, 4) = 0.5*sngl((sqrt(max(0D0,dhs2**2-4.D0*dhs1*dhs3)))/abs(dhs1)-dhs2/dhs1)
  If (dhm(jt+1)+dhm(4)*dble(p(in(jr)+2,4))<=0.D0) Goto 550
  p(in(jt)+2, 4) = (p(i,5)**2-sngl(dhm(1))-sngl(dhm(jr+1))*p(in(jr)+2,4))/(sngl(dhm(jt+1))+sngl(dhm(4))*p(in(jr)+2,4))
  If (p(in(jr)+2,4)>p(in(jr)+2,3)) Then
    p(in(jr)+2, 4) = p(in(jr)+2, 3)
    p(in(jr)+2, jt) = 1.
    in(jr) = in(jr) + 4*js
    If (js*in(jr)>js*in(4*jr)) Goto 550
    If (four(in(1),in(2))<=1E-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
    End If
    Goto 710
  Else If (p(in(jt)+2,4)>p(in(jt)+2,3)) Then
    p(in(jt)+2, 4) = p(in(jt)+2, 3)
    p(in(jt)+2, jt) = 0.
    in(jt) = in(jt) + 4*js
    Goto 710
  End If
  770 Do j = 1, 4
    p(i, j) = p(i, j) + p(in(1)+2, 4)*p(in(1), j) + p(in(2)+2, 4)*p(in(2), j)
    p(n+nrs, j) = p(n+nrs, j) - p(i, j)
  End Do
  If (p(i,4)<=0.) Goto 550
  kfl(jt) = -kfl(3)
  pmq(jt) = pmq(3)
  px(jt) = -px(3)
  py(jt) = -py(3)
  gam(jt) = gam(3)
  If (in(3)/=in(3*jt+3)) Then
    Do j = 1, 4
      p(in(3*jt+3), j) = p(in(3), j)
      p(in(3*jt+3)+1, j) = p(in(3)+1, j)
    End Do
  End If
  Do jq = 1, 2
    in(3*jt+jq) = in(jq)
    p(in(jq)+2, 3) = p(in(jq)+2, 3) - p(in(jq)+2, 4)
    p(in(jq)+2, jt) = p(in(jq)+2, jt) - js*(3-2*jq)*p(in(jq)+2, 4)
  End Do
  Goto 670
  810 i = i + 1
  k(i, 1) = 1
  k(i, 3) = ie(jr)
  k(i, 4) = 0
  k(i, 5) = 0
  Call lukfdi(kfl(jr), -kfl(3), kfldmp, k(i,2))
  If (k(i,2)==0) Goto 550
  p(i, 5) = ulmass(k(i,2))
  pr(jr) = p(i, 5)**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2
  jq = 1
  If (p(in(4)+2,3)*p(in(5)+2,3)*four(in(4),in(5))<p(in(7),3)*p(in(8),3)*four(in(7),in(8))) jq = 2
  dhc12 = dble(four(in(3*jq+1),in(3*jq+2)))
  dhr1 = dble(four(n+nrs,in(3*jq+2)))/dhc12
  dhr2 = dble(four(n+nrs,in(3*jq+1)))/dhc12
  If (in(4)/=in(7) .Or. in(5)/=in(8)) Then
    px(3-jq) = -four(n+nrs, in(3*jq+3)) - px(jq)
    py(3-jq) = -four(n+nrs, in(3*jq+3)+1) - py(jq)
    pr(3-jq) = p(i+(jt+jq-3)**2-1, 5)**2 + (px(3-jq)+(2*jq-3)*js*px(3))**2 + (py(3-jq)+(2*jq-3)*js*py(3))**2
  End If
  wrem2 = wrem2 + (px(1)+px(2))**2 + (py(1)+py(2))**2
  fd = (sqrt(pr(1))+sqrt(pr(2)))/sqrt(wrem2)
  If (mju(1)+mju(2)/=0 .And. i==isav+2 .And. fd>=1.) Goto 180
  If (fd>=1.) Goto 550
  fa = wrem2 + pr(jt) - pr(jr)
  If (mstj(11)==2) prev = 0.5*fd**parj(37+mstj(11))
  If (mstj(11)/=2) prev = 0.5*exp(max(-100.,log(fd)*parj(37+mstj(11))*(pr(1)+pr(2))**2))
  fb = sign(sqrt(max(0.,fa**2-4.*wrem2*pr(jt))), js*(rlu(0)-prev))
  kfl1a = iabs(kfl(1))
  kfl2a = iabs(kfl(2))
  If (max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),mod(kfl2a/1000,10))>=6) fb = sign(sqrt(max(0.,fa**2-4.*wrem2*pr(jt))), float(js))
  Do j = 1, 4
    p(i-1, j) = (px(jt)+px(3))*p(in(3*jq+3), j) + (py(jt)+py(3))*p(in(3*jq+3)+1, j) + 0.5*(sngl(dhr1)*(fa+fb)*p(in(3*jq+1),j)+sngl(dhr2)*(fa-fb)*p(in(3*jq+2),j))/wrem2
    p(i, j) = p(n+nrs, j) - p(i-1, j)
  End Do
  n = i - nrs + 1
  Do i = nsav + 1, nsav + np
    im = k(i, 3)
    k(im, 1) = k(im, 1) + 10
    If (mstu(16)/=2) Then
      k(im, 4) = nsav + 1
      k(im, 5) = nsav + 1
    Else
      k(im, 4) = nsav + 2
      k(im, 5) = n
    End If
  End Do
  nsav = nsav + 1
  k(nsav, 1) = 11
  k(nsav, 2) = 92
  k(nsav, 3) = ip
  k(nsav, 4) = nsav + 1
  k(nsav, 5) = n
  Do j = 1, 4
    p(nsav, j) = sngl(dps(j))
    v(nsav, j) = v(ip, j)
  End Do
  p(nsav, 5) = sqrt(sngl(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2)))
  v(nsav, 5) = 0.
  Do i = nsav + 1, n
    Do j = 1, 5
      k(i, j) = k(i+nrs-1, j)
      p(i, j) = p(i+nrs-1, j)
      v(i, j) = 0.
    End Do
  End Do
  Do i = nsav + 1, n
    Do j = 1, 5
      k(i-nsav+n, j) = k(i, j)
      p(i-nsav+n, j) = p(i, j)
    End Do
  End Do
  i1 = nsav
  Do i = n + 1, 2*n - nsav
    If (k(i,3)/=ie(1)) Goto 880
    i1 = i1 + 1
    Do j = 1, 5
      k(i1, j) = k(i, j)
      p(i1, j) = p(i, j)
    End Do
    If (mstu(16)/=2) k(i1, 3) = nsav
  880 End Do
  Do i = 2*n - nsav, n + 1, -1
    If (k(i,3)==ie(1)) Goto 900
    i1 = i1 + 1
    Do j = 1, 5
      k(i1, j) = k(i, j)
      p(i1, j) = p(i, j)
    End Do
    If (mstu(16)/=2) k(i1, 3) = nsav
  900 End Do
  Call ludbrb(nsav+1, n, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))
  Do i = nsav + 1, n
    Do j = 1, 4
      v(i, j) = v(ip, j)
    End Do
  End Do
  Return
End Subroutine lustrf
