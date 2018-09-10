Subroutine bbarfs(lbb1, lbb2, ei1, ei2, iblock, iseed)
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Common /rndf77/nseed
  Save
  rd = ranart(nseed)
  wsum = 0.
  Do i = 1, nstate
    wsum = wsum + weight(i)
    If (rd<=(wsum/wtot)) Then
      ifs = i
      ei1 = ppbm(i, 1)
      ei2 = ppbm(i, 2)
      Goto 10
    End If
  End Do
  10 Continue
  If (ifs==1) Then
    iblock = 1801
    lbb1 = -1
    lbb2 = 1
  Else If (ifs==2) Then
    If (ranart(nseed)<=0.5) Then
      iblock = 18021
      lbb1 = -1
      lbb2 = 2
    Else
      iblock = 18022
      lbb1 = 1
      lbb2 = -2
    End If
  Else If (ifs==3) Then
    iblock = 1803
    lbb1 = -2
    lbb2 = 2
  Else If (ifs==4 .Or. ifs==5) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
      If (ifs==4) Then
        iblock = 18041
        lbb1 = -1
      Else
        iblock = 18051
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.25) Then
        lbb2 = 6
      Else If (rd2<=0.5) Then
        lbb2 = 7
      Else If (rd2<=0.75) Then
        lbb2 = 8
      Else
        lbb2 = 9
      End If
    Else
      If (ifs==4) Then
        iblock = 18042
        lbb1 = 1
      Else
        iblock = 18052
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.25) Then
        lbb2 = -6
      Else If (rd2<=0.5) Then
        lbb2 = -7
      Else If (rd2<=0.75) Then
        lbb2 = -8
      Else
        lbb2 = -9
      End If
    End If
  Else If (ifs==6 .Or. ifs==7) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
      If (ifs==6) Then
        iblock = 18061
        lbb1 = -1
      Else
        iblock = 18071
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 10
      Else
        lbb2 = 11
      End If
    Else
      If (ifs==6) Then
        iblock = 18062
        lbb1 = 1
      Else
        iblock = 18072
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -10
      Else
        lbb2 = -11
      End If
    End If
  Else If (ifs==8) Then
    iblock = 1808
    rd1 = ranart(nseed)
    If (rd1<=0.25) Then
      lbb1 = 6
    Else If (rd1<=0.5) Then
      lbb1 = 7
    Else If (rd1<=0.75) Then
      lbb1 = 8
    Else
      lbb1 = 9
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.25) Then
      lbb2 = -6
    Else If (rd2<=0.5) Then
      lbb2 = -7
    Else If (rd2<=0.75) Then
      lbb2 = -8
    Else
      lbb2 = -9
    End If
  Else If (ifs==9 .Or. ifs==10) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
      If (ifs==9) Then
        iblock = 18091
        lbb1 = -1
      Else
        iblock = 18101
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 12
      Else
        lbb2 = 13
      End If
    Else
      If (ifs==9) Then
        iblock = 18092
        lbb1 = 1
      Else
        iblock = 18102
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -12
      Else
        lbb2 = -13
      End If
    End If
  Else If (ifs==11 .Or. ifs==12) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
      rd1 = ranart(nseed)
      If (rd1<=0.25) Then
        lbb1 = -6
      Else If (rd1<=0.5) Then
        lbb1 = -7
      Else If (rd1<=0.75) Then
        lbb1 = -8
      Else
        lbb1 = -9
      End If
      If (ifs==11) Then
        iblock = 18111
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = 10
        Else
          lbb2 = 11
        End If
      Else
        iblock = 18121
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = 12
        Else
          lbb2 = 13
        End If
      End If
    Else
      rd1 = ranart(nseed)
      If (rd1<=0.25) Then
        lbb1 = 6
      Else If (rd1<=0.5) Then
        lbb1 = 7
      Else If (rd1<=0.75) Then
        lbb1 = 8
      Else
        lbb1 = 9
      End If
      If (ifs==11) Then
        iblock = 18112
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = -10
        Else
          lbb2 = -11
        End If
      Else
        iblock = 18122
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = -12
        Else
          lbb2 = -13
        End If
      End If
    End If
  Else If (ifs==13) Then
    iblock = 1813
    rd1 = ranart(nseed)
    If (rd1<=0.5) Then
      lbb1 = 10
    Else
      lbb1 = 11
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.5) Then
      lbb2 = -10
    Else
      lbb2 = -11
    End If
  Else If (ifs==14) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
      iblock = 18141
      rd1 = ranart(nseed)
      If (rd1<=0.5) Then
        lbb1 = -10
      Else
        lbb1 = -11
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 12
      Else
        lbb2 = 13
      End If
    Else
      iblock = 18142
      rd1 = ranart(nseed)
      If (rd1<=0.5) Then
        lbb1 = 10
      Else
        lbb1 = 11
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -12
      Else
        lbb2 = -13
      End If
    End If
  Else If (ifs==15) Then
    iblock = 1815
    rd1 = ranart(nseed)
    If (rd1<=0.5) Then
      lbb1 = 12
    Else
      lbb1 = 13
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.5) Then
      lbb2 = -12
    Else
      lbb2 = -13
    End If
  Else
  End If
  Return
End Subroutine bbarfs
