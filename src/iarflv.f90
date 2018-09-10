Function iarflv(ipdg)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  If (ipdg==-1114) Then
    iarflv = -6
    Return
  End If
  If (ipdg==-2114) Then
    iarflv = -7
    Return
  End If
  If (ipdg==-2214) Then
    iarflv = -8
    Return
  End If
  If (ipdg==-2224) Then
    iarflv = -9
    Return
  End If
  If (ipdg==-2212) Then
    iarflv = -1
    Return
  End If
  If (ipdg==-2112) Then
    iarflv = -2
    Return
  End If
  If (ipdg==221) Then
    iarflv = 0
    Return
  End If
  If (ipdg==2212) Then
    iarflv = 1
    Return
  End If
  If (ipdg==2112) Then
    iarflv = 2
    Return
  End If
  If (ipdg==-211) Then
    iarflv = 3
    Return
  End If
  If (ipdg==111) Then
    iarflv = 4
    Return
  End If
  If (ipdg==211) Then
    iarflv = 5
    Return
  End If
  If (ipdg==1114) Then
    iarflv = 6
    Return
  End If
  If (ipdg==2114) Then
    iarflv = 7
    Return
  End If
  If (ipdg==2214) Then
    iarflv = 8
    Return
  End If
  If (ipdg==2224) Then
    iarflv = 9
    Return
  End If
  If (ipdg==3122) Then
    iarflv = 14
    Return
  End If
  If (ipdg==-3122) Then
    iarflv = -14
    Return
  End If
  If (ipdg==3112) Then
    iarflv = 15
    Return
  End If
  If (ipdg==-3112) Then
    iarflv = -15
    Return
  End If
  If (ipdg==3212) Then
    iarflv = 16
    Return
  End If
  If (ipdg==-3212) Then
    iarflv = -16
    Return
  End If
  If (ipdg==3222) Then
    iarflv = 17
    Return
  End If
  If (ipdg==-3222) Then
    iarflv = -17
    Return
  End If
  If (ipdg==-321) Then
    iarflv = 21
    Return
  End If
  If (ipdg==321) Then
    iarflv = 23
    Return
  End If
  If (ipdg==311) Then
    iarflv = 23
    Return
  End If
  If (ipdg==-311) Then
    iarflv = 21
    Return
  End If
  If (ipdg==310 .Or. ipdg==130) Then
    r = ranart(nseed)
    If (r>0.5) Then
      iarflv = 23
    Else
      iarflv = 21
    End If
    Return
  End If
  If (ipdg==-213) Then
    iarflv = 25
    Return
  End If
  If (ipdg==113) Then
    iarflv = 26
    Return
  End If
  If (ipdg==213) Then
    iarflv = 27
    Return
  End If
  If (ipdg==223) Then
    iarflv = 28
    Return
  End If
  If (ipdg==333) Then
    iarflv = 29
    Return
  End If
  If (ipdg==323) Then
    iarflv = 30
    Return
  End If
  If (ipdg==-323) Then
    iarflv = -30
    Return
  End If
  If (ipdg==313) Then
    iarflv = 30
    Return
  End If
  If (ipdg==-313) Then
    iarflv = -30
    Return
  End If
  If (ipdg==331) Then
    iarflv = 31
    Return
  End If
  If (ipdg==3312) Then
    iarflv = 40
    Return
  End If
  If (ipdg==-3312) Then
    iarflv = -40
    Return
  End If
  If (ipdg==3322) Then
    iarflv = 41
    Return
  End If
  If (ipdg==-3322) Then
    iarflv = -41
    Return
  End If
  If (ipdg==3334) Then
    iarflv = 45
    Return
  End If
  If (ipdg==-3334) Then
    iarflv = -45
    Return
  End If
  If (ipdg==6666) Then
    iarflv = 44
    Return
  End If
  If (ipdg==42 .Or. ipdg==-42) Then
    iarflv = ipdg
    Return
  End If
  iarflv = ipdg + 10000
  Return
End Function iarflv
