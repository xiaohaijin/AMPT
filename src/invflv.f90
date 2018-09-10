Function invflv(iart)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  If (iart==-6) Then
    invflv = -1114
    Return
  End If
  If (iart==-7) Then
    invflv = -2114
    Return
  End If
  If (iart==-8) Then
    invflv = -2214
    Return
  End If
  If (iart==-9) Then
    invflv = -2224
    Return
  End If
  If (iart==-1) Then
    invflv = -2212
    Return
  End If
  If (iart==-2) Then
    invflv = -2112
    Return
  End If
  If (iart==0) Then
    invflv = 221
    Return
  End If
  If (iart==1) Then
    invflv = 2212
    Return
  End If
  If (iart==2) Then
    invflv = 2112
    Return
  End If
  If (iart==3) Then
    invflv = -211
    Return
  End If
  If (iart==4) Then
    invflv = 111
    Return
  End If
  If (iart==5) Then
    invflv = 211
    Return
  End If
  If (iart==6) Then
    invflv = 1114
    Return
  End If
  If (iart==7) Then
    invflv = 2114
    Return
  End If
  If (iart==8) Then
    invflv = 2214
    Return
  End If
  If (iart==9) Then
    invflv = 2224
    Return
  End If
  If (iart==14) Then
    invflv = 3122
    Return
  End If
  If (iart==-14) Then
    invflv = -3122
    Return
  End If
  If (iart==15) Then
    invflv = 3112
    Return
  End If
  If (iart==-15) Then
    invflv = -3112
    Return
  End If
  If (iart==16) Then
    invflv = 3212
    Return
  End If
  If (iart==-16) Then
    invflv = -3212
    Return
  End If
  If (iart==17) Then
    invflv = 3222
    Return
  End If
  If (iart==-17) Then
    invflv = -3222
    Return
  End If
  If (iart==21) Then
    invflv = -321
    Return
  End If
  If (iart==23) Then
    invflv = 321
    Return
  End If
  If (iart==22) Then
    invflv = 130
    Return
  End If
  If (iart==24) Then
    invflv = 310
    Return
  End If
  If (iart==25) Then
    invflv = -213
    Return
  End If
  If (iart==26) Then
    invflv = 113
    Return
  End If
  If (iart==27) Then
    invflv = 213
    Return
  End If
  If (iart==28) Then
    invflv = 223
    Return
  End If
  If (iart==29) Then
    invflv = 333
    Return
  End If
  If (iart==30) Then
    invflv = 323
    If (ranart(nseed)>0.5) invflv = 313
    Return
  End If
  If (iart==-30) Then
    invflv = -323
    If (ranart(nseed)>0.5) invflv = -313
    Return
  End If
  If (iart==31) Then
    invflv = 331
    Return
  End If
  If (iart==32) Then
    invflv = 777
    Return
  End If
  If (iart==40) Then
    invflv = 3312
    Return
  End If
  If (iart==-40) Then
    invflv = -3312
    Return
  End If
  If (iart==41) Then
    invflv = 3322
    Return
  End If
  If (iart==-41) Then
    invflv = -3322
    Return
  End If
  If (iart==45) Then
    invflv = 3334
    Return
  End If
  If (iart==-45) Then
    invflv = -3334
    Return
  End If
  If (iart==44) Then
    invflv = 6666
    Return
  End If
  If (iart==42) Then
    invflv = 42
    Return
  Else If (iart==-42) Then
    invflv = -42
    Return
  End If
  invflv = iart - 10000
  Return
End Function invflv
