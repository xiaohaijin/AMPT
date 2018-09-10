Subroutine m1535(lb1, lb2, srt, x1535)
  Save
  s0 = 2.424
  x1535 = 0.
  If (srt<=s0) Return
  sigma = 2.*0.102*(srt-s0)/(0.058+(srt-s0)**2)
  If ((lb1*lb2==18 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==6 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==8 .And. (lb1==1 .Or. lb2==1))) Then
     x1535 = sigma
     Return
  End If
  If (lb1*lb2==7) Then
     x1535 = 3.*sigma
     Return
  End If
  If ((lb1*lb2==11) .Or. (lb1*lb2==20 .And. (lb1==2 .Or. lb2==2))) Then
     x1535 = sigma
     Return
  End If
  If ((lb1*lb2==10 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==22 .And. (lb1==2 .Or. lb2==2))) x1535 = 3.*sigma
  Return
End Subroutine m1535
