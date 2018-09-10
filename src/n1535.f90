Subroutine n1535(lb1, lb2, srt, x1535)
  Save
  s0 = 2.424
  x1535 = 0.
  If (srt<=s0) Return
  sigma = 2.*0.102*(srt-s0)/(0.058+(srt-s0)**2)
  If ((lb1*lb2==1) .Or. (lb1==2 .And. lb2==2)) Then
     x1535 = sigma
     Return
  End If
  If (lb1*lb2==2) Then
     x1535 = 3.*sigma
     Return
  End If
  If ((lb1*lb2==63 .And. (lb1==7 .Or. lb2==7)) .Or. (lb1*lb2==64 .And. (lb1==8 .Or. lb2==8)) .Or. (lb1*lb2==48 .And. (lb1==6 .Or. lb2==6)) .Or. (lb1*lb2==49 .And. (lb1==7 .Or. lb2==7))) Then
     x1535 = sigma
     Return
  End If
  If ((lb1*lb2==54 .And. (lb1==6 .Or. lb2==6)) .Or. (lb1*lb2==56 .And. (lb1==7 .Or. lb2==7))) Then
     x1535 = 3.*sigma
     Return
  End If
  If ((lb1==10 .And. lb2==10) .Or. (lb1==11 .And. lb2==11)) x1535 = sigma
  If (lb1*lb2==110 .And. (lb1==10 .Or. lb2==10)) x1535 = 3.*sigma
  Return
End Subroutine n1535
