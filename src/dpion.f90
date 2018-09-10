Real Function dpion(em1, em2, lb1, lb2, srt)
  Dimension arrayj(19), arrayl(19), arraym(19), arrayw(19), arrayb(19)
  Save
  Data arrayj/0.5, 1.5, 0.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/1, 2, 0, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.44, 1.52, 1.535, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.2, 0.125, 0.15, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.15, 0.25, 0., 0.05, 0.575, 0.125, 0.379, 0.10, 0.10, 0.062, 0.45, 0.60, 0.6984, 0.05, 0.25, 0.089, 0.19, 0.2, 0.13/
  pi = 3.1415926
  amn = 0.94
  amp = 0.14
  xs = 0
  Do ir = 1, 19
    branch = 0.
    If (ir<=8) Then
      If (((lb1*lb2==5*7 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*8 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*7 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*8 .And. (lb1==5 .Or. lb2==5)))) branch = 1./6.
      If ((iabs(lb1*lb2)==4*7 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*8 .And. (lb1==4 .Or. lb2==4))) branch = 1./3.
      If (((lb1*lb2==5*6 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*9 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*6 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*9 .And. (lb1==5 .Or. lb2==5)))) branch = 1./2.
    Else
      If (((lb1*lb2==5*8 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==5*6 .And. (lb1==5 .Or. lb2==5))) .Or. ((lb1*lb2==-3*8 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-3*6 .And. (lb1==3 .Or. lb2==3)))) branch = 2./5.
      If (((lb1*lb2==3*9 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==3*7 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-5*9 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-5*7 .And. (lb1==5 .Or. lb2==5)))) branch = 2./5.
      If (((lb1*lb2==5*7 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*8 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*7 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*8 .And. (lb1==5 .Or. lb2==5)))) branch = 8./15.
      If ((iabs(lb1*lb2)==4*7 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*8 .And. (lb1==4 .Or. lb2==4))) branch = 1./15.
      If ((iabs(lb1*lb2)==4*9 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*6 .And. (lb1==4 .Or. lb2==4))) branch = 3./5.
    End If
    xs0 = fd2(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), em1, em2, srt)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  dpion = xs
  Return
End Function dpion
