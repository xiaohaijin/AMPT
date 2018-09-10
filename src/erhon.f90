Real Function erhon(em1, em2, lb1, lb2, srt)
  Dimension arrayj(19), arrayl(19), arraym(19), arrayw(19), arrayb(19)
  Save
  Data arrayj/0.5, 1.5, 0.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/1, 2, 0, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.44, 1.52, 1.535, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.2, 0.125, 0.15, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.15, 0.20, 0.05, 0.175, 0.025, 0.125, 0.1, 0.20, 0.53, 0.34, 0.05, 0.07, 0.15, 0.45, 0.45, 0.058, 0.08, 0.12, 0.08/
  pi = 3.1415926
  xs = 0
  Do ir = 1, 19
    If (ir<=8) Then
      If (((lb1*lb2==27 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==25*2 .And. (lb1==2 .Or. lb2==2))) .Or. ((lb1*lb2==-25 .And. (lb1==-1 .Or. lb2==-1)) .Or. (lb1*lb2==-27*2 .And. (lb1==-2 .Or. lb2==-2)))) branch = 0.
      If ((iabs(lb1*lb2)==26 .And. (iabs(lb1)==1 .Or. iabs(lb2)==1)) .Or. (iabs(lb1*lb2)==26*2 .And. (iabs(lb1)==2 .Or. iabs(lb2)==2))) branch = 1./3.
      If (((lb1*lb2==27*2 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==25 .And. (lb1==1 .Or. lb2==1))) .Or. ((lb1*lb2==-25*2 .And. (lb1==-2 .Or. lb2==-2)) .Or. (lb1*lb2==-27 .And. (lb1==-1 .Or. lb2==-1)))) branch = 2./3.
    Else
      If (((lb1*lb2==27 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==25*2 .And. (lb1==2 .Or. lb2==2))) .Or. ((lb1*lb2==-25 .And. (lb1==-1 .Or. lb2==-1)) .Or. (lb1*lb2==-27*2 .And. (lb1==-2 .Or. lb2==-2)))) branch = 1.
      If ((iabs(lb1*lb2)==26 .And. (iabs(lb1)==1 .Or. iabs(lb2)==1)) .Or. (iabs(lb1*lb2)==26*2 .And. (iabs(lb1)==2 .Or. iabs(lb2)==2))) branch = 2./3.
      If (((lb1*lb2==27*2 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==25 .And. (lb1==1 .Or. lb2==1))) .Or. ((lb1*lb2==-25*2 .And. (lb1==-2 .Or. lb2==-2)) .Or. (lb1*lb2==-27 .And. (lb1==-1 .Or. lb2==-1)))) branch = 1./3.
    End If
    xs0 = fdr(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), srt, em1, em2)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  erhon = xs
  Return
End Function erhon
