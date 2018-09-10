Subroutine getnp
  Parameter (maxstr=150001)
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Save
  If (natt==0) Then
     npart1 = 0
     npart2 = 0
     Return
  End If
  pzproj = sqrt(hint1(6)**2-hint1(8)**2)
  pztarg = sqrt(hint1(7)**2-hint1(9)**2)
  epsipz = 0.01
  epsipt = 1E-6
  nspec1 = 0
  nspec2 = 0
  Do i = 1, natt
     If ((katt(i,1)==2112 .Or. katt(i,1)==2212) .And. abs(patt(i,1))<=epsipt .And. abs(patt(i,2))<=epsipt) Then
        If (patt(i,3)>amax1(0.,pzproj-epsipz)) Then
           nspec1 = nspec1 + 1
        Else If (patt(i,3)<(-pztarg+epsipz)) Then
           nspec2 = nspec2 + 1
        End If
     End If
  End Do
  npart1 = ihnt2(1) - nspec1
  npart2 = ihnt2(3) - nspec2
  Return
End Subroutine getnp
