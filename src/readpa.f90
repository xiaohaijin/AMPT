Subroutine readpa
  Implicit Double Precision (A-H, O-Z)
  External ran1
  Character *50 str
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Common /para5/iconfg, iordsc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /rndm1/number
  Common /rndm2/iff
  Common /rndm3/iseedp
  Save
  iseed = iseedp
  Open (25, File='../data/zpc.res', Status='unknown')
  If (ioscar==1) Then
     Open (26, File='../data/parton.oscar', Status='unknown')
     Open (19, File='../data/hadron.oscar', Status='unknown')
  End If
  xmp = 0D0
  cutof2 = 4.5D0*(alpha/xmu)**2
  rscut2 = 0.01D0
  nsevt = 1
  nevnt = 1
  nsbrun = 1
  iftflg = 0
  ireflg = 1
  If (ireflg==0) Then
     Open (27, File='zpc.inp', Status='UNKNOWN')
  End If
  igeflg = 0
  ibstfg = 0
  iconfg = 1
  iordsc = 11
  v1 = 0.2D0
  v2 = 0.2D0
  v3 = 0.2D0
  size1 = 1.5D0
  size2 = 1.5D0
  size3 = 0.7D0
  If (size1==0D0 .Or. size2==0D0 .Or. size3==0D0) Then
     If (size1/=0D0 .Or. size2/=0D0 .Or. size3/=0D0 .Or. v1/=0D0 .Or. v2/=0D0 .Or. v3/=0D0) Then
        Print *, 'to get rid of space division:'
        Print *, 'set all sizes and vs to 0'
        Stop 'chker'
     End If
  End If
  size = min(size1, size2, size3)
  iff = -1
  isedng = -iseed
  a = ran1(isedng)
  irused = 2
  Do i = 1, irused - 1
     iseed2 = 2
     a = ran1(iseed2)
  End Do
  If (iconfg==2 .Or. iconfg==3) Then
     v1 = 0D0
     v2 = 0D0
  End If
  If (iconfg==4 .Or. iconfg==5) Then
     v1 = 0D0
     v2 = 0D0
     v3 = 0D0
  End If
  Close (5)
  Return
End Subroutine readpa
