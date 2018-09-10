Subroutine zpcmn
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Save
  Do i = 1, nevnt
     ievt = i
     Call inievt
     Do j = 1, nsbrun
        isbrun = j
        Call inirun
3000    Continue
        Call zpcrun(*4000)
        Call zpca1
        Goto 3000
4000    Continue
        Call zpca2
     End Do
  End Do
  Call zpcou
  Call zpstrg
  Return
End Subroutine zpcmn
