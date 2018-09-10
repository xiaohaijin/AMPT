Subroutine zpca2
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /para5/iconfg, iordsc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /ilist6/t, iopern, icolln
  Common /rndm1/number
  Common /rndm2/iff
  Common /rndm3/iseedp
  Common /arevt/iaevt, iarun, miss
  Save
  If (iconfg<=3) Then
     Call zpca2a
  Else
     Call zpca2b
  End If
  If (ioscar==1) Then
     Call zpca2c
  End If
  Write (25, *) ' Event ', iaevt, ', run ', iarun
  Write (25, *) '    number of operations = ', iopern
  Write (25, *) '    number of collisions between particles = ', icolln
  Write (25, *) '    freezeout time=', t
  Write (25, *) '    ending at the ', number, 'th random number'
  Write (25, *) '    ending collision iff=', iff
  Return
End Subroutine zpca2
