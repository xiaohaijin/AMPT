Subroutine arini
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Save
  iflg = iapar2(1)
  Goto (200, 200, 300) iflg
  Print *, 'IAPAR2(1) must be 1, 2, or 3'
  Stop
  200 Return
  300 Call arini1
  Call artord
  Return
End Subroutine arini
