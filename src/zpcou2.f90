Subroutine zpcou2
  Implicit Double Precision (A-H, O-Z)
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ana1/ts(12)
  Common /ana3/em(4, 4, 12)
  Save
  Open (28, File='ana4/em.dat', Status='unknown')
  vol = 1000.D0*size1*size2*size3
  ntotal = nevnt*nsbrun
  Do ian = 1, 12
     Write (28, *) '*** for time ', ts(ian), 'fm(s)'
     Do i = 1, 4
        Write (28, *) em(i, 1, ian)/vol/ntotal, em(i, 2, ian)/vol/ntotal, em(i, 3, ian)/vol/ntotal, em(i, 4, ian)/vol/ntotal
     End Do
  End Do
  Return
End Subroutine zpcou2
