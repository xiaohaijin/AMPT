Subroutine coulin(masspr, massta, num)
  Integer zta, zpr
  Parameter (maxstr=150001)
  Common /ee/id(maxstr), lb(maxstr)
  Common /zz/zta, zpr
  Save
  mass = massta + masspr
  Do irun = 1, num
     Do i = 1 + (irun-1)*mass, zta + (irun-1)*mass
        lb(i) = 1
     End Do
     Do i = zta + 1 + (irun-1)*mass, massta + (irun-1)*mass
        lb(i) = 2
     End Do
     Do i = massta + 1 + (irun-1)*mass, massta + zpr + (irun-1)*mass
        lb(i) = 1
     End Do
     Do i = massta + zpr + 1 + (irun-1)*mass, massta + masspr + (irun-1)*mass
        lb(i) = 2
     End Do
  End Do
  Return
End Subroutine coulin
