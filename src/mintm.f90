Subroutine mintm(i, j, tmin, nc)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /aurec1/jxa, jya, jza
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  Logical allok
  Call isco(i, j, allok, tm, t1, t2)
  If (allok .And. tm<tmin) Then
     tmin = tm
     ct(i) = t1
     nc = j
     If (iconfg==3 .Or. iconfg==5) Then
        dgxa(i) = -jxa*10D0*size1
        dgya(i) = -jya*10D0*size2
        If (iconfg==5) Then
           dgza(i) = -jza*10D0*size3
        End If
     End If
  End If
  Return
End Subroutine mintm
