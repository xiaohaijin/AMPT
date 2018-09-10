Subroutine ud2(i, j, t, tmin, nc)
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
  If (allok) Then
     If (tm<tmin) Then
        tmin = tm
        ct(j) = t2
        nc = i
        If (iconfg==3 .Or. iconfg==5) Then
           dgxa(j) = jxa*10D0*size1
           dgya(j) = jya*10D0*size2
           If (iconfg==5) Then
              dgza(j) = jza*10D0*size3
           End If
        End If
     End If
     If (tm<=ot(i)) Then
        ct(i) = t1
        icels0 = icels(i)
        i1 = icels0/10000
        i2 = (icels0-i1*10000)/100
        i3 = icels0 - i1*10000 - i2*100
        Call wallc(i, i1, i2, i3, t, tmin1)
        Call fixtim(i, t, tmin1, tm, j)
        If (iconfg==3 .Or. iconfg==5) Then
           dgxa(i) = -jxa*10D0*size1
           dgya(i) = -jya*10D0*size2
           If (iconfg==5) Then
              dgza(i) = -jza*10D0*size3
           End If
        End If
     End If
     If (tm>ot(i) .And. next(i)==j) Then
        ct(i) = t1
        Call reor(t, tm, i, j)
     End If
  Else If (next(i)==j) Then
     tm = tlarge
     Call reor(t, tm, i, j)
  End If
  Return
End Subroutine ud2
