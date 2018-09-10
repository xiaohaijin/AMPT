Subroutine dchcel(l, i, j, k, t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist2/icell, icel(10, 10, 10)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  If (i==11 .Or. j==11 .Or. k==11) Then
     If (.Not. (i==11 .And. j==11 .And. k==11)) Stop 'cerr'
     m = icell
  Else
     m = icel(i, j, k)
  End If
  If (m==0) Return
  If (next(m)==l) Then
     tm = tlarge
     last0 = 0
     Call reor(t, tm, m, last0)
  End If
  n = nic(m)
  If (n==0) Return
  Do While (n/=m)
     If (next(n)==l) Then
        tm = tlarge
        last0 = 0
        Call reor(t, tm, n, last0)
     End If
     n = nic(n)
  End Do
  Return
End Subroutine dchcel
