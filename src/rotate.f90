Subroutine rotate(px0, py0, pz0, px, py, pz)
  Save
  pr0 = sqrt(px0**2+py0**2+pz0**2)
  If (pr0==0) pr0 = 0.00000001
  c2 = pz0/pr0
  If (px0==0.0 .And. py0==0.0) Then
    t2 = 0.0
  Else
    t2 = atan2(py0, px0)
  End If
  scheck = 1.0 - c2**2
  If (scheck<0) Then
    Write (99, *) 'scheck45: ', scheck
    scheck = 0.
  End If
  s2 = sqrt(scheck)
  ct2 = cos(t2)
  st2 = sin(t2)
  pr = sqrt(px**2+py**2+pz**2)
  If (pr==0) pr = 0.0000001
  c1 = pz/pr
  If (px==0 .And. py==0) Then
    t1 = 0.
  Else
    t1 = atan2(py, px)
  End If
  scheck = 1.0 - c1**2
  If (scheck<0) Then
    Write (99, *) 'scheck46: ', scheck
    scheck = 0.
  End If
  s1 = sqrt(scheck)
  ct1 = cos(t1)
  st1 = sin(t1)
  ss = c2*s1*ct1 + s2*c1
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  pz = pr*(c1*c2-s1*s2*ct1)
  Return
End Subroutine rotate
