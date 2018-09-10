Function romg(x)
  Dimension fr(0:1000)
  Save
  Data i0/0/
  If (i0/=0) Goto 100
  Do i = 1, 1001
    xr = (i-1)*0.01
    fr(i-1) = omg0(xr)
  End Do
  100 i0 = 1
  If (x>=10.0) Then
    romg = 0.0
    Return
  End If
  ix = int(x*100)
  romg = (fr(ix)*((ix+1)*0.01-x)+fr(ix+1)*(x-ix*0.01))/0.01
  Return
End Function romg
