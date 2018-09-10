Function ptr(ptmax, iseed)
  Common /table/xarray(0:1000), earray(0:1000)
  Common /rndf77/nseed
  Save
  ptr = 0.
  If (ptmax<=1.E-02) Then
    ptr = ptmax
    Return
  End If
  If (ptmax>2.01) ptmax = 2.01
  tryial = ptdis(ptmax)/ptdis(2.01)
  xt = ranart(nseed)*tryial
  Do ie = 1, 200
    If (earray(ie)==xt) Then
      ptr = xarray(ie)
      Return
    End If
    If (xarray(ie-1)<=0.00001) Goto 50
    If (xarray(ie)<=0.00001) Goto 50
    If (earray(ie-1)<=0.00001) Goto 50
    If (earray(ie)<=0.00001) Goto 50
    If (earray(ie)>xt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      ptr = exp(ymin+(alog(xt)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (ptr>ptmax) ptr = ptmax
      Return
    End If
  50 End Do
  Return
End Function ptr
