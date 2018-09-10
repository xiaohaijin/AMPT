Function hirnd2(i, xmin, xmax)
  Common /hijhb/rr(10, 201), xx(10, 201)
  Common /rndf77/nseed
  Save
  If (xmin<xx(i,1)) xmin = xx(i, 1)
  If (xmax>xx(i,201)) xmax = xx(i, 201)
  jmin = 1 + int(200*(xmin-xx(i,1))/(xx(i,201)-xx(i,1)))
  jmax = 1 + int(200*(xmax-xx(i,1))/(xx(i,201)-xx(i,1)))
  rx = rr(i, jmin) + (rr(i,jmax)-rr(i,jmin))*ranart(nseed)
  jl = 0
  ju = 202
  10 If (ju-jl>1) Then
    jm = (ju+jl)/2
    If ((rr(i,201)>rr(i,1)) .Eqv. (rx>rr(i,jm))) Then
      jl = jm
    Else
      ju = jm
    End If
    Goto 10
  End If
  j = jl
  If (j<1) j = 1
  If (j>=201) j = 200
  hirnd2 = (xx(i,j)+xx(i,j+1))/2.0
  Return
End Function hirnd2
