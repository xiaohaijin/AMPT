Subroutine zpca1a(i)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
  Common /ana1/ts(12)
  Save
  If (iconfg==1) Then
     t1 = fts(i)
     t2 = ft(i)
     ipic = 11
  Else If (iconfg==2 .Or. iconfg==3) Then
     t1 = taus(i)
     t2 = tau(i)
     ipic = 12
  Else If (iconfg==4 .Or. iconfg==5) Then
     t1 = fts(i)
     t2 = ft(i)
     ipic = 12
  End If
  If (iconfg<=3) Then
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           rapi = raps(i)
           et = dsqrt(pxs(i)**2+pys(i)**2+xmp**2)
           Call zpca1b(rapi, et, ian)
        End If
     End Do
  Else
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           p0 = es(i)
           p1 = pxs(i)
           p2 = pys(i)
           p3 = pzs(i)
           Call zpca1c(p0, p1, p2, p3, ian)
        End If
     End Do
  End If
  Return
End Subroutine zpca1a
