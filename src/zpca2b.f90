Subroutine zpca2b
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /ana1/ts(12)
  Save
  Do i = 1, ichkpt
     t1 = ft(i)
     t2 = tlarge
     ipic = 12
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           p0 = e(i)
           p1 = px(i)
           p2 = py(i)
           p3 = pz(i)
           Call zpca1c(p0, p1, p2, p3, ian)
        End If
     End Do
  End Do
  Return
End Subroutine zpca2b
