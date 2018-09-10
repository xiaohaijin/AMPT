Subroutine newpos(t, i)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Save
  dt1 = ct(i) - ft(i)
  gx(i) = gx(i) + vx(i)*dt1
  gy(i) = gy(i) + vy(i)*dt1
  gz(i) = gz(i) + vz(i)*dt1
  ft(i) = ct(i)
  If (iconfg<=3) Then
     If (ft(i)<=abs(gz(i))) Then
        eta(i) = 1000000.D0
     Else
        eta(i) = 0.5D0*log((ft(i)+gz(i))/(ft(i)-gz(i)))
     End If
     If (eta(i)<1000000.D0) Then
        tau(i) = ft(i)/cosh(eta(i))
     Else
        tau(i) = 1D-10
     End If
  End If
  Return
End Subroutine newpos
