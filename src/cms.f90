Subroutine cms(i1, i2, px1cm, py1cm, pz1cm, srt)
  Parameter (maxstr=150001)
  Double Precision px1, py1, pz1, px2, py2, pz2, em1, em2, e1, e2, s, etotal, p1beta, transf, dbetax, dbetay, dbetaz, dgamma, scheck
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /bg/betax, betay, betaz, gamma
  Save
  px1 = dble(p(1,i1))
  py1 = dble(p(2,i1))
  pz1 = dble(p(3,i1))
  px2 = dble(p(1,i2))
  py2 = dble(p(2,i2))
  pz2 = dble(p(3,i2))
  em1 = dble(e(i1))
  em2 = dble(e(i2))
  e1 = dsqrt(em1**2+px1**2+py1**2+pz1**2)
  e2 = dsqrt(em2**2+px2**2+py2**2+pz2**2)
  s = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
  If (s<=0) s = 0D0
  srt = sngl(dsqrt(s))
  etotal = e1 + e2
  dbetax = (px1+px2)/etotal
  dbetay = (py1+py2)/etotal
  dbetaz = (pz1+pz2)/etotal
  scheck = 1.D0 - dbetax**2 - dbetay**2 - dbetaz**2
  If (scheck<=0D0) Then
     Write (99, *) 'scheck1: ', scheck
     Stop
  End If
  dgamma = 1.D0/dsqrt(scheck)
  p1beta = px1*dbetax + py1*dbetay + pz1*dbetaz
  transf = dgamma*(dgamma*p1beta/(dgamma+1D0)-e1)
  px1cm = sngl(dbetax*transf+px1)
  py1cm = sngl(dbetay*transf+py1)
  pz1cm = sngl(dbetaz*transf+pz1)
  betax = sngl(dbetax)
  betay = sngl(dbetay)
  betaz = sngl(dbetaz)
  gamma = sngl(dgamma)
  Return
End Subroutine cms
