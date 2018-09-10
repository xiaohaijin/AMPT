Subroutine decomp(px0, py0, pz0, xm0, i, itq1)
  Implicit Double Precision (D)
  Double Precision enenew, pxnew, pynew, pznew
  Double Precision de0, beta2, gam, ptwo, px0, py0, pz0, xm0
  Common /lor/enenew, pxnew, pynew, pznew
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /decom/ptwo(2, 5)
  Common /rndf77/nseed
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Save
  dcth = dble(ranart(nseed))*2.D0 - 1.D0
  dphi = dble(ranart(nseed)*hipr1(40))*2.D0
  If (iembed>=1 .And. iembed<=4) Then
     If (i==(natt-2*nsembd) .Or. i==(natt-2*nsembd-1)) Then
        dcth = 0.D0
        dphi = dble(phidecomp)
     End If
  End If
  ds = xm0**2
  dpcm = dsqrt((ds-(ptwo(1,5)+ptwo(2,5))**2)*(ds-(ptwo(1,5)-ptwo(2,5))**2)/ds/4D0)
  dpz = dpcm*dcth
  dpx = dpcm*dsqrt(1.D0-dcth**2)*dcos(dphi)
  dpy = dpcm*dsqrt(1.D0-dcth**2)*dsin(dphi)
  de1 = dsqrt(ptwo(1,5)**2+dpcm**2)
  de2 = dsqrt(ptwo(2,5)**2+dpcm**2)
  de0 = dsqrt(px0**2+py0**2+pz0**2+xm0**2)
  dbex = px0/de0
  dbey = py0/de0
  dbez = pz0/de0
  beta2 = dbex**2 + dbey**2 + dbez**2
  gam = 1.D0/dsqrt(1.D0-beta2)
  If (beta2>=0.9999999999999D0) Then
     Write (6, *) '1', dbex, dbey, dbez, beta2, gam
  End If
  Call lorenz(de1, dpx, dpy, dpz, -dbex, -dbey, -dbez)
  ptwo(1, 1) = pxnew
  ptwo(1, 2) = pynew
  ptwo(1, 3) = pznew
  ptwo(1, 4) = enenew
  Call lorenz(de2, -dpx, -dpy, -dpz, -dbex, -dbey, -dbez)
  ptwo(2, 1) = pxnew
  ptwo(2, 2) = pynew
  ptwo(2, 3) = pznew
  ptwo(2, 4) = enenew
  Return
End Subroutine decomp
