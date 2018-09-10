Subroutine pyxtot
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension bcs(5, 8), bcb(2, 5), bcc(3)
  Data ((bcs(i,j),j=1,8), i=1, 5)/41.74, 0.66, 0.0000, 337., 0.0, 0.0, -39.3, 0.48, 41.66, 0.60, 0.0000, 306., 0.0, 0.0, -34.6, 0.51, 41.36, 0.63, 0.0000, 299., 7.3, 0.5, -40.4, 0.47, 41.68, 0.63, 0.0083, 330., 0.0, 0.0, -39.0, 0.48, 41.13, 0.59, 0.0074, 278., 10.5, 0.5, -41.2, 0.46/
  Data ((bcb(i,j),j=1,5), i=1, 2)/10.79, -0.049, 0.040, 21.5, 1.23, 9.92, -0.027, 0.013, 18.9, 1.07/
  Data bcc/2.0164346, -0.5590311, 0.0376279/
  nfit = min(5, max(1,mstp(31)))
  sigp = bcs(nfit, 1) + bcs(nfit, 2)*(-0.25*paru(1)**2*(1.-0.25*bcs(nfit,3)*paru(1)**2)+(1.+0.5*bcs(nfit,3)*paru(1)**2)*(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)*(log(vint(2)/bcs(nfit,4)))**4)/((1.-0.25*bcs(nfit,3)*paru(1)**2)**2+2.*bcs(nfit,3)*(1.+0.25*bcs(nfit,3)*paru(1)**2)*(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)**2*(log(vint(2)/bcs(nfit,4)))**4) + bcs(nfit, 5)*vint(2)**(bcs(nfit,6)-1.)*sin(0.5*paru(1)*bcs(nfit,6))
  sigm = -bcs(nfit, 7)*vint(2)**(bcs(nfit,8)-1.)*cos(0.5*paru(1)*bcs(nfit,8))
  refp = bcs(nfit, 2)*paru(1)*log(vint(2)/bcs(nfit,4))/((1.-0.25*bcs(nfit,3)*paru(1)**2)**2+2.*bcs(nfit,3)*(1.+0.25*bcs(nfit,3)*paru(1)**2)+(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)**2*(log(vint(2)/bcs(nfit,4)))**4) - bcs(nfit, 5)*vint(2)**(bcs(nfit,6)-1.)*cos(0.5*paru(1)*bcs(nfit,6))
  refm = -bcs(nfit, 7)*vint(2)**(bcs(nfit,8)-1.)*sin(0.5*paru(1)*bcs(nfit,8))
  sigma = sigp - isign(1, mint(11)*mint(12))*sigm
  rho = (refp-isign(1,mint(11)*mint(12))*refm)/sigma
  nfit = 1
  If (mstp(31)>=4) nfit = 2
  bp = bcb(nfit, 1) + bcb(nfit, 2)*log(vint(2)) + bcb(nfit, 3)*(log(vint(2)))**2
  bm = bcb(nfit, 4) + bcb(nfit, 5)*log(vint(2))
  b = bp - isign(1, mint(11)*mint(12))*sigm/sigp*(bm-bp)
  vint(121) = b
  c = -0.5*bcc(2)/bcc(3)*(1.-sqrt(max(0.,1.+4.*bcc(3)/bcc(2)**2*(1.E-03*vint(1)-bcc(1)))))
  vint(122) = c
  sigel = sigma**2*(1.+rho**2)/(16.*paru(1)*paru(5)*b)
  sigsd = 2.*0.68*(1.+36./vint(2))*log(0.6+0.1*vint(2))
  sigdd = sigsd**2/(3.*sigel)
  signd = sigma - sigdd - sigsd - sigel
  If (iabs(mint(11))==211 .And. iabs(mint(12))==211) Then
    sigma = 4./9.*sigma
    sigdd = 4./9.*sigdd
    sigsd = 4./9.*sigsd
    sigel = 4./9.*sigel
    signd = 4./9.*signd
  Else If (iabs(mint(11))==211 .Or. iabs(mint(12))==211) Then
    sigma = 2./3.*sigma
    sigdd = 2./3.*sigdd
    sigsd = 2./3.*sigsd
    sigel = 2./3.*sigel
    signd = 2./3.*signd
  End If
  vint(101) = sigma
  vint(102) = sigel
  vint(103) = sigsd
  vint(104) = sigdd
  vint(106) = signd
  xsec(95, 1) = signd
  Return
End Subroutine pyxtot
