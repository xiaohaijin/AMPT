Subroutine vegas(fxn, avgi, sd, chi2a)
  Implicit Double Precision (A-H, O-Z)
  Common /bveg1/xl(10), xu(10), acc, ndim, ncall, itmx, nprn
  Common /bveg2/xi(50, 10), si, si2, swgt, schi, ndo, it
  Common /bveg3/f, ti, tsi
  External fxn
  Dimension d(50, 10), di(50, 10), xin(50), r(50), dx(10), dt(10), x(10), kg(10), ia(10)
  Real qran(10)
  Save
  Data ndmx/50/, alph/1.5D0/, one/1.D0/, mds/ -1/
  ndo = 1
  Do j = 1, ndim
    xi(1, j) = one
  End Do
  Entry vegas1(fxn, avgi, sd, chi2a)
  it = 0
  si = 0.D0
  si2 = si
  swgt = si
  schi = si
  Entry vegas2(fxn, avgi, sd, chi2a)
  nd = ndmx
  ng = 1
  If (mds==0) Goto 2
  ng = int((real(ncall)/2.)**(1./real(ndim)))
  mds = 1
  If ((2*ng-ndmx)<0) Goto 2
  mds = -1
  npg = ng/ndmx + 1
  nd = ng/npg
  ng = npg*nd
  2 k = ng**ndim
  npg = ncall/k
  If (npg<2) npg = 2
  calls = npg*k
  dxg = one/ng
  dv2g = (calls*dxg**ndim)**2/npg/npg/(npg-one)
  xnd = nd
  ndm = nd - 1
  dxg = dxg*xnd
  xjac = one/calls
  Do j = 1, ndim
    dx(j) = xu(j) - xl(j)
    xjac = xjac*dx(j)
  End Do
  If (nd==ndo) Goto 8
  rc = ndo/xnd
  Do j = 1, ndim
    k = 0
    xn = 0.D0
    dr = xn
    i = k
    4 k = k + 1
    dr = dr + one
    xo = xn
    xn = xi(k, j)
    5 If (rc>dr) Goto 4
    i = i + 1
    dr = dr - rc
    xin(i) = xn - (xn-xo)*dr
    If (i<ndm) Goto 5
    Do i = 1, ndm
      xi(i, j) = xin(i)
    End Do
    xi(nd, j) = one
  End Do
  ndo = nd
  8 Continue
  Entry vegas3(fxn, avgi, sd, chi2a)
  9 it = it + 1
  ti = 0.D0
  tsi = ti
  Do j = 1, ndim
    kg(j) = 1
    Do i = 1, nd
      d(i, j) = ti
      di(i, j) = ti
    End Do
  End Do
  11 fb = 0.D0
  f2b = fb
  k = 0
  12 k = k + 1
  Call aran9(qran, ndim)
  wgt = xjac
  Do j = 1, ndim
    xn = dble(float(kg(j))-qran(j))*dxg + one
    ia(j) = int(xn)
    If (ia(j)>1) Goto 13
    xo = xi(ia(j), j)
    rc = (xn-ia(j))*xo
    Goto 14
    13 xo = xi(ia(j), j) - xi(ia(j)-1, j)
    rc = xi(ia(j)-1, j) + (xn-ia(j))*xo
    14 x(j) = xl(j) + rc*dx(j)
    wgt = wgt*xo*xnd
  End Do
  f = wgt
  f = f*fxn(x, wgt)
  f2 = f*f
  fb = fb + f
  f2b = f2b + f2
  Do j = 1, ndim
    di(ia(j), j) = di(ia(j), j) + f
    If (mds>=0) d(ia(j), j) = d(ia(j), j) + f2
  End Do
  If (k<npg) Goto 12
  f2b = dsqrt(f2b*npg)
  f2b = (f2b-fb)*(f2b+fb)
  ti = ti + fb
  tsi = tsi + f2b
  If (mds>=0) Goto 18
  Do j = 1, ndim
    d(ia(j), j) = d(ia(j), j) + f2b
  End Do
  18 k = ndim
  19 kg(k) = mod(kg(k), ng) + 1
  If (kg(k)/=1) Goto 11
  k = k - 1
  If (k>0) Goto 19
  tsi = tsi*dv2g
  ti2 = ti*ti
  wgt = ti2/(tsi+1.0D-37)
  si = si + ti*wgt
  si2 = si2 + ti2
  swgt = swgt + wgt
  swgt = swgt + 1.0D-37
  si2 = si2 + 1.0D-37
  schi = schi + ti2*wgt
  avgi = si/swgt
  sd = swgt*it/si2
  chi2a = sd*(schi/swgt-avgi*avgi)/dble(float(it)-.999)
  sd = dsqrt(one/sd)
  If (nprn==0) Goto 21
  tsi = dsqrt(tsi)
  21 Do j = 1, ndim
    xo = d(1, j)
    xn = d(2, j)
    d(1, j) = (xo+xn)/2.D0
    dt(j) = d(1, j)
    Do i = 2, ndm
      d(i, j) = xo + xn
      xo = xn
      xn = d(i+1, j)
      d(i, j) = (d(i,j)+xn)/3.D0
      dt(j) = dt(j) + d(i, j)
    End Do
    d(nd, j) = (xn+xo)/2.D0
    dt(j) = dt(j) + d(nd, j)
  End Do
  Do j = 1, ndim
    rc = 0.D0
    Do i = 1, nd
      r(i) = 0.D0
      If (dt(j)>=1.0D18) Then
        Write (6, *) '************** A SINGULARITY >1.0D18'
      End If
      If (d(i,j)<=1.0D-18) Goto 24
      xo = dt(j)/d(i, j)
      r(i) = ((xo-one)/xo/dlog(xo))**alph
      24 rc = rc + r(i)
    End Do
    rc = rc/xnd
    k = 0
    xn = 0.D0
    dr = xn
    i = k
    25 k = k + 1
    dr = dr + r(k)
    xo = xn
    xn = xi(k, j)
    26 If (rc>dr) Goto 25
    i = i + 1
    dr = dr - rc
    xin(i) = xn - (xn-xo)*dr/(r(k)+1.0D-30)
    If (i<ndm) Goto 26
    Do i = 1, ndm
      xi(i, j) = xin(i)
    End Do
    xi(nd, j) = one
  End Do
  If (it<itmx .And. acc*dabs(avgi)<sd) Goto 9
  Return
End Subroutine vegas
