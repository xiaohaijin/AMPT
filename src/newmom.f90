Subroutine newmom(t)
  Implicit Double Precision (A-H, O-Z)
  Parameter (hbarc=0.197327054D0)
  Parameter (maxptn=400001)
  Parameter (pi=3.14159265358979D0)
  Common /para1/mul
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para5/iconfg, iordsc
  Common /para6/centy
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /aurec1/jxa, jya, jza
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  Common /lor/enenew, pxnew, pynew, pznew
  Common /cprod/xn1, xn2, xn3
  Common /rndm2/iff
  Common /anim/nevent, isoft, isflag, izpc
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  Save
  If (isoft==5) Then
     If (ifrz(iscat)==1 .Or. ifrz(jscat)==1) Then
        last(iscat) = jscat
        last(jscat) = iscat
        Return
     End If
  End If
  iff = -iff
  If (iconfg==2 .Or. iconfg==4) Then
     icels1 = icels(iscat)
     i1 = icels1/10000
     j1 = (icels1-i1*10000)/100
     icels2 = icels(jscat)
     i2 = icels2/10000
     j2 = (icels2-i2*10000)/100
     If (iconfg==4) Then
        k1 = icels1 - i1*10000 - j1*100
        k2 = icels2 - i2*10000 - j2*100
     End If
  End If
  px1 = px(iscat)
  py1 = py(iscat)
  pz1 = pz(iscat)
  e1 = e(iscat)
  x1 = gx(iscat)
  y1 = gy(iscat)
  z1 = gz(iscat)
  t1 = ft(iscat)
  px2 = px(jscat)
  py2 = py(jscat)
  pz2 = pz(jscat)
  e2 = e(jscat)
  If (iconfg==1) Then
     x2 = gx(jscat)
     y2 = gy(jscat)
     z2 = gz(jscat)
  Else If (iconfg==2 .Or. iconfg==4) Then
     If (i1-i2>5) Then
        x2 = gx(jscat) + 10D0*size1
     Else If (i1-i2<-5) Then
        x2 = gx(jscat) - 10D0*size1
     Else
        x2 = gx(jscat)
     End If
     If (j1-j2>5) Then
        y2 = gy(jscat) + 10D0*size2
     Else If (j1-j2<-5) Then
        y2 = gy(jscat) - 10D0*size2
     Else
        y2 = gy(jscat)
     End If
     If (iconfg==4) Then
        If (k1-k2>5) Then
           z2 = gz(jscat) + 10D0*size3
        Else If (k1-k2<-5) Then
           z2 = gz(jscat) - 10D0*size3
        Else
           z2 = gz(jscat)
        End If
     Else
        z2 = gz(jscat)
     End If
  Else If (iconfg==3 .Or. iconfg==5) Then
     x2 = gx(jscat) + dgxa(jscat)
     y2 = gy(jscat) + dgya(jscat)
     If (iconfg==5) Then
        z2 = gz(jscat) + dgza(jscat)
     Else
        z2 = gz(jscat)
     End If
  End If
  t2 = ft(jscat)
  rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
  bex = (px1+px2)/(e1+e2)
  bey = (py1+py2)/(e1+e2)
  bez = (pz1+pz2)/(e1+e2)
  Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
  px1 = pxnew
  py1 = pynew
  pz1 = pznew
  e1 = enenew
  pp2 = pxnew**2 + pynew**2 + pznew**2
  Call getht(iscat, jscat, pp2, that)
  theta = dacos(that/(2D0*pp2)+1D0)
  theta = dble(iff)*theta
  Call lorenz(t1, x1, y1, z1, bex, bey, bez)
  x1 = pxnew
  y1 = pynew
  z1 = pznew
  Call lorenz(t2, x2, y2, z2, bex, bey, bez)
  x2 = pxnew
  y2 = pynew
  z2 = pznew
  Call cropro(x1-x2, y1-y2, z1-z2, px1, py1, pz1)
  Call xnormv(xn1, xn2, xn3)
  Call zprota(xn1, xn2, xn3, theta, px1, py1, pz1)
  px2 = -px1
  py2 = -py1
  pz2 = -pz1
  e2 = dsqrt(px2**2+py2**2+pz2**2+xmass(jscat)**2)
  Call lorenz(e1, px1, py1, pz1, -bex, -bey, -bez)
  px(iscat) = pxnew
  py(iscat) = pynew
  pz(iscat) = pznew
  e(iscat) = enenew
  Call lorenz(e2, px2, py2, pz2, -bex, -bey, -bez)
  px(jscat) = pxnew
  py(jscat) = pynew
  pz(jscat) = pznew
  e(jscat) = enenew
  vx(iscat) = px(iscat)/e(iscat)
  vy(iscat) = py(iscat)/e(iscat)
  vz(iscat) = pz(iscat)/e(iscat)
  vx(jscat) = px(jscat)/e(jscat)
  vy(jscat) = py(jscat)/e(jscat)
  vz(jscat) = pz(jscat)/e(jscat)
  last(iscat) = jscat
  last(jscat) = iscat
  If (iconfg<=3) Then
     If (e(iscat)<=abs(pz(iscat))) Then
        rap(iscat) = 1000000.D0
     Else
        rap(iscat) = 0.5D0*log((e(iscat)+pz(iscat))/(e(iscat)-pz(iscat)))
     End If
     If (e(jscat)<=abs(pz(jscat))) Then
        rap(jscat) = 1000000.D0
     Else
        rap(jscat) = 0.5D0*log((e(jscat)+pz(jscat))/(e(jscat)-pz(jscat)))
     End If
     rap1 = rap(iscat)
     rap2 = rap(jscat)
     If ((rap1<centy+0.5D0 .And. rap1>centy-0.5D0)) Then
     End If
     If ((rap2<centy+0.5D0 .And. rap2>centy-0.5D0)) Then
     End If
  End If
  Return
End Subroutine newmom
