Subroutine genei
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /para5/iconfg, iordsc
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /lor/enenew, pxnew, pynew, pznew
  Common /rndm3/iseedp
  Save
  External ran1
  iseed = iseedp
  incmul = 4000
  temp = 0.5D0
  etamin = -5D0
  etamax = 5D0
  r0 = 5D0
  tau0 = 0.1D0
  deta = etamax - etamin
  Do i = mul + 1, mul + incmul
     ityp0(i) = 21
     xmass0(i) = xmp
     Call energy(e, temp)
     Call momntm(px, py, pz, e)
     e = dsqrt(e**2+xmp**2)
     If (iconfg<=3) Then
        eta(i) = etamin + deta*ran1(iseed)
        bex = 0D0
        bey = 0D0
        bez = -tanh(eta(i))
        Call lorenz(e, px, py, pz, bex, bey, bez)
        px0(i) = pxnew
        py0(i) = pynew
        pz0(i) = pznew
        e0(i) = enenew
     Else
        px0(i) = px
        py0(i) = py
        pz0(i) = pz
        e0(i) = e
     End If
  End Do
  Do i = mul + 1, mul + incmul
     If (iconfg<=3) Then
        gz0(i) = tau0*sinh(eta(i))
        ft0(i) = tau0*cosh(eta(i))
        If (iconfg==1) Then
           Call posit1(x, y, r0)
           gx0(i) = x + px0(i)*ft0(i)/e0(i)
           gy0(i) = y + py0(i)*ft0(i)/e0(i)
        Else If (iconfg==2 .Or. iconfg==3) Then
           Call posit2(x, y)
           gx0(i) = x
           gy0(i) = y
        End If
     Else
        ft0(i) = 0D0
        Call posit3(x, y, z)
        gx0(i) = x
        gy0(i) = y
        gz0(i) = z
     End If
  End Do
  mul = mul + incmul
  If (mul>=maxptn .Or. mul==0) Then
     Print *, 'event', ievt, 'has', mul, 'number of gluon', 'adjusting counting is necessary'
     Stop 'adarr'
  End If
  Return
End Subroutine genei
