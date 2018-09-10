Subroutine frztm(nevnt, idd)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Dimension tsf(31)
  Common /para1/mul
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  Common /frzout/xnprod(30), etprod(30), xnfrz(30), etfrz(30), dnprod(30), detpro(30), dnfrz(30), detfrz(30)
  Save
  Data tsf/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30/
  If (idd==0) Then
     Do ii = 1, 30
        xnprod(ii) = 0D0
        etprod(ii) = 0D0
        xnfrz(ii) = 0D0
        etfrz(ii) = 0D0
        dnprod(ii) = 0D0
        detpro(ii) = 0D0
        dnfrz(ii) = 0D0
        detfrz(ii) = 0D0
     End Do
     Open (86, File='ana1/production.dat', Status='UNKNOWN')
     Open (87, File='ana1/freezeout.dat', Status='UNKNOWN')
  Else If (idd==1) Then
     Do ip = 1, mul
        Do ii = 1, 30
           eth0 = dsqrt(px0(ip)**2+py0(ip)**2+xmass0(ip)**2)
           eth2 = dsqrt(px5(ip)**2+py5(ip)**2+xmass5(ip)**2)
           If (ft0(ip)<tsf(ii+1)) Then
              xnprod(ii) = xnprod(ii) + 1D0
              etprod(ii) = etprod(ii) + eth0
              If (ft0(ip)>=tsf(ii)) Then
                 dnprod(ii) = dnprod(ii) + 1D0
                 detpro(ii) = detpro(ii) + eth0
              End If
           End If
           If (ft5(ip)<tsf(ii+1)) Then
              xnfrz(ii) = xnfrz(ii) + 1D0
              etfrz(ii) = etfrz(ii) + eth2
              If (ft5(ip)>=tsf(ii)) Then
                 dnfrz(ii) = dnfrz(ii) + 1D0
                 detfrz(ii) = detfrz(ii) + eth2
              End If
           End If
        End Do
     End Do
  Else If (idd==2) Then
     Write (86, *) '       t,       np,       dnp/dt,      etp ' // ' detp/dt'
     Write (87, *) '       t,       nf,       dnf/dt,      etf ' // ' detf/dt'
     Do ii = 1, 30
        xnp = xnprod(ii)/dble(nevnt)
        xnf = xnfrz(ii)/dble(nevnt)
        etp = etprod(ii)/dble(nevnt)
        etf = etfrz(ii)/dble(nevnt)
        dxnp = dnprod(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        dxnf = dnfrz(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        detp = detpro(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        detf = detfrz(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        Write (86, 200) tsf(ii+1), xnp, dxnp, etp, detp
        Write (87, 200) tsf(ii+1), xnf, dxnf, etf, detf
     End Do
  End If
  Return
200 Format (2X, F9.2, 4(2X,F10.2))
End Subroutine frztm
