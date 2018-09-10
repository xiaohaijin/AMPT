Subroutine locldr(icall, drlocl)
  Implicit Double Precision (A-H, O-Z)
  Dimension ftp0(3), pxp0(3), pyp0(3), pzp0(3), pep0(3)
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  Common /prtn23/gxp0(3), gyp0(3), gzp0(3), ft0fom
  Common /lor/enenew, pxnew, pynew, pznew
  Save
  If (icall==2) Then
     etot = pep(1) + pep(2)
     bex = (pxp(1)+pxp(2))/etot
     bey = (pyp(1)+pyp(2))/etot
     bez = (pzp(1)+pzp(2))/etot
     Do j = 1, 2
        beta2 = bex**2 + bey**2 + bez**2
        gam = 1.D0/dsqrt(1.D0-beta2)
        If (beta2>=0.9999999999999D0) Then
           Write (6, *) '4', pxp(1), pxp(2), pyp(1), pyp(2), pzp(1), pzp(2), pep(1), pep(2), pmp(1), pmp(2), dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2+pmp(1)**2)/pep(1), dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2)/pep(1)
           Write (6, *) '4a', pxp(1) + pxp(2), pyp(1) + pyp(2), pzp(1) + pzp(2), etot
           Write (6, *) '4b', bex, bey, bez, beta2, gam
        End If
        Call lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez)
        gxp0(j) = pxnew
        gyp0(j) = pynew
        gzp0(j) = pznew
        ftp0(j) = enenew
        Call lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez)
        pxp0(j) = pxnew
        pyp0(j) = pynew
        pzp0(j) = pznew
        pep0(j) = enenew
     End Do
     If (ftp0(1)>=ftp0(2)) Then
        ilate = 1
        iearly = 2
     Else
        ilate = 2
        iearly = 1
     End If
     ft0fom = ftp0(ilate)
     dt0 = ftp0(ilate) - ftp0(iearly)
     gxp0(iearly) = gxp0(iearly) + pxp0(iearly)/pep0(iearly)*dt0
     gyp0(iearly) = gyp0(iearly) + pyp0(iearly)/pep0(iearly)*dt0
     gzp0(iearly) = gzp0(iearly) + pzp0(iearly)/pep0(iearly)*dt0
     drlocl = dsqrt((gxp0(ilate)-gxp0(iearly))**2+(gyp0(ilate)-gyp0(iearly))**2+(gzp0(ilate)-gzp0(iearly))**2)
  Else If (icall==3) Then
     etot = pep(1) + pep(2) + pep(3)
     bex = (pxp(1)+pxp(2)+pxp(3))/etot
     bey = (pyp(1)+pyp(2)+pyp(3))/etot
     bez = (pzp(1)+pzp(2)+pzp(3))/etot
     beta2 = bex**2 + bey**2 + bez**2
     gam = 1.D0/dsqrt(1.D0-beta2)
     If (beta2>=0.9999999999999D0) Then
        Write (6, *) '5', bex, bey, bez, beta2, gam
     End If
     Do j = 1, 3
        Call lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez)
        gxp0(j) = pxnew
        gyp0(j) = pynew
        gzp0(j) = pznew
        ftp0(j) = enenew
        Call lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez)
        pxp0(j) = pxnew
        pyp0(j) = pynew
        pzp0(j) = pznew
        pep0(j) = enenew
     End Do
     If (ftp0(1)>ftp0(2)) Then
        ilate = 1
        If (ftp0(3)>ftp0(1)) ilate = 3
     Else
        ilate = 2
        If (ftp0(3)>=ftp0(2)) ilate = 3
     End If
     ft0fom = ftp0(ilate)
     If (ilate==1) Then
        imin = 2
        imax = 3
        istep = 1
     Else If (ilate==2) Then
        imin = 1
        imax = 3
        istep = 2
     Else If (ilate==3) Then
        imin = 1
        imax = 2
        istep = 1
     End If
     Do iearly = imin, imax, istep
        dt0 = ftp0(ilate) - ftp0(iearly)
        gxp0(iearly) = gxp0(iearly) + pxp0(iearly)/pep0(iearly)*dt0
        gyp0(iearly) = gyp0(iearly) + pyp0(iearly)/pep0(iearly)*dt0
        gzp0(iearly) = gzp0(iearly) + pzp0(iearly)/pep0(iearly)*dt0
     End Do
  End If
  Return
End Subroutine locldr
