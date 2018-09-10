Subroutine hoscar
  Parameter (maxstr=150001, amn=0.939457, amp=0.93828)
  Character *8 code, reffra, frame
  Character *25 amptvn
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Common /lastt/itimeh, bimp
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  Common /oscar1/iap, izp, iat, izt
  Common /oscar2/frame, amptvn
  Save
  Data nff/0/
  If (nff==0) Then
     Write (19, 101) 'OSCAR1997A'
     Write (19, 111) 'final_id_p_x'
     code = 'AMPT'
     If (frame=='CMS') Then
        reffra = 'nncm'
        xmp = (amp*izp+amn*(iap-izp))/iap
        xmt = (amp*izt+amn*(iat-izt))/iat
        ebeam = (efrm**2-xmp**2-xmt**2)/2./xmt
     Else If (frame=='LAB') Then
        reffra = 'lab'
        ebeam = efrm
     Else
        reffra = 'unknown'
        ebeam = 0.
     End If
     ntestp = 1
     Write (19, 102) code, amptvn, iap, izp, iat, izt, reffra, ebeam, ntestp
     nff = 1
     ievent = 1
     phi = 0.
     If (frame=='CMS') Write (19, 112) efrm
  End If
  Write (19, 103) ievent, nlast, bimp, phi
  Do i = 1, nlast
     ene = sqrt(plast(1,i)**2+plast(2,i)**2+plast(3,i)**2+plast(4,i)**2)
     Write (19, 104) i, invflv(lblast(i)), plast(1, i), plast(2, i), plast(3, i), ene, plast(4, i), xlast(1, i), xlast(2, i), xlast(3, i), xlast(4, i)
  End Do
  ievent = ievent + 1
  Return
101 Format (A10)
111 Format (A12)
102 Format (A4, 1X, A20, 1X, '(', I3, ',', I3, ')+(', I3, ',', I3, ')', 2X, A4, 2X, E10.4, 2X, I8)
103 Format (I10, 2X, I10, 2X, F8.3, 2X, F8.3)
104 Format (I10, 2X, I10, 2X, 9(E12.6,2X))
112 Format ('# Center-of-mass energy/nucleon-pair is', F12.3, 'GeV')
End Subroutine hoscar
