Subroutine zpca2c
  Implicit Double Precision (A-H, O-Z)
  Character *8 code, versn
  Character *4 reffra
  Integer aproj, zproj, atarg, ztarg, event
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Save
  Data nff/0/
  If (nff==0) Then
     Write (26, 101) 'OSCAR1997A'
     Write (26, 101) 'final_id_p_x'
     code = 'ZPC'
     versn = '1.0.1'
     aproj = -1
     zproj = -1
     atarg = -1
     ztarg = -1
     reffra = 'cm'
     ebeam = 0D0
     ntestp = 1
     Write (26, 102) code, versn, aproj, zproj, atarg, ztarg, reffra, ebeam, ntestp
     nff = 1
     event = 1
     bimp = 0D0
     phi = 0D0
  End If
  Write (26, 103) event, mul, bimp, phi
  Do i = 1, mul
     Write (26, 104) i, ityp(i), px(i), py(i), pz(i), e(i), xmass(i), gx(i), gy(i), gz(i), ft(i)
  End Do
  event = event + 1
  Return
101 Format (A12)
102 Format (2(A8,2X), '(', I3, ',', I6, ')+(', I3, ',', I6, ')', 2X, A4, 2X, E10.4, 2X, I8)
103 Format (I10, 2X, I10, 2X, F8.3, 2X, F8.3)
104 Format (I10, 2X, I10, 2X, 9(E12.6,2X))
End Subroutine zpca2c
