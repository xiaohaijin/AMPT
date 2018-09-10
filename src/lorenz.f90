Subroutine lorenz(energy, px, py, pz, bex, bey, bez)
  Implicit Double Precision (A-H, O-Z)
  Common /lor/enenew, pxnew, pynew, pznew
  Save
  beta2 = bex**2 + bey**2 + bez**2
  If (beta2==0D0) Then
     enenew = energy
     pxnew = px
     pynew = py
     pznew = pz
  Else
     If (beta2>0.999999999999999D0) Then
        beta2 = 0.999999999999999D0
        Print *, 'beta2=0.999999999999999'
     End If
     gam = 1.D0/dsqrt(1.D0-beta2)
     enenew = gam*(energy-bex*px-bey*py-bez*pz)
     pxnew = -gam*bex*energy + (1.D0+(gam-1.D0)*bex**2/beta2)*px + (gam-1.D0)*bex*bey/beta2*py + (gam-1.D0)*bex*bez/beta2*pz
     pynew = -gam*bey*energy + (gam-1.D0)*bex*bey/beta2*px + (1.D0+(gam-1.D0)*bey**2/beta2)*py + (gam-1.D0)*bey*bez/beta2*pz
     pznew = -gam*bez*energy + (gam-1.D0)*bex*bez/beta2*px + (gam-1.D0)*bey*bez/beta2*py + (1.D0+(gam-1.D0)*bez**2/beta2)*pz
  End If
  Return
End Subroutine lorenz
