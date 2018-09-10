Subroutine pyspli(kf, kflin, kflch, kflsp)
  Dimension kfl(3)
  kfa = iabs(kf)
  kfs = isign(1, kf)
  kfl(1) = mod(kfa/1000, 10)
  kfl(2) = mod(kfa/100, 10)
  kfl(3) = mod(kfa/10, 10)
  kflr = kflin*kfs
  kflch = 0
  If (kfl(1)==0) Then
    kfl(2) = kfl(2)*(-1)**kfl(2)
    kfl(3) = -kfl(3)*(-1)**iabs(kfl(2))
    If (kflr==kfl(2)) Then
      kflsp = kfl(3)
    Else If (kflr==kfl(3)) Then
      kflsp = kfl(2)
    Else If (iabs(kflr)==21 .And. rlu(0)>0.5) Then
      kflsp = kfl(2)
      kflch = kfl(3)
    Else If (iabs(kflr)==21) Then
      kflsp = kfl(3)
      kflch = kfl(2)
    Else If (kflr*kfl(2)>0) Then
      Call lukfdi(-kflr, kfl(2), kfdump, kflch)
      kflsp = kfl(3)
    Else
      Call lukfdi(-kflr, kfl(3), kfdump, kflch)
      kflsp = kfl(2)
    End If
  Else
    nagr = 0
    Do j = 1, 3
      If (kflr==kfl(j)) nagr = nagr + 1
    End Do
    If (nagr>=1) Then
      ragr = 0.00001 + (nagr-0.00002)*rlu(0)
      iagr = 0
      Do j = 1, 3
        If (kflr==kfl(j)) ragr = ragr - 1.
        If (iagr==0 .And. ragr<=0.) iagr = j
      End Do
    Else
      iagr = int(1.00001+2.99998*rlu(0))
    End If
    id1 = 1
    If (iagr==1) id1 = 2
    If (iagr==1 .And. kfl(3)>kfl(2)) id1 = 3
    id2 = 6 - iagr - id1
    ksp = 3
    If (mod(kfa,10)==2 .And. kfl(1)==kfl(2)) Then
      If (iagr/=3 .And. rlu(0)>0.25) ksp = 1
    Else If (mod(kfa,10)==2 .And. kfl(2)>=kfl(3)) Then
      If (iagr/=1 .And. rlu(0)>0.25) ksp = 1
    Else If (mod(kfa,10)==2) Then
      If (iagr==1) ksp = 1
      If (iagr/=1 .And. rlu(0)>0.75) ksp = 1
    End If
    kflsp = 1000*kfl(id1) + 100*kfl(id2) + ksp
    If (kflin==21) Then
      kflch = kfl(iagr)
    Else If (nagr==0 .And. kflr>0) Then
      Call lukfdi(-kflr, kfl(iagr), kfdump, kflch)
    Else If (nagr==0) Then
      Call lukfdi(10000+kflsp, -kflr, kfdump, kflch)
      kflsp = kfl(iagr)
    End If
  End If
  kflch = kflch*kfs
  kflsp = kflsp*kfs
  Return
End Subroutine pyspli
