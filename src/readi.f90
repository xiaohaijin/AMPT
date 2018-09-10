Subroutine readi
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Double Precision field(9)
  Common /para1/mul
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Save
  Do i = 1, maxptn
     If (ievt/=1 .And. i==1) Then
        ityp0(i) = ntyp
        gx0(1) = field(1)
        gy0(1) = field(2)
        gz0(1) = field(3)
        ft0(1) = field(4)
        px0(1) = field(5)
        py0(1) = field(6)
        pz0(1) = field(7)
        e0(1) = field(8)
        xmass0(i) = field(9)
        mul = 1
     Else
900     Read (27, *, End=1000) neve, ntyp, field
        If (neve<nsevt) Goto 900
        If (neve>nsevt+ievt-1) Goto 1000
        ityp0(i) = ntyp
        gx0(i) = field(1)
        gy0(i) = field(2)
        gz0(i) = field(3)
        ft0(i) = field(4)
        px0(i) = field(5)
        py0(i) = field(6)
        pz0(i) = field(7)
        e0(i) = field(8)
        xmass0(i) = field(9)
        mul = mul + 1
     End If
  End Do
1000 Continue
  Return
End Subroutine readi
