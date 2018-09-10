Subroutine ftime
  Implicit Double Precision (A-H, O-Z)
  External ftime1
  Parameter (maxptn=400001)
  Common /para1/mul
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /par1/formt
  Common /anim/nevent, isoft, isflag, izpc
  Common /rndm3/iseedp
  Save
  iseed = iseedp
  Do i = 1, maxptn
     ct(i) = 0D0
     ot(i) = 0D0
  End Do
  tlarge = 1000000.D0
  If (iftflg==0) Then
     If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        Do i = 1, mul
           If (ft0(i)>tlarge) ft0(i) = tlarge
        End Do
        Goto 150
     Else
        Do i = 1, maxptn
           ft0(i) = tlarge
        End Do
        Do i = 1, mul
           xmt2 = px0(i)**2 + py0(i)**2 + xmp**2
           formt = xmt2/e0(i)
           ft0(i) = ftime1(iseed)
           If (ft0(i)>tlarge) ft0(i) = tlarge
        End Do
     End If
  End If
150 Continue
  If (mul>1) Then
     Call index1(maxptn, mul, ft0, indx)
  Else
     indx(1) = 1
  End If
  Return
End Subroutine ftime
