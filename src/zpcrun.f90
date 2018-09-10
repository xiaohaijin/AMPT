Subroutine zpcrun(*)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (tend1=250D0)
  Parameter (tend2=6.1D0)
  Common /para1/mul
  Common /para5/iconfg, iordsc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  Common /ilist6/t, iopern, icolln
  Common /ana1/ts(12)
  Common /anim/nevent, isoft, isflag, izpc
  Common /arevt/iaevt, iarun, miss
  Save
  If (mod(ictype,2)==0) Then
     Call savrec(iscat)
     Call savrec(jscat)
  End If
  Call getict(t1)
  If (iconfg==1 .And. t1>tlarge/2D0) Return 1
  If (iconfg==2 .Or. iconfg==3) Then
     If (t1>tend1) Return 1
  End If
  If (iconfg==4 .Or. iconfg==5) Then
     If (t1>tend2) Return 1
  End If
  If (isoft==5) Then
     Call local(t1)
  End If
  iopern = iopern + 1
  t = t1
  If (mod(ictype,2)==0) Then
     icolln = icolln + 1
  End If
  If (iconfg==1 .Or. iconfg==2 .Or. iconfg==4) Then
     If (ictype==1 .Or. ictype==2 .Or. ictype==5 .Or. ictype==6) Then
        Call celasn
     End If
  End If
  If (ictype/=1) Then
     iscat0 = iscat
     jscat0 = jscat
     iscat = max0(iscat0, jscat0)
     jscat = min0(iscat0, jscat0)
     If (jscat/=0) Then
        If (next(jscat)/=iscat) Then
           Print *, 'iscat=', iscat, 'jscat=', jscat, 'next(', jscat, ')=', next(jscat)
           If (ct(iscat)<tlarge/2D0) Stop 'tterr'
           If (ct(jscat)<tlarge/2D0) Stop 'tterr'
        End If
     End If
     niscat = iscat
     njscat = jscat
     If (icsta(iscat)/=0) Call cellre(niscat, t)
     If (jscat/=0) Then
        If (icsta(jscat)/=0) Call cellre(njscat, t)
     End If
     If (mod(ictype,2)==0) Then
        If (ioscar==3) Then
           Write (95, *) 'event,miss,iscat,jscat=', iaevt, miss, iscat, jscat
           If (dmax1(abs(gx(iscat)),abs(gy(iscat)),abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))<9999) Then
              Write (95, 200) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 200) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           Else
              Write (95, 201) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 201) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           End If
        End If
        Call scat(t, iscat, jscat)
        If (ioscar==3) Then
           If (dmax1(abs(gx(iscat)),abs(gy(iscat)),abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))<9999) Then
              Write (95, 200) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 200) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           Else
              Write (95, 201) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 201) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           End If
        End If
     End If
  End If
  Call ulist(t)
  If (ifmpt<=mul) Then
     If (ictype/=0 .And. ictype/=3 .And. ictype/=4) Then
        ichkpt = ichkpt + 1
        ifmpt = ifmpt + 1
     End If
  End If
  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
End Subroutine zpcrun
