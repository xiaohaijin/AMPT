Subroutine pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Parameter (aka=0.498, ala=1.1157, pimass=0.140, aphi=1.02)
  Parameter (arho=0.77)
  Save
  xphi = 0.0
  xphin = 0.0
  xphid = 0.0
  If ((lb1>=3 .And. lb1<=5) .Or. (lb2>=3 .And. lb2<=5)) Then
    If ((iabs(lb1)>=1 .And. iabs(lb1)<=2) .Or. (iabs(lb2)>=1 .And. iabs(lb2)<=2)) Then
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        sig = 0.0235*srrt**(-0.519)
        xphin = sig*1.*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        sig = 0.0235*srrt**(-0.519)
        xphid = sig*4.*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    Else
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphin = sig*(1./4.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphid = sig*1.*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    End If
  Else
    If ((iabs(lb1)>=1 .And. iabs(lb1)<=2) .Or. (iabs(lb2)>=1 .And. iabs(lb2)<=2)) Then
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0166*srrt**(-0.786)
        Else
          sig = 0.0189*srrt**(-0.277)
        End If
        xphin = sig*(1./3.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0166*srrt**(-0.786)
        Else
          sig = 0.0189*srrt**(-0.277)
        End If
        xphid = sig*(4./3.)*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    Else
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphin = sig*(1./12.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphid = sig*(1./3.)*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    End If
  End If
  xphin = xphin/10.
  xphid = xphid/10.
  xphi = xphin + xphid
  Return
End Subroutine pibphi
