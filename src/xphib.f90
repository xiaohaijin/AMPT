Subroutine xphib(lb1, lb2, em1, em2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, sigp)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Parameter (aka=0.498, ala=1.1157, pimass=0.140, aphi=1.02)
  Parameter (arho=0.77)
  Save
  sigp = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  srrt = srt - (em1+em2)
  xsk1 = 8.00
  If (srt>(ap1+amn)) Then
    xsk2 = 0.0235*srrt**(-0.519)
  End If
  If (srt>(ap1+am0)) Then
    If (srrt<0.7) Then
      xsk3 = 0.0119*srrt**(-0.534)
    Else
      xsk3 = 0.0130*srrt**(-0.304)
    End If
  End If
  If (srt>(arho+amn)) Then
    If (srrt<0.7) Then
      xsk4 = 0.0166*srrt**(-0.786)
    Else
      xsk4 = 0.0189*srrt**(-0.277)
    End If
  End If
  If (srt>(arho+am0)) Then
    If (srrt<0.7) Then
      xsk5 = 0.0119*srrt**(-0.534)
    Else
      xsk5 = 0.0130*srrt**(-0.304)
    End If
  End If
  If ((lb1>=1 .And. lb1<=2) .Or. (lb2>=1 .And. lb2<=2)) Then
    If (srt>(aka+ala)) Then
      xsk6 = 1.715/((srrt+3.508)**2-12.138)
    End If
  End If
  sigp = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6
  Return
End Subroutine xphib
