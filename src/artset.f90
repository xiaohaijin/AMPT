Subroutine artset
  Parameter (amu=0.9383, nxymax=10001)
  Double Precision dpcoal, drcoal, ecritl
  Integer zta, zpr
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  Common /zz/zta, zpr
  Common /run/num
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /input3/plab, elab, zeropt, b0, bi, bm, dencut, cycbox
  Common /imulst/iperts
  Common /coal/dpcoal, drcoal, ecritl
  Common /anim/nevent, isoft, isflag, izpc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Common /xyembed/nxyjet, xyjet(nxymax, 2)
  Save
  ecritl = 1.D0
  massta = 1
  zta = 1
  masspr = 1
  zpr = 1
  plab = 14.6
  iplab = 2
  If (iplab==2) Then
    elab = sqrt(plab**2+amu**2) - amu
  Else
    elab = plab
  End If
  elab = elab*1000.
  zeropt = 0.
  iseed = 700721
  iperts = 0
  manyb = 1
  b0 = 1
  bi = 0
  bm = 0
  icoll = -1
  num = 1
  insys = 1
  ipot = 3
  mode = 0
  If (icoll==-1) ipot = 0
  dx = 2.73
  dy = 2.73
  dz = 2.73
  dpx = 0.6
  dpy = 0.6
  dpz = 0.6
  iavoid = 1
  imomen = 1
  If (icoll==-1) imomen = 3
  nfreq = 10
  icflow = 0
  icrho = 0
  icou = 0
  kpoten = 0
  kmul = 1
  dencut = 15
  cycbox = 0
  If (ioscar==2 .Or. ioscar==3) Then
    Open (92, File='../data/parton-initial-afterPropagation.dat', Status='UNKNOWN')
  End If
  If (ioscar==3) Then
    Open (95, File='../data/parton-collisionsHistory.dat', Status='UNKNOWN')
    Open (96, File='../data/minijet-initial-beforePropagation.dat', Status='UNKNOWN')
    If (isoft==4 .Or. isoft==5) Then
      Open (85, File='../data/parton-after-coalescence.dat', Status='UNKNOWN')
    End If
  End If
  Open (94, File='../data/npart-xy.dat', Status='UNKNOWN')
  If (iembed==3 .Or. iembed==4) Then
    Open (97, File='embed-jet-xy.txt', Status='UNKNOWN')
    Read (97, *) nxyjet
    If (nevent>nxyjet) Then
      If (nxyjet>nxymax) Then
        Print *, 'Too many lines in embed-jet-xy.txt:             increase value of the parameter nxymax'
        Stop
      Else If (nxyjet<=0) Then
        Print *, 'Check number of entries in embed-jet-xy.txt'
        Stop
      End If
      Do ixy = 1, nxyjet
        Read (97, *) xyjet(ixy, 1), xyjet(ixy, 2)
      End Do
    End If
  End If
  Return
End Subroutine artset
