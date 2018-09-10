Subroutine exchge(isg, ipi, jsg, ipj)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxstr=150001)
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  Save
  k1 = k1sgs(isg, ipi)
  k2 = k2sgs(isg, ipi)
  px = pxsgs(isg, ipi)
  py = pysgs(isg, ipi)
  pz = pzsgs(isg, ipi)
  pe = pesgs(isg, ipi)
  pm = pmsgs(isg, ipi)
  gx = gxsgs(isg, ipi)
  gy = gysgs(isg, ipi)
  gz = gzsgs(isg, ipi)
  ft = ftsgs(isg, ipi)
  k1sgs(isg, ipi) = k1sgs(jsg, ipj)
  k2sgs(isg, ipi) = k2sgs(jsg, ipj)
  pxsgs(isg, ipi) = pxsgs(jsg, ipj)
  pysgs(isg, ipi) = pysgs(jsg, ipj)
  pzsgs(isg, ipi) = pzsgs(jsg, ipj)
  pesgs(isg, ipi) = pesgs(jsg, ipj)
  pmsgs(isg, ipi) = pmsgs(jsg, ipj)
  gxsgs(isg, ipi) = gxsgs(jsg, ipj)
  gysgs(isg, ipi) = gysgs(jsg, ipj)
  gzsgs(isg, ipi) = gzsgs(jsg, ipj)
  ftsgs(isg, ipi) = ftsgs(jsg, ipj)
  k1sgs(jsg, ipj) = k1
  k2sgs(jsg, ipj) = k2
  pxsgs(jsg, ipj) = px
  pysgs(jsg, ipj) = py
  pzsgs(jsg, ipj) = pz
  pesgs(jsg, ipj) = pe
  pmsgs(jsg, ipj) = pm
  gxsgs(jsg, ipj) = gx
  gysgs(jsg, ipj) = gy
  gzsgs(jsg, ipj) = gz
  ftsgs(jsg, ipj) = ft
  Return
End Subroutine exchge
