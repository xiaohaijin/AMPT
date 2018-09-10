Real Function dirct3(srt)
  Dimension arrayj(17), arrayl(17), arraym(17), arrayw(17), arrayb(17)
  Save
  Data arrayj/1.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/2, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.52, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.125, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.55, 0.6, 0.375, 0.6, 0.1, 0.15, 0.15, 0.05, 0.35, 0.3, 0.15, 0.1, 0.1, 0.22, 0.2, 0.09, 0.4/
  pi = 3.1415926
  amn = 0.938
  amp = 0.138
  xs = 0
  branch = 1./3.
  Do ir = 1, 17
    If (ir>8) branch = 2./3.
    xs0 = fd1(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), srt)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  dirct3 = xs
  Return
End Function dirct3
