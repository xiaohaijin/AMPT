Block Data ppbdat
  Parameter (amp=0.93828, amn=0.939457, am0=1.232, am1440=1.44, am1535=1.535)
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  Data thresh/1.87656, 1.877737, 1.878914, 2.17028, 2.171457, 2.37828, 2.379457, 2.464, 2.47328, 2.474457, 2.672, 2.767, 2.88, 2.975, 3.07/
  Data (ppbm(i,1), i=1, 15)/amp, amp, amn, amp, amn, amp, amn, am0, amp, amn, am0, am0, am1440, am1440, am1535/
  Data (ppbm(i,2), i=1, 15)/amp, amn, amn, am0, am0, am1440, am1440, am0, am1535, am1535, am1440, am1535, am1440, am1535, am1535/
  Data factr2/0, 1, 1.17E-01, 3.27E-03, 3.58E-05, 1.93E-07/
  Data niso/1, 2, 1, 16, 16, 4, 4, 64, 4, 4, 32, 32, 4, 8, 4/
End Block Data ppbdat
