Real Function fmassn(dmass)
  Save
  am0 = 1.44
  fmassn = am0*w1440(dmass)/((dmass**2-am0**2)**2+am0**2*w1440(dmass)**2)
  Return
End Function fmassn
