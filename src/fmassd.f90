Real Function fmassd(dmass)
  Save
  am0 = 1.232
  fmassd = am0*width(dmass)/((dmass**2-am0**2)**2+am0**2*width(dmass)**2)
  Return
End Function fmassd
