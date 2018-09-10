Subroutine lorntz(ilo, b, pi, pj)
  Dimension pi(4), pj(4), b(3)
  Save
  bb = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
  deno3 = sqrt(1.-bb)
  If (deno3==0.) deno3 = 1.E-10
  gam = 1./deno3
  ga = gam*gam/(gam+1.)
  If (ilo==1) Goto 100
  pib = pi(1)*b(1) + pi(2)*b(2) + pi(3)*b(3)
  pjb = pj(1)*b(1) + pj(2)*b(2) + pj(3)*b(3)
  Do i = 1, 3
    pi(i) = pi(i) + b(i)*(ga*pib-gam*pi(4))
    pj(i) = pj(i) + b(i)*(ga*pjb-gam*pj(4))
  End Do
  pi(4) = gam*(pi(4)-pib)
  pj(4) = gam*(pj(4)-pjb)
  Return
  100 Continue
  pib = pi(1)*b(1) + pi(2)*b(2) + pi(3)*b(3)
  pjb = pj(1)*b(1) + pj(2)*b(2) + pj(3)*b(3)
  Do i = 1, 3
    pi(i) = pi(i) + b(i)*(ga*pib+gam*pi(4))
    pj(i) = pj(i) + b(i)*(ga*pjb+gam*pj(4))
  End Do
  pi(4) = gam*(pi(4)+pib)
  pj(4) = gam*(pj(4)+pjb)
  Return
End Subroutine lorntz
