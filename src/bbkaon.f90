Subroutine bbkaon(ic, srt, px, py, pz, ana, plx, ply, plz, ala, pkx, pky, pkz, icou1)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save
  pi = 3.1415962
  icou1 = 0
  aka = 0.498
  ala = 1.116
  If (ic==2 .Or. ic==4) ala = 1.197
  ana = 0.939
  If (ic>2) Then
    dmax = srt - aka - ala - 0.02
    dm1 = rmass(dmax, iseed)
    ana = dm1
  End If
  t1 = aka + ana + ala
  t2 = ana + ala - aka
  If (srt<=t1) Then
    icou1 = -1
    Return
  End If
  pmax = sqrt((srt**2-t1**2)*(srt**2-t2**2))/(2.*srt)
  If (pmax==0.) pmax = 1.E-09
  ntry = 0
  1 pk = pmax*ranart(nseed)
  ntry = ntry + 1
  prob = fkaon(pk, pmax)
  If ((prob<ranart(nseed)) .And. (ntry<=40)) Goto 1
  cs = 1. - 2.*ranart(nseed)
  ss = sqrt(1.-cs**2)
  fai = 2.*3.14*ranart(nseed)
  pkx = pk*ss*cos(fai)
  pky = pk*ss*sin(fai)
  pkz = pk*cs
  ek = sqrt(aka**2+pk**2)
  eln = srt - ek
  If (eln<=0) Then
    icou1 = -1
    Return
  End If
  bx = -pkx/eln
  by = -pky/eln
  bz = -pkz/eln
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
    Write (99, *) 'scheck44: ', scheck
    Stop
  End If
  ga = 1./sqrt(scheck)
  elnc = eln/ga
  pn2 = ((elnc**2+ana**2-ala**2)/(2.*elnc))**2 - ana**2
  If (pn2<=0.) pn2 = 1.E-09
  pn = sqrt(pn2)
  csn = 1. - 2.*ranart(nseed)
  ssn = sqrt(1.-csn**2)
  fain = 2.*3.14*ranart(nseed)
  px = pn*ssn*cos(fain)
  py = pn*ssn*sin(fain)
  pz = pn*csn
  en = sqrt(ana**2+pn2)
  plx = -px
  ply = -py
  plz = -pz
  pbeta = px*bx + py*by + pz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  px = bx*trans0 + px
  py = by*trans0 + py
  pz = bz*trans0 + pz
  el = sqrt(ala**2+plx**2+ply**2+plz**2)
  pbeta = plx*bx + ply*by + plz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+el)
  plx = bx*trans0 + plx
  ply = by*trans0 + ply
  plz = bz*trans0 + plz
  Return
End Subroutine bbkaon
