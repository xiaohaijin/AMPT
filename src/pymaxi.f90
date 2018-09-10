Subroutine pymaxi
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Common /pyint6/proc(0:200)
  Character proc*28
  Save /pyint6/
  Character cvar(4)*4
  Dimension npts(4), mvarpt(200, 4), vintpt(200, 30), sigspt(200), narel(6), wtrel(6), wtmat(6, 6), coefu(6), iaccmx(4), sigsmx(4), sigssm(3)
  Data cvar/'tau ', 'tau''', 'y*  ', 'cth '/
  vint(143) = 1.
  vint(144) = 1.
  xsec(0, 1) = 0.
  Do isub = 1, 200
    If (isub>=91 .And. isub<=95) Then
      xsec(isub, 1) = vint(isub+11)
      If (msub(isub)/=1) Goto 350
      Goto 340
    Else If (isub==96) Then
      If (mint(43)/=4) Goto 350
      If (msub(95)/=1 .And. mstp(81)<=0 .And. mstp(131)<=0) Goto 350
    Else If (isub==11 .Or. isub==12 .Or. isub==13 .Or. isub==28 .Or. isub==53 .Or. isub==68) Then
      If (msub(isub)/=1 .Or. msub(95)==1) Goto 350
    Else
      If (msub(isub)/=1) Goto 350
    End If
    mint(1) = isub
    istsb = iset(isub)
    If (isub==96) istsb = 2
    If (mstp(122)>=2) Write (mstu(11), 1000) isub
    mint(72) = 0
    kfr1 = 0
    If (istsb==1 .Or. istsb==3) Then
      kfr1 = kfpr(isub, 1)
    Else If (isub>=71 .And. isub<=77) Then
      kfr1 = 25
    End If
    If (kfr1/=0) Then
      taur1 = pmas(kfr1, 1)**2/vint(2)
      gamr1 = pmas(kfr1, 1)*pmas(kfr1, 2)/vint(2)
      mint(72) = 1
      mint(73) = kfr1
      vint(73) = taur1
      vint(74) = gamr1
    End If
    If (isub==141) Then
      kfr2 = 23
      taur2 = pmas(kfr2, 1)**2/vint(2)
      gamr2 = pmas(kfr2, 1)*pmas(kfr2, 2)/vint(2)
      mint(72) = 2
      mint(74) = kfr2
      vint(75) = taur2
      vint(76) = gamr2
    End If
    sqm3 = 0.
    sqm4 = 0.
    mint(71) = 0
    vint(71) = ckin(3)
    If (istsb==2 .Or. istsb==4) Then
      If (kfpr(isub,1)/=0) sqm3 = pmas(kfpr(isub,1), 1)**2
      If (kfpr(isub,2)/=0) sqm4 = pmas(kfpr(isub,2), 1)**2
      If (min(sqm3,sqm4)<ckin(6)**2) mint(71) = 1
      If (mint(71)==1) vint(71) = max(ckin(3), ckin(5))
      If (isub==96 .And. mstp(82)<=1) vint(71) = parp(81)
      If (isub==96 .And. mstp(82)>=2) vint(71) = 0.08*parp(82)
    End If
    vint(63) = sqm3
    vint(64) = sqm4
    npts(1) = 2 + 2*mint(72)
    If (mint(43)==1 .And. (istsb==1 .Or. istsb==2)) npts(1) = 1
    npts(2) = 1
    If (mint(43)>=2 .And. (istsb==3 .Or. istsb==4)) npts(2) = 2
    npts(3) = 1
    If (mint(43)==4) npts(3) = 3
    npts(4) = 1
    If (istsb==2 .Or. istsb==4) npts(4) = 5
    ntry = npts(1)*npts(2)*npts(3)*npts(4)
    Do j = 1, 20
      coef(isub, j) = 0.
    End Do
    coef(isub, 1) = 1.
    coef(isub, 7) = 0.5
    coef(isub, 8) = 0.5
    coef(isub, 10) = 1.
    coef(isub, 15) = 1.
    mcth = 0
    mtaup = 0
    cth = 0.
    taup = 0.
    sigsam = 0.
    Call pyklim(1)
    nacc = 0
    Do itry = 1, ntry
      If (mod(itry-1,npts(2)*npts(3)*npts(4))==0) Then
        mtau = 1 + (itry-1)/(npts(2)*npts(3)*npts(4))
        Call pykmap(1, mtau, 0.5)
        If (istsb==3 .Or. istsb==4) Call pyklim(4)
      End If
      If ((istsb==3 .Or. istsb==4) .And. mod(itry-1,npts(3)*npts(4))==0) Then
        mtaup = 1 + mod((itry-1)/(npts(3)*npts(4)), npts(2))
        Call pykmap(4, mtaup, 0.5)
      End If
      If (mod(itry-1,npts(3)*npts(4))==0) Call pyklim(2)
      If (mod(itry-1,npts(4))==0) Then
        myst = 1 + mod((itry-1)/npts(4), npts(3))
        Call pykmap(2, myst, 0.5)
        Call pyklim(3)
      End If
      If (istsb==2 .Or. istsb==4) Then
        mcth = 1 + mod(itry-1, npts(4))
        Call pykmap(3, mcth, 0.5)
      End If
      If (isub==96) vint(25) = vint(21)*(1.-vint(23)**2)
      mint(51) = 0
      Call pyklim(0)
      If (mint(51)==1) Goto 120
      nacc = nacc + 1
      mvarpt(nacc, 1) = mtau
      mvarpt(nacc, 2) = mtaup
      mvarpt(nacc, 3) = myst
      mvarpt(nacc, 4) = mcth
      Do j = 1, 30
        vintpt(nacc, j) = vint(10+j)
      End Do
      Call pysigh(nchn, sigs)
      sigspt(nacc) = sigs
      If (sigs>sigsam) sigsam = sigs
      If (mstp(122)>=2) Write (mstu(11), 1100) mtau, mtaup, myst, mcth, vint(21), vint(22), vint(23), vint(26), sigs
    120 End Do
    If (sigsam==0.) Then
      Write (mstu(11), 1200) isub
      Stop
    End If
    taumin = vint(11)
    taumax = vint(31)
    atau1 = log(taumax/taumin)
    atau2 = (taumax-taumin)/(taumax*taumin)
    If (npts(1)>=3) Then
      atau3 = log(taumax/taumin*(taumin+taur1)/(taumax+taur1))/taur1
      atau4 = (atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1))/gamr1
    End If
    If (npts(1)>=5) Then
      atau5 = log(taumax/taumin*(taumin+taur2)/(taumax+taur2))/taur2
      atau6 = (atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2))/gamr2
    End If
    ystmin = 0.5*log(taumin)
    ystmax = -ystmin
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))
    Do ivar = 1, 4
      If (npts(ivar)==1) Goto 230
      If (isub==96 .And. ivar==4) Goto 230
      nbin = npts(ivar)
      Do j1 = 1, nbin
        narel(j1) = 0
        wtrel(j1) = 0.
        coefu(j1) = 0.
        Do j2 = 1, nbin
          wtmat(j1, j2) = 0.
        End Do
      End Do
      Do iacc = 1, nacc
        ibin = mvarpt(iacc, ivar)
        narel(ibin) = narel(ibin) + 1
        wtrel(ibin) = wtrel(ibin) + sigspt(iacc)
        If (ivar==1) Then
          tau = vintpt(iacc, 11)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (atau1/atau2)/tau
          If (nbin>=3) Then
            wtmat(ibin, 3) = wtmat(ibin, 3) + (atau1/atau3)/(tau+taur1)
            wtmat(ibin, 4) = wtmat(ibin, 4) + (atau1/atau4)*tau/((tau-taur1)**2+gamr1**2)
          End If
          If (nbin>=5) Then
            wtmat(ibin, 5) = wtmat(ibin, 5) + (atau1/atau5)/(tau+taur2)
            wtmat(ibin, 6) = wtmat(ibin, 6) + (atau1/atau6)*tau/((tau-taur2)**2+gamr2**2)
          End If
        Else If (ivar==2) Then
          tau = vintpt(iacc, 11)
          taup = vintpt(iacc, 16)
          taupmn = vintpt(iacc, 6)
          taupmx = vintpt(iacc, 26)
          ataup1 = log(taupmx/taupmn)
          ataup2 = ((1.-tau/taupmx)**4-(1.-tau/taupmn)**4)/(4.*tau)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (ataup1/ataup2)*(1.-tau/taup)**3/taup
        Else If (ivar==3) Then
          yst = vintpt(iacc, 12)
          wtmat(ibin, 1) = wtmat(ibin, 1) + (ayst0/ayst1)*(yst-ystmin)
          wtmat(ibin, 2) = wtmat(ibin, 2) + (ayst0/ayst1)*(ystmax-yst)
          wtmat(ibin, 3) = wtmat(ibin, 3) + (ayst0/ayst3)/cosh(yst)
        Else
          rm34 = 2.*sqm3*sqm4/(vintpt(iacc,11)*vint(2))**2
          rsqm = 1. + rm34
          cthmax = sqrt(1.-4.*vint(71)**2/(taumax*vint(2)))
          cthmin = -cthmax
          If (cthmax>0.9999) rm34 = max(rm34, 2.*vint(71)**2/(taumax*vint(2)))
          acth1 = cthmax - cthmin
          acth2 = log(max(rm34,rsqm-cthmin)/max(rm34,rsqm-cthmax))
          acth3 = log(max(rm34,rsqm+cthmax)/max(rm34,rsqm+cthmin))
          acth4 = 1./max(rm34, rsqm-cthmax) - 1./max(rm34, rsqm-cthmin)
          acth5 = 1./max(rm34, rsqm+cthmin) - 1./max(rm34, rsqm+cthmax)
          cth = vintpt(iacc, 13)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (acth1/acth2)/max(rm34, rsqm-cth)
          wtmat(ibin, 3) = wtmat(ibin, 3) + (acth1/acth3)/max(rm34, rsqm+cth)
          wtmat(ibin, 4) = wtmat(ibin, 4) + (acth1/acth4)/max(rm34, rsqm-cth)**2
          wtmat(ibin, 5) = wtmat(ibin, 5) + (acth1/acth5)/max(rm34, rsqm+cth)**2
        End If
      End Do
      If (mstp(122)>=2) Write (mstu(11), 1300) cvar(ivar)
      msolv = 1
      Do ibin = 1, nbin
        If (mstp(122)>=2) Write (mstu(11), 1400)(wtmat(ibin,ired), ired=1, nbin), wtrel(ibin)
        If (narel(ibin)==0) msolv = 0
      End Do
      If (msolv==0) Then
        Do ibin = 1, nbin
          coefu(ibin) = 1.
        End Do
      Else
        Do ired = 1, nbin - 1
          Do ibin = ired + 1, nbin
            rqt = wtmat(ibin, ired)/wtmat(ired, ired)
            wtrel(ibin) = wtrel(ibin) - rqt*wtrel(ired)
            Do icoe = ired, nbin
              wtmat(ibin, icoe) = wtmat(ibin, icoe) - rqt*wtmat(ired, icoe)
            End Do
          End Do
        End Do
        Do ired = nbin, 1, -1
          Do icoe = ired + 1, nbin
            wtrel(ired) = wtrel(ired) - wtmat(ired, icoe)*coefu(icoe)
          End Do
          coefu(ired) = wtrel(ired)/wtmat(ired, ired)
        End Do
      End If
      coefsu = 0.
      Do ibin = 1, nbin
        coefu(ibin) = max(0., coefu(ibin))
        coefsu = coefsu + coefu(ibin)
      End Do
      If (ivar==1) ioff = 0
      If (ivar==2) ioff = 14
      If (ivar==3) ioff = 6
      If (ivar==4) ioff = 9
      If (coefsu>0.) Then
        Do ibin = 1, nbin
          coef(isub, ioff+ibin) = parp(121)/nbin + (1.-parp(121))*coefu(ibin)/coefsu
        End Do
      Else
        Do ibin = 1, nbin
          coef(isub, ioff+ibin) = 1./nbin
        End Do
      End If
      If (mstp(122)>=2) Write (mstu(11), 1500) cvar(ivar), (coef(isub,ioff+ibin), ibin=1, nbin)
    230 End Do
    Do j = 1, 4
      iaccmx(j) = 0
      sigsmx(j) = 0.
    End Do
    nmax = 0
    Do iacc = 1, nacc
      Do j = 1, 30
        vint(10+j) = vintpt(iacc, j)
      End Do
      Call pysigh(nchn, sigs)
      ieq = 0
      Do imv = 1, nmax
        If (abs(sigs-sigsmx(imv))<1E-4*(sigs+sigsmx(imv))) ieq = imv
      End Do
      If (ieq==0) Then
        Do imv = nmax, 1, -1
          iin = imv + 1
          If (sigs<=sigsmx(imv)) Goto 280
          iaccmx(imv+1) = iaccmx(imv)
          sigsmx(imv+1) = sigsmx(imv)
        End Do
        iin = 1
        280 iaccmx(iin) = iacc
        sigsmx(iin) = sigs
        If (nmax<=1) nmax = nmax + 1
      End If
    End Do
    If (mstp(122)>=2) Write (mstu(11), 1600)
    sigsam = sigsmx(1)
    Do imax = 1, nmax
      iacc = iaccmx(imax)
      mtau = mvarpt(iacc, 1)
      mtaup = mvarpt(iacc, 2)
      myst = mvarpt(iacc, 3)
      mcth = mvarpt(iacc, 4)
      vtau = 0.5
      vyst = 0.5
      vcth = 0.5
      vtaup = 0.5
      Do irpt = 1, 2
        Do ivar = 1, 4
          If (npts(ivar)==1) Goto 310
          If (ivar==1) vvar = vtau
          If (ivar==2) vvar = vtaup
          If (ivar==3) vvar = vyst
          If (ivar==4) vvar = vcth
          If (ivar==1) mvar = mtau
          If (ivar==2) mvar = mtaup
          If (ivar==3) mvar = myst
          If (ivar==4) mvar = mcth
          If (irpt==1) vdel = 0.1
          If (irpt==2) vdel = max(0.01, min(0.05,vvar-0.02,0.98-vvar))
          If (irpt==1) vmar = 0.02
          If (irpt==2) vmar = 0.002
          imov0 = 1
          If (irpt==1 .And. ivar==1) imov0 = 0
          Do imov = imov0, 8
            If (imov==0) Then
              inew = 2
              vnew = vvar
            Else If (imov==1) Then
              inew = 3
              vnew = vvar + vdel
            Else If (imov==2) Then
              inew = 1
              vnew = vvar - vdel
            Else If (sigssm(3)>=max(sigssm(1),sigssm(2)) .And. vvar+2.*vdel<1.-vmar) Then
              vvar = vvar + vdel
              sigssm(1) = sigssm(2)
              sigssm(2) = sigssm(3)
              inew = 3
              vnew = vvar + vdel
            Else If (sigssm(1)>=max(sigssm(2),sigssm(3)) .And. vvar-2.*vdel>vmar) Then
              vvar = vvar - vdel
              sigssm(3) = sigssm(2)
              sigssm(2) = sigssm(1)
              inew = 1
              vnew = vvar - vdel
            Else If (sigssm(3)>=sigssm(1)) Then
              vdel = 0.5*vdel
              vvar = vvar + vdel
              sigssm(1) = sigssm(2)
              inew = 2
              vnew = vvar
            Else
              vdel = 0.5*vdel
              vvar = vvar - vdel
              sigssm(3) = sigssm(2)
              inew = 2
              vnew = vvar
            End If
            If (ivar==1) Then
              vtau = vnew
              Call pykmap(1, mtau, vtau)
              If (istsb==3 .Or. istsb==4) Call pyklim(4)
            End If
            If (ivar<=2 .And. (istsb==3 .Or. istsb==4)) Then
              If (ivar==2) vtaup = vnew
              Call pykmap(4, mtaup, vtaup)
            End If
            If (ivar<=2) Call pyklim(2)
            If (ivar<=3) Then
              If (ivar==3) vyst = vnew
              Call pykmap(2, myst, vyst)
              Call pyklim(3)
            End If
            If (istsb==2 .Or. istsb==4) Then
              If (ivar==4) vcth = vnew
              Call pykmap(3, mcth, vcth)
            End If
            If (isub==96) vint(25) = vint(21)*(1.-vint(23)**2)
            Call pysigh(nchn, sigs)
            sigssm(inew) = sigs
            If (sigs>sigsam) sigsam = sigs
            If (mstp(122)>=2) Write (mstu(11), 1700) imax, ivar, mvar, imov, vnew, vint(21), vint(22), vint(23), vint(26), sigs
          End Do
        310 End Do
      End Do
      If (imax==1) sigs11 = sigsam
    End Do
    xsec(isub, 1) = 1.05*sigsam
    340 If (isub/=96) xsec(0, 1) = xsec(0, 1) + xsec(isub, 1)
  350 End Do
  If (mstp(122)>=1) Then
    Write (mstu(11), 1800)
    Write (mstu(11), 1900)
    Do isub = 1, 200
      If (msub(isub)/=1 .And. isub/=96) Goto 360
      If (isub==96 .And. mint(43)/=4) Goto 360
      If (isub==96 .And. msub(95)/=1 .And. mstp(81)<=0) Goto 360
      If (msub(95)==1 .And. (isub==11 .Or. isub==12 .Or. isub==13 .Or. isub==28 .Or. isub==53 .Or. isub==68)) Goto 360
      Write (mstu(11), 2000) isub, proc(isub), xsec(isub, 1)
    360 End Do
    Write (mstu(11), 2100)
  End If
  Return
  1000 Format (/1X, 'Coefficient optimization and maximum search for ', 'subprocess no', I4/1X, 'Coefficient modes     tau', 10X, 'y*', 9X, 'cth', 9X, 'tau''', 7X, 'sigma')
  1100 Format (1X, 4I4, F12.8, F12.6, F12.7, F12.8, 1P, E12.4)
  1200 Format (1X, 'Error: requested subprocess ', I3, ' has vanishing ', 'cross-section.'/1X, 'Execution stopped!')
  1300 Format (1X, 'Coefficients of equation system to be solved for ', A4)
  1400 Format (1X, 1P, 7E11.3)
  1500 Format (1X, 'Result for ', A4, ':', 6F9.4)
  1600 Format (1X, 'Maximum search for given coefficients'/2X, 'MAX VAR ', 'MOD MOV   VNEW', 7X, 'tau', 7X, 'y*', 8X, 'cth', 7X, 'tau''', 7X, 'sigma')
  1700 Format (1X, 4I4, F8.4, F11.7, F9.3, F11.6, F11.7, 1P, E12.4)
  1800 Format (/1X, 8('*'), 1X, 'PYMAXI: summary of differential ', 'cross-section maximum search', 1X, 8('*'))
  1900 Format (/11X, 58('=')/11X, 'I', 38X, 'I', 17X, 'I'/11X, 'I  ISUB  ', 'Subprocess name', 15X, 'I  Maximum value  I'/11X, 'I', 38X, 'I', 17X, 'I'/11X, 58('=')/11X, 'I', 38X, 'I', 17X, 'I')
  2000 Format (11X, 'I', 2X, I3, 3X, A28, 2X, 'I', 2X, 1P, E12.4, 3X, 'I')
  2100 Format (11X, 'I', 38X, 'I', 17X, 'I'/11X, 58('='))
End Subroutine pymaxi
