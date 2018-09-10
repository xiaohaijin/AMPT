!> @file
!!      AMPT模型的主驱动程序，该subroutine主要读取外部参数和进行驱动。
!> @date
!!      2018 年07月04日
Program ampt
  Double Precision xmp, xmu, alpha, rscut2, cutof2, dshadow
  Double Precision smearp, smearh, dpcoal, drcoal, ecritl
  Character frame*8, proj*8, targ*8
  Character *25 amptvn
  !>==============================================================
  !> @param
  !!       eatt 所产生的粒子的总能量。
  !> @param
  !!       jatt 当前事件所产生的硬散射的数量。
  !> @param
  !!       natt 产生的粒子总数。
  !> @param
  !!       nt,np 靶核和弹核的核子数量。
  !> @param
  !!       n0,n01,n10,n11 参加反应的核子的状态统计数。
  !>==============================================================
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  !>==============================================================
  !> @param
  !!       HIPR1,IHPR2 其中HIPR1，IHPR2包含一些输入参数信息。
  !> @param
  !!       HINT1,IHNT2 包含当前事件的一些额外信息。
  !>==============================================================
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arout/iout
  Common /arevt/iaevt, iarun, miss
  Common /smearz/smearp, smearh
  Common /rndf77/nseed
  Common /anim/nevent, isoft, isflag, izpc
  !>==============================================================
  !> @param
  !!       coal  部分子的组合半径
  !>==============================================================
  Common /coal/dpcoal, drcoal, ecritl
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !>==============================================================
  !> @param
  !!       部分子散射的初始值
  !>==============================================================
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /para8/idpert, npertd, idxsec
  Common /rndm3/iseedp
  !>==============================================================
  !> @param
  !!       强子散射的初始值
  !>==============================================================
  Common /run/num
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /oscar1/iap, izp, iat, izt
  Common /oscar2/frame, amptvn
  Common /resdcy/nsav, iksdcy
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Common /cmsflag/dshadow, ishadow
  !>==============================================================
  !> @param
  !!     允许反应平面进行旋转
  !>==============================================================
  Common /phihj/iphirp, phirp
  External hidata, pydata, ludata, ardata, ppbdat, zpcbdt
  Save
  !>==============================================================
  !> 打开输入文件，读取外部参数。
  !>==============================================================
  Open (24, File='../input/input.ampt', Status='UNKNOWN')
  Open (12, File='../data/version', Status='UNKNOWN')
  Read (24, *) efrm
  Read (24, 111) frame
  Read (24, 111) proj
  Read (24, 111) targ
  Read (24, *) iap
  Read (24, *) izp
  Read (24, *) iat
  Read (24, *) izt
  Read (24, *) nevnt
  Read (24, *) bmin
  Read (24, *) bmax
  Read (24, *) isoft
  Read (24, *) ntmax
  Read (24, *) dt
  Read (24, *) parj(41)
  Read (24, *) parj(42)
  Read (24, *) ipop
  If (ipop==1) ihpr2(11) = 3
  Read (24, *) parj(5)
  Read (24, *) ihpr2(6)
  Read (24, *) ihpr2(4)
  Read (24, *) hipr1(14)
  Read (24, *) hipr1(8)
  Read (24, *) xmu
  Read (24, *) izpc
  Read (24, *) alpha
  Read (24, *) dpcoal
  Read (24, *) drcoal
  Read (24, *) ihjsed
  Read (24, *) nseed
  Read (24, *) iseedp
  Read (24, *) iksdcy
  Read (24, *) iphidcy
  Read (24, *) ipi0dcy
  Read (24, *) ioscar
  Read (24, *) idpert
  Read (24, *) npertd
  Read (24, *) idxsec
  Read (24, *) pttrig
  Read (24, *) maxmiss
  Read (24, *) ihpr2(2)
  Read (24, *) ihpr2(5)
  Read (24, *) iembed
  Read (24, *) pxqembd, pyqembd
  Read (24, *) xembd, yembd
  Read (24, *) nsembd, psembd, tmaxembd
  Read (24, *) ishadow
  Read (24, *) dshadow
  Read (24, *) iphirp
  Close (24)


  !>==============================================================
  !> 版本检测
  !>==============================================================
  If (isoft==1) Then
     amptvn = '1.26t7 (Default)'
  Else If (isoft==4) Then
     amptvn = '2.26t7 (StringMelting)'
  Else
     amptvn = 'Test-Only'
  End If


  !>==============================================================
  !> 时间型随机数的种子
  !>==============================================================
  If (ihjsed==11) Then
     Print *, '# Read in NSEED in HIJING at run time (e.g. 20030819):'
  End If
  Read (*, *) nseedr
  If (ihjsed==11) Then
     nseed = nseedr
  End If
  If (ihjsed==11) Then
     Print *, '#   read in: ', nseed
     Write (12, *) '# Read in NSEED in HIJING at run time:', nseed
  End If


  Close (12)


  !>==============================================================
  !> 随机数控制
  !>==============================================================
  nseed = 2*nseed + 1
  Call srand(nseed)


  !>==============================================================
  !> 初始设定
  !>==============================================================
  ihpr2(10) = 1
  arpar1(1) = 0.7
  smearp = 0D0
  iamax = max(iap, iat)
  smearh = 1.2D0*iamax**0.3333D0/(dble(efrm)/2/0.938D0)
  nevent = nevnt


  !>==============================================================
  !> 输出结果的文件
  !>==============================================================
  Open (16, File='../data/ampt.dat', Status='UNKNOWN')
  Open (14, File='../data/zpc.dat', Status='UNKNOWN')


  !>==============================================================
  !> 初始化HIJING.
  !>==============================================================
  Call hijset(efrm, frame, proj, targ, iap, izp, iat, izt)
  Call artset
  Call inizpc
  Do j = 1, nevnt
     iaevt = j
     Do k = 1, num
        iarun = k
        If (iaevt==nevnt .And. iarun==num) Then
           iout = 1
        End If
        Print *, ' EVENT ', j, ', RUN ', k
        imiss = 0
100     Call hijing(frame, bmin, bmax)
        iaint2(1) = natt
        If (j==-2) Then
           Write (98, *) hipr1
           Write (98, *) ' '
           Write (98, *) ihpr2
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=1, 20)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=21, 40)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=41, 60)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=61, 80)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=81, 100)
           Write (98, *) ' '
           Write (98, *) ihnt2
        End If
        Call getnp
        If (ihpr2(20)==0) Goto 2000
        If (natt==0) Then
           imiss = imiss + 1
           If (imiss<=20) Then
              Write (6, *) 'repeated event: natt=0,j,imiss=', j, imiss
              Goto 100
           Else
              Write (6, *) 'missed event: natt=0,j=', j
              Goto 2000
           End If
        End If
        Call arini
        Call arini2(k)
     End Do
     Call artan1
     Call artmn
     Call artan2
2000 End Do
  Call artout(nevnt)
  Stop
111 Format (A8)
End Program ampt
