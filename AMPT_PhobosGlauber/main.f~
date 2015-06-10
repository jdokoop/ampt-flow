c.....driver program for A Multi-Phase Transport model
      PROGRAM AMPT
c
      double precision xmp, xmu, alpha, rscut2, cutof2
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER FRAME*8, PROJ*8, TARG*8
      character*25 amptvn
      COMMON /HMAIN1/ NATT, EATT, JATT, NT, NP, N0, N01, N10, N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
      common/phidcy/iphidcy
      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************
      OPEN (24, FILE = 'input.ampt', STATUS = 'UNKNOWN')
      OPEN (12, FILE = 'ana/version', STATUS = 'UNKNOWN')
      READ (24, *) EFRM
c     format-read characters (for ALPHA compilers):
      READ (24, 111) FRAME
      READ (24, 111) PROJ
      READ (24, 111) TARG
      READ (24, *) IAP
      READ (24, *) IZP
      READ (24, *) IAT
      READ (24, *) IZT
      READ (24, *) NEVNT
      READ (24, *) BMIN
      READ (24, *) BMAX
c     flag to select default AMPT or string melting:
      READ (24, *) isoft
c     read initialization value for hadron cascade:
      READ (24, *) NTMAX
      READ (24, *) DT
c     parj(41) and (42) are a and b parameters in Lund string fragmentation:
      READ (24, *) PARJ(41)
      READ (24, *) PARJ(42)
c     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and 
c     increase the net-baryon stopping in rapidity (value HIJING is 1):
      READ (24, *) ipop
      if(ipop.eq.1) IHPR2(11)=3
c     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
      READ (24, *) PARJ(5)
c     shadowing flag in HIJING:
      READ (24, *) IHPR2(6)
c     quenching flag in HIJING:
      READ (24, *) IHPR2(4)
c     quenching rate when quenching flag is on (=1.0 GeV/fm):
      READ (24, *) HIPR1(14)
c     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. 
      READ (24, *) HIPR1(8)
c     read initialization value for parton cascade:
      READ (24, *) xmu
      READ (24, *) izpc
      READ (24, *) alpha
c     quark coalescence radii in momentum and space for string melting:
      READ (24, *) dpcoal
      READ (24, *) drcoal
c     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
      READ (24, *) ihjsed
c     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
      READ (24, *) nseed
      READ (24, *) iseedp
      READ (24, *) iksdcy
      READ (24, *) iphidcy
c     flag for OSCAR output for final partons and hadrons:
      READ (24, *) ioscar
clin-5/2008 disabled flags for perturbative treatment of deuterons:
      idpert=0
      npertd=1
      idxsec=1
c
      CLOSE (24)
 111  format(a8)
c
      if(isoft.eq.1) then
         amptvn = '1.21 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.21 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
      WRITE(12,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#               Version ',a20,             '     #'/1X,
     &10X,'#                10/01/2008                      #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')
c     when ihjsed=11: use environment variable at run time for HIJING nseed:
      if(ihjsed.eq.11) then
         PRINT *,
     1 '# Read in NSEED in HIJING at run time (e.g. 20030819):'
      endif
      READ (*, *) nseedr
      if(ihjsed.eq.11) then
         nseed=nseedr
      endif
c     an odd number is needed for the random number generator:
      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
      if(ihjsed.eq.11) then      
         PRINT *, '#   read in: ', nseed
         WRITE(12,*) '# Read in NSEED in HIJING at run time:',nseed
      endif
      CLOSE(12)
c     9/26/03 random number generator for f77 compiler:
      CALL SRAND(NSEED)
c
c.....turn on warning messages
      IHPR2(10) = 1
c     string formation time:
      ARPAR1(1) = 0.7
c     smearp is the smearing halfwidth on parton z0, 
c     set to 0 for now to avoid overflow in eta.
c     smearh is the smearing halfwidth on string production point z0.
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
      nevent=NEVNT
c
c     AMPT momentum and space info at freezeout:
      OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      OPEN (14, FILE = 'ana/zpc.dat', STATUS = 'UNKNOWN')
      open (77, FILE = 'ana/geometry.dat', STATUS = 'UNKNOWN')
ctest off for resonance (phi, K*) studies:
c      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
c      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
      CALL INIZPC
c
       DO 2000 J = 1, NEVNT
          IAEVT = J
          DO 1000 K = 1, NUM
             IARUN = K
             IF (IAEVT .EQ. NEVNT .AND. IARUN .EQ. NUM) THEN
                IOUT = 1
             END IF
             PRINT *, ' EVENT ', J, ', RUN ', K
             imiss=0
 100         CALL HIJING(FRAME, BMIN, BMAX)
             IAINT2(1) = NATT             
c     evaluate Npart (from primary NN collisions) for both proj and targ:
             call getnp
c     switch for final parton fragmentation:
             IF (IHPR2(20) .EQ. 0) GOTO 2000
c     In the unlikely case of no interaction (even after loop of 20 in HIJING),
c     still repeat the event to get an interaction 
c     (this may have an additional "trigger" effect):
             if(NATT.eq.0) then
                imiss=imiss+1
                if(imiss.le.20) then
                   write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
                   goto 100
                else
                   write(6,*) 'missed event: natt=0,j=',j
                   goto 2000
                endif
             endif
c.....ART initialization and run
             CALL ARINI
             CALL ARINI2(K)
 1000     CONTINUE
c
          CALL ARTAN1
          CALL HJANA3
          CALL ARTMN
          CALL HJANA4
          CALL ARTAN2
 2000  CONTINUE
c
       CALL ARTOUT(NEVNT)
c
       STOP
       END
c     FYI: free file unit numbers are 9,10,20-23,59-88.
