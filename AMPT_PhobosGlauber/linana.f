c.................... linana.f
c=======================================================================
c     10/26/01 update freezeout positions in case of interactions:
      subroutine hbtout(nnew,nt,ntmax)
c
      PARAMETER  (MAXSTR=150001,MAXR=1)
clin-5/2008 give tolerance to regular particles (perturbative probability 1):
      PARAMETER  (oneminus=0.99999,oneplus=1.00001)
      dimension lastkp(MAXSTR), newkp(MAXSTR),xnew(3)
      common /para7/ ioscar
cc      SAVE /para7/
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
cc      SAVE /hbt/
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      COMMON   /AA/  R(3,MAXSTR)
cc      SAVE /AA/
      COMMON   /BB/  P(3,MAXSTR)
cc      SAVE /BB/
      COMMON   /CC/  E(MAXSTR)
cc      SAVE /CC/
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
      common /lastt/itimeh,bimp
cc      SAVE /lastt/
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
      COMMON /AREVT/ IAEVT, IARUN
cc      SAVE /AREVT/
      common/snn/efrm,npart1,npart2
cc      SAVE /snn/
      COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
cc      SAVE /HJGLBR/
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
clin-12/14/03:
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      EXTERNAL IARFLV, INVFLV
      common /para8/ idpert,npertd,idxsec
      SAVE   
c
      do 1001 i=1,max0(nlast,nnew)
         lastkp(i)=0
 1001 continue
      do 1002 i=1,nnew
         newkp(i)=0
 1002 continue
c     for each of the particles, search the freezeout record (common /hbt/) 
c     to find & keep those which do not have interactions during this timestep:
      do 100 ip=1,nnew
         do 1004 iplast=1,nlast
            if(p(1,ip).eq.plast(1,iplast).and.
     1           p(2,ip).eq.plast(2,iplast).and.
     2           p(3,ip).eq.plast(3,iplast).and.
     3           e(ip).eq.plast(4,iplast).and.
     4           lb(ip).eq.lblast(iplast).and.
     5      dpertp(ip).eq.dplast(iplast).and.lastkp(iplast).eq.0) then
clin-5/2008 modified below to the above in case we have perturbative particles:
c     5           lastkp(iplast).eq.0) then
               deltat=nt*dt-xlast(4,iplast)
               ene=sqrt(plast(1,iplast)**2+plast(2,iplast)**2
     1              +plast(3,iplast)**2+plast(4,iplast)**2)
c     xnew gives the coordinate if a particle free-streams to current time:
               do 1003 ii=1,3
                  xnew(ii)=xlast(ii,iplast)+plast(ii,iplast)/ene*deltat
 1003          continue
                  dr=sqrt((r(1,ip)-xnew(1))**2+(r(2,ip)-xnew(2))**2
     1              +(r(3,ip)-xnew(3))**2)
c     find particles with dp=0 and dr<0.01, considered to be those 
c     without any interactions during this timestep, 
c     thus keep their last positions and time:
               if(dr.le.0.01) then
                  lastkp(iplast)=1
                  newkp(ip)=1
ctest off: write collision info
c                  write(9,*) 'nt,ip,px,x=',nt,ip,p(1,ip),r(1,ip)
                  goto 100
               endif
            endif
 1004    continue
 100  continue
c     for current particles with interactions, fill their current info in 
c     the freezeout record (if that record entry needs not to be kept):
      do 150 ip=1,nnew
         if(newkp(ip).eq.0) then
            do 1005 iplast=1,nnew
               if(lastkp(iplast).eq.0) then
ctest off: write collision info
c                  if(nt.ge.(ntmax-5)) then
c                     write(95,*) 'nt,lb(ip)=',nt,lb(ip)
c                  write(95,*) '  last p=',plast(1,iplast),
c     1 plast(2,iplast),plast(3,iplast),plast(4,iplast)
c                  write(95,*) '  after p=',p(1,ip),p(2,ip),p(3,ip),e(ip)
c                  write(95,*) '  after x=',r(1,ip),r(2,ip),r(3,ip)
c                  endif
                  xlast(1,iplast)=r(1,ip)
                  xlast(2,iplast)=r(2,ip)
                  xlast(3,iplast)=r(3,ip)
                  xlast(4,iplast)=nt*dt
c
                  if(nt.eq.ntmax) then
c     freezeout time for decay daughters at the last timestep 
c     needs to include the decay time of the parent:
                     if(tfdcy(ip).gt.(ntmax*dt+0.001)) then
                        xlast(4,iplast)=tfdcy(ip)
c     freezeout time for particles unformed at the next-to-last timestep 
c     needs to be their formation time instead of (ntmax*dt):
                     elseif(ftsv(ip).gt.((ntmax-1)*dt)) then
                        xlast(4,iplast)=ftsv(ip)
                     endif
                  endif
                  plast(1,iplast)=p(1,ip)
                  plast(2,iplast)=p(2,ip)
                  plast(3,iplast)=p(3,ip)
                  plast(4,iplast)=e(ip)
                  lblast(iplast)=lb(ip)
                  lastkp(iplast)=1
clin-5/2008:
                  dplast(iplast)=dpertp(ip)
                  goto 150
               endif
 1005       continue
         endif
 150  continue
c     if the current particle list is shorter than the freezeout record,
c     condense the last-collision record by filling new record from 1 to nnew, 
c     and label these entries as keep:
      if(nnew.lt.nlast) then
         do 170 iplast=1,nlast
            if(lastkp(iplast).eq.0) then
               do 1006 ip2=iplast+1,nlast
                  if(lastkp(ip2).eq.1) then
                     xlast(1,iplast)=xlast(1,ip2)
                     xlast(2,iplast)=xlast(2,ip2)
                     xlast(3,iplast)=xlast(3,ip2)
                     xlast(4,iplast)=xlast(4,ip2)
                     plast(1,iplast)=plast(1,ip2)
                     plast(2,iplast)=plast(2,ip2)
                     plast(3,iplast)=plast(3,ip2)
                     plast(4,iplast)=plast(4,ip2)
                     lblast(iplast)=lblast(ip2)
                     lastkp(iplast)=1
clin-5/2008:
                     dplast(iplast)=dplast(ip2)
                     goto 170
                  endif
 1006          continue
            endif
 170     continue
      endif
      nlast=nnew
ctest off look inside each NT timestep (for debugging purpose):
c      do ip=1,nlast
c         if(nt.eq.5000) then
c            write(95,*) ' p ',nt,ip,INVFLV(lblast(ip)),plast(1,ip),
c     1           plast(2,ip),plast(3,ip),plast(4,ip),dplast(ip)
c            write(95,*) '  x ',nt,ip,INVFLV(lblast(ip)),xlast(1,ip),
c     1           xlast(2,ip),xlast(3,ip),xlast(4,ip),dplast(ip)
c         endif
c      enddo
c
      if(nt.eq.ntmax) then
clin-5/2008 find final number of perturbative particles (deuterons only):
         ndpert=0
         do ip=1,nlast
            if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            else
               ndpert=ndpert+1
            endif
         enddo
c
c         write(16,190) IAEVT,IARUN,nlast,bimp,npart1,npart2,
c     1 NELP,NINP,NELT,NINTHJ
         write(16,190) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
     1 NELP,NINP,NELT,NINTHJ
clin-5/2008 write out perturbatively-produced particles (deuterons only):
         if(idpert.eq.1.or.idpert.eq.2)
     1        write(90,190) IAEVT,IARUN,ndpert,bimp,npart1,npart2,
     2        NELP,NINP,NELT,NINTHJ
         do 1007 ip=1,nlast
clin-12/14/03   No formation time for spectator projectile or target nucleons,
c     see ARINI1 in 'amptsub.f':
            if(plast(1,ip).eq.0.and.plast(2,ip).eq.0
     1           .and.(sqrt(plast(3,ip)**2+plast(4,ip)**2)*2/HINT1(1))
     2           .gt.0.99.and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
clin-5/2008 perturbatively-produced particles (currently only deuterons) 
c     are written to ana/ampt_pert.dat (without the column for the mass); 
c     ana/ampt.dat has regularly-produced particles (including deuterons);
c     these two sets of deuteron data are close to each other(but not the same 
c     because of the bias from triggering the perturbative production); 
c     ONLY use one data set for analysis to avoid double-counting:
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
                  write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),plast(4,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 1.E-20*sqrt(plast(3,ip)**2+plast(4,ip)**2)
     4                 /plast(4,ip)
clin-12/14/03-end
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,250) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 1.E-20*sqrt(plast(3,ip)**2+plast(4,ip)**2)
     4                 /plast(4,ip),dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            elseif(amax1(abs(xlast(1,ip)),abs(xlast(2,ip)),
     1              abs(xlast(3,ip)),abs(xlast(4,ip))).lt.9999) then
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
            write(90,250) INVFLV(lblast(ip)),plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            else
c     change format for large numbers:
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,201) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,251) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
              endif
           endif
 1007    continue
         if(ioscar.eq.1) call hoscar
      endif
 190  format(3(i7),f10.4,5x,6(i4))
 200  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
 250  format(I5,2(1x,f8.3),1x,f10.3,2(1x,f7.1),1x,f8.2,1x,f7.2,1x,e9.3)
 251  format(I5,2(1x,f8.3),1x,f10.3,4(1x,e8.2),1x,e9.3)
c     
        return
        end

c=======================================================================
        SUBROUTINE decomp(px0,py0,pz0,xm0)
c
        IMPLICIT DOUBLE PRECISION(D)  
        DOUBLE PRECISION  enenew, pxnew, pynew, pznew
        DOUBLE PRECISION  de0, beta2, gam
        common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        common /decom/ptwo(2,5)
cc      SAVE /decom/
        COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        SAVE   
c
        dcth=dble(RANART(NSEED))*2.d0-1.d0
        dPHI=dble(RANART(NSEED)*HIPR1(40))*2.d0
        ds=dble(xm0)**2
        dpcm=dsqrt((ds-dble(ptwo(1,5)+ptwo(2,5))**2)
     1 *(ds-dble(ptwo(1,5)-ptwo(2,5))**2)/ds/4d0)
        dpz=dpcm*dcth
        dpx=dpcm*dsqrt(1.d0-dcth**2)*dcos(dphi)
        dpy=dpcm*dsqrt(1.d0-dcth**2)*dsin(dphi)
        de1=dsqrt(dble(ptwo(1,5))**2+dpcm**2)
        de2=dsqrt(dble(ptwo(2,5))**2+dpcm**2)
c
      de0=dsqrt(dble(px0)**2+dble(py0)**2+dble(pz0)**2+dble(xm0)**2)
        dbex=dble(px0)/de0
        dbey=dble(py0)/de0
        dbez=dble(pz0)/de0
c     boost the reference frame up by beta (pznew=gam(pz+beta e)):
      beta2 = dbex ** 2 + dbey ** 2 + dbez ** 2
      gam = 1.d0 / dsqrt(1.d0 - beta2)
      if(beta2.ge.0.9999999999999d0) then
         write(6,*) '1',dbex,dbey,dbez,beta2,gam
      endif
c
      call lorenz(de1,dpx,dpy,dpz,-dbex,-dbey,-dbez)
        ptwo(1,1)=sngl(pxnew)
        ptwo(1,2)=sngl(pynew)
        ptwo(1,3)=sngl(pznew)
        ptwo(1,4)=sngl(enenew)
c
      call lorenz(de2,-dpx,-dpy,-dpz,-dbex,-dbey,-dbez)
        ptwo(2,1)=sngl(pxnew)
        ptwo(2,2)=sngl(pynew)
        ptwo(2,3)=sngl(pznew)
        ptwo(2,4)=sngl(enenew)
c
      RETURN
      END

c=======================================================================
      SUBROUTINE HTOP
c
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXIDL=4001)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      dimension it(4)
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
      COMMON/HMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
      COMMON /PARA1/ MUL
cc      SAVE /PARA1/
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
      common /decom/ptwo(2,5)
cc      SAVE /decom/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &     GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &     PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &     XMN(MAXIDL)
cc      SAVE /NOPREC/
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
c     7/20/01: use double precision
c     otherwise sometimes beta>1 and gamma diverge in lorenz():
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
      DOUBLE PRECISION  vxp0,vyp0,vzp0
      common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
cc      SAVE /precpa/
      common /para7/ ioscar
      SAVE   
c
        npar=0
        nnozpc=0
clin-5b/2008 calculate the number of hadrons to be converted to q/qbar:
       if((isoft.eq.4.or.isoft.eq.5).and.ioscar.eq.2) then
           nsmbbbar=0
           nsmmeson=0
           do i=1,natt
              id=ITYPAR(i)
              idabs=iabs(id)
              i2=MOD(idabs/10,10)
              if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
     1             .ge.(HINT1(1)/2*0.99).and.
     2             ((id.eq.2112).or.(id.eq.2212))) then
c     proj or targ nucleons without interactions, do not enter ZPC:
              elseif(idabs.gt.1000.and.i2.ne.0) then
c     baryons to be converted to q/qbar:
                 nsmbbbar=nsmbbbar+1
              elseif((idabs.gt.100.and.idabs.lt.1000)
     1                .or.idabs.gt.10000) then
c     mesons to be converted to q/qbar:
                 nsmmeson=nsmmeson+1
              endif
           enddo
           write(92,*) 3*nsmbbbar+2*nsmmeson
           write(92,*) ' is the total # of initial partons after string 
     1 melting'
           write(92,*) 'String melting converts ',nsmbbbar, ' baryons &'
     1, nsmmeson, ' mesons'
           write(92,*) 'Total # of initial particles= ',natt
           write(92,*) 'Total # of initial particles (gamma,e,muon,...) 
     1 not entering ZPC= ',natt-nsmbbbar-nsmmeson
        endif
clin-5b/2008-over
        do 100 i=1,natt
           id=ITYPAR(i)
           idabs=iabs(id)
           i4=MOD(idabs/1000,10)
           i3=MOD(idabs/100,10)
           i2=MOD(idabs/10,10)
           i1=MOD(idabs,10)
           rnum=RANART(NSEED)
           ftime=0.197*PEAR(i)/(PXAR(i)**2+PYAR(i)**2+XMAR(i)**2)
           inozpc=0
           it(1)=0
           it(2)=0
           it(3)=0
           it(4)=0
c
           if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
     1 .ge.(HINT1(1)/2*0.99).and.((id.eq.2112).or.(id.eq.2212))) then
c     proj or targ nucleons without interactions, do not enter ZPC:
              inozpc=1
           elseif(idabs.gt.1000.and.i2.ne.0) then
c     baryons:
              if(((i4.eq.1.or.i4.eq.2).and.i4.eq.i3)
     1 .or.(i4.eq.3.and.i3.eq.3)) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 endif
              elseif(i4.eq.1.or.i4.eq.2) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 endif
              elseif(i4.ge.3) then
                 it(1)=i4
                 if(i3.lt.i2) then
                    it(2)=i2*1000+i3*100+1
                 else
                    it(2)=i3*1000+i2*100+3
                 endif
              endif
c       antibaryons:
              if(id.lt.0) then
                 it(1)=-it(1)
                 it(2)=-it(2)
              endif
c     isoft=4or5 decompose diquark flavor it(2) to two quarks it(3)&(4):
              if(isoft.eq.4.or.isoft.eq.5) then
                 it(3)=MOD(it(2)/1000,10)
                 it(4)=MOD(it(2)/100,10)
              endif

           elseif((idabs.gt.100.and.idabs.lt.1000)
     1 .or.idabs.gt.10000) then
c     mesons:
              if(i3.eq.i2) then
                 if(i3.eq.1.or.i3.eq.2) then
                    if(rnum.le.0.5) then
                       it(1)=1
                       it(2)=-1
                    else
                       it(1)=2
                       it(2)=-2
                    endif
                 else
                    it(1)=i3
                    it(2)=-i3
                 endif
              else
                 if((isign(1,id)*(-1)**i3).eq.1) then
                    it(1)=i3
                    it(2)=-i2
                 else
                    it(1)=i2
                    it(2)=-i3
                 endif
              endif
           else
c     save other particles (leptons and photons) outside of ZPC:
              inozpc=1
           endif
c
           if(inozpc.eq.1) then
              NJSGS(i)=0
              nnozpc=nnozpc+1
              itypn(nnozpc)=ITYPAR(i)
              pxn(nnozpc)=PXAR(i)
              pyn(nnozpc)=PYAR(i)
              pzn(nnozpc)=PZAR(i)
              een(nnozpc)=PEAR(i)
              xmn(nnozpc)=XMAR(i)
              gxn(nnozpc)=GXAR(i)
              gyn(nnozpc)=GYAR(i)
              gzn(nnozpc)=GZAR(i)
              ftn(nnozpc)=FTAR(i)
           else
              NJSGS(i)=2
              ptwo(1,5)=ulmass(it(1))
              ptwo(2,5)=ulmass(it(2))
              call decomp(patt(i,1),patt(i,2),patt(i,3),XMAR(i))
              ipamax=2
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) ipamax=1
              do 1001 ipar=1,ipamax
                 npar=npar+1
                 ityp0(npar)=it(ipar)
                 px0(npar)=dble(ptwo(ipar,1))
                 py0(npar)=dble(ptwo(ipar,2))
                 pz0(npar)=dble(ptwo(ipar,3))
                 e0(npar)=dble(ptwo(ipar,4))
                 xmass0(npar)=dble(ptwo(ipar,5))
                 gx0(npar)=dble(GXAR(i))
                 gy0(npar)=dble(GYAR(i))
                 gz0(npar)=dble(GZAR(i))
                 ft0(npar)=dble(ftime)
                 lstrg0(npar)=i
                 lpart0(npar)=ipar
                 vxp0(npar)=dble(patt(i,1)/patt(i,4))
                 vyp0(npar)=dble(patt(i,2)/patt(i,4))
                 vzp0(npar)=dble(patt(i,3)/patt(i,4))
clin-5b/2008:
                 if((isoft.eq.4.or.isoft.eq.5).and.ioscar.eq.2) then
                    if(dmax1(abs(gx0(npar)),abs(gy0(npar)),
     1                   abs(gz0(npar)),abs(ft0(npar))).lt.9999) then
                       write(92,200) ityp0(npar),px0(npar),py0(npar),
     1                      pz0(npar),xmass0(npar),gx0(npar),gy0(npar),
     2                      gz0(npar),ft0(npar)
                    else
                       write(92,201) ityp0(npar),px0(npar),py0(npar),
     1                      pz0(npar),xmass0(npar),gx0(npar),gy0(npar),
     2                      gz0(npar),ft0(npar)
                    endif
                 endif
c
 1001     continue
 200      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
c
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) then
                 NJSGS(i)=3
                 xmdq=ptwo(2,5)
                 ptwo(1,5)=ulmass(it(3))
                 ptwo(2,5)=ulmass(it(4))
c     8/19/02 avoid actual argument in common blocks of DECOMP:
c                 call decomp(ptwo(2,1),ptwo(2,2),ptwo(2,3),xmdq)
             ptwox=ptwo(2,1)
             ptwoy=ptwo(2,2)
             ptwoz=ptwo(2,3)
             call decomp(ptwox,ptwoy,ptwoz,xmdq)
c
                 do 1002 ipar=1,2
                    npar=npar+1
                    ityp0(npar)=it(ipar+2)
                    px0(npar)=dble(ptwo(ipar,1))
                    py0(npar)=dble(ptwo(ipar,2))
                    pz0(npar)=dble(ptwo(ipar,3))
                    e0(npar)=dble(ptwo(ipar,4))
                    xmass0(npar)=dble(ptwo(ipar,5))
                    gx0(npar)=dble(GXAR(i))
                    gy0(npar)=dble(GYAR(i))
                    gz0(npar)=dble(GZAR(i))
                    ft0(npar)=dble(ftime)
                    lstrg0(npar)=i
                    lpart0(npar)=ipar+1
                    vxp0(npar)=dble(patt(i,1)/patt(i,4))
                    vyp0(npar)=dble(patt(i,2)/patt(i,4))
                    vzp0(npar)=dble(patt(i,3)/patt(i,4))
clin-5b/2008:
                    if((isoft.eq.4.or.isoft.eq.5).and.ioscar.eq.2) then
                       if(dmax1(abs(gx0(npar)),abs(gy0(npar)),
     1                      abs(gz0(npar)),abs(ft0(npar))).lt.9999) then
                          write(92,200) ityp0(npar),px0(npar),py0(npar),
     1                         pz0(npar),xmass0(npar),gx0(npar),
     2                         gy0(npar),gz0(npar),ft0(npar)
                       else
                          write(92,201) ityp0(npar),px0(npar),py0(npar),
     1                         pz0(npar),xmass0(npar),gx0(npar),
     2                         gy0(npar),gz0(npar),ft0(npar)
                       endif
                    endif
c
 1002        continue
              endif
c
           endif
 100        continue
      MUL=NPAR
c      
clin-5b/2008:
      if((isoft.eq.4.or.isoft.eq.5).and.ioscar.eq.2) then
         if((natt-nsmbbbar-nsmmeson).ne.nnozpc) 
     1        write(20,*) 'Problem with the total # of initial particles
     2 (gamma,e,muon,...) not entering ZPC'
         if((3*nsmbbbar+2*nsmmeson).ne.npar) 
     1        write(20,*) 'Problem with the total # of initial partons
     2 after string melting'
      endif
c
      RETURN
      END

c=======================================================================
      SUBROUTINE PTOH
c
      PARAMETER (MAXSTR=150001)
      DOUBLE PRECISION  gxp,gyp,gzp,ftp,pxp,pyp,pzp,pep,pmp
      DOUBLE PRECISION  gxp0,gyp0,gzp0,ft0fom,drlocl
      DOUBLE PRECISION  enenew, pxnew, pynew, pznew, beta2, gam
      DOUBLE PRECISION  ftavg0,gxavg0,gyavg0,gzavg0,bex,bey,bez
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      DOUBLE PRECISION  xmdiag,px1,py1,pz1,e1,px2,py2,pz2,e2,
     1     px3,py3,pz3,e3,xmpair,etot
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
cc      SAVE /loclco/
      COMMON/HMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
cc      SAVE /ARPRNT/
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
cc      SAVE /ARPRC/
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
      common /prtn23/ gxp0(3),gyp0(3),gzp0(3),ft0fom
cc      SAVE /prtn23/
      common /nzpc/nattzp
cc      SAVE /nzpc/
      common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
cc      SAVE /LUDAT1/ 
      dimension xmdiag(MAXSTR),indx(MAXSTR),ndiag(MAXSTR)
      SAVE   
c
      call coales
c     obtain particle mass here without broadening by Breit-Wigner width:
      mstj24=MSTJ(24)
      MSTJ(24)=0
        nuudd=0
        npich=0
        nrhoch=0
      ppi0=1.
      prho0=0.
c     determine hadron flavor except (pi0,rho0,eta,omega):
      DO 1001 ISG = 1, NSG
           if(NJSGS(ISG).ne.0) then
              NATT=NATT+1
              K1=K2SGS(ISG,1)
              k1abs=iabs(k1)
              PX1=PXSGS(ISG,1)
              PY1=PYSGS(ISG,1)
              PZ1=PZSGS(ISG,1)
              K2=K2SGS(ISG,2)
              k2abs=iabs(k2)
              PX2=PXSGS(ISG,2)
              PY2=PYSGS(ISG,2)
              PZ2=PZSGS(ISG,2)
c     5/02/01 try lowest spin states as first choices, 
c     i.e. octet baryons and pseudoscalar mesons (ibs=2*baryonspin+1):
              e1=PESGS(ISG,1)
              e2=PESGS(ISG,2)
              xmpair=dsqrt((e1+e2)**2-(px1+px2)**2-(py1+py2)**2
     1 -(pz1+pz2)**2)
              ibs=2
              imspin=0
              if(k1.eq.-k2.and.iabs(k1).le.2.
     1           and.NJSGS(ISG).eq.2) then
               nuudd=nuudd+1
               xmdiag(nuudd)=xmpair
               ndiag(nuudd)=natt
            endif
              K3=0
              if((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) then
               K3=K2SGS(ISG,3)
               k3abs=iabs(k3)
               PX3=PXSGS(ISG,3)
               PY3=PYSGS(ISG,3)
               PZ3=PZSGS(ISG,3)
               e3=PESGS(ISG,3)
               xmpair=dsqrt((e1+e2+e3)**2-(px1+px2+px3)**2
     1              -(py1+py2+py3)**2-(pz1+pz2+pz3)**2)
              endif
c*****     isoft=3 baryon decomposition is different:
              if(isoft.eq.3.and.
     1           (k1abs.gt.1000.or.k2abs.gt.1000)) then
               if(k1abs.gt.1000) then
                  kdq=k1abs
                  kk=k2abs
               else
                  kdq=k2abs
                  kk=k1abs
               endif
               ki=MOD(kdq/1000,10)
               kj=MOD(kdq/100,10)
               if(MOD(kdq,10).eq.1) then
                  idqspn=0
               else
                  idqspn=1
               endif
c
               if(kk.gt.ki) then
                  ktemp=kk
                  kk=kj
                  kj=ki
                  ki=ktemp
               elseif(kk.gt.kj) then
                  ktemp=kk
                  kk=kj
                  kj=ktemp
               endif
c     
               if(ki.ne.kj.and.ki.ne.kk.and.kj.ne.kk) then
                  if(idqspn.eq.0) then
                     kf=1000*ki+100*kk+10*kj+ibs
                  else
                     kf=1000*ki+100*kj+10*kk+ibs
                  endif
               elseif(ki.eq.kj.and.ki.eq.kk) then
c     can only be decuplet baryons:
                  kf=1000*ki+100*kj+10*kk+4
               else
                  kf=1000*ki+100*kj+10*kk+ibs
               endif
c     form a decuplet baryon if the q+diquark mass is closer to its mass 
c     (and if the diquark has spin 1):
cc     for now only include Delta, which is present in ART:
cc                 if(idqspn.eq.1.and.MOD(kf,10).eq.2) then
               if(kf.eq.2112.or.kf.eq.2212) then
                  if(abs(sngl(xmpair)-ULMASS(kf)).gt.
     1                 abs(sngl(xmpair)-ULMASS(kf+2))) kf=kf+2
               endif
               if(k1.lt.0) kf=-kf
clin-6/22/01 isoft=4or5 baryons:
              elseif((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) 
     1              then
               if(k1abs.gt.k2abs) then
                  ki=k1abs
                  kk=k2abs
               else
                  ki=k2abs
                  kk=k1abs
               endif
               if(k3abs.gt.ki) then
                  kj=ki
                  ki=k3abs
               elseif(k3abs.lt.kk) then
                  kj=kk
                  kk=k3abs
               else
                  kj=k3abs
               endif
c     
               if(ki.eq.kj.and.ki.eq.kk) then
c     can only be decuplet baryons (Delta-,++, Omega):
                  ibs=4
                  kf=1000*ki+100*kj+10*kk+ibs
               elseif(ki.ne.kj.and.ki.ne.kk.and.kj.ne.kk) then
c     form Lambda or Sigma according to 3-quark mass, 
c     for now neglect decuplet (Sigma*0 etc) which is absent in ART:
                  ibs=2
                  kf1=1000*ki+100*kj+10*kk+ibs
                  kf2=1000*ki+100*kk+10*kj+ibs
                  kf=kf1
                  if(abs(sngl(xmpair)-ULMASS(kf1)).gt.
     1                 abs(sngl(xmpair)-ULMASS(kf2))) kf=kf2
               else
                  ibs=2
                  kf=1000*ki+100*kj+10*kk+ibs
cc     for now only include Delta0,+ as decuplets, which are present in ART:
                  if(kf.eq.2112.or.kf.eq.2212) then
                     if(abs(sngl(xmpair)-ULMASS(kf)).gt.
     1                    abs(sngl(xmpair)-ULMASS(kf+2))) kf=kf+2
                  endif
               endif
               if(k1.lt.0) kf=-kf
c*****     mesons:
              else
               if(k1abs.eq.k2abs) then
                  if(k1abs.le.2) then
c     treat diagonal mesons later in the subroutine:
                     kf=0
                  elseif(k1abs.le.3) then
c     do not form eta', only form phi from s-sbar, since no eta' in ART:
                     kf=333
                  else
                     kf=100*k1abs+10*k1abs+2*imspin+1
                  endif
               else
                  if(k1abs.gt.k2abs) then
                     kmax=k1abs
                     kmin=k2abs
                  elseif(k1abs.lt.k2abs) then
                     kmax=k2abs
                     kmin=k1abs
                  endif
                  kf=(100*kmax+10*kmin+2*imspin+1)
     1                 *isign(1,k1+k2)*(-1)**kmax
c     form a vector meson if the q+qbar mass is closer to its mass:
                  if(MOD(iabs(kf),10).eq.1) then
                     if(abs(sngl(xmpair)-ULMASS(iabs(kf))).gt.
     1                    abs(sngl(xmpair)-ULMASS(iabs(kf)+2))) 
     2                    kf=(iabs(kf)+2)*isign(1,kf)
                  endif
               endif
              endif
              ITYPAR(NATT)=kf
              KATT(NATT,1)=kf
            if(iabs(kf).eq.211) then
               npich=npich+1
            elseif(iabs(kf).eq.213) then
               nrhoch=nrhoch+1
            endif
           endif
 1001   CONTINUE
c     assume Npi0=(Npi+ + Npi-)/2, Nrho0=(Nrho+ + Nrho-)/2 on the average:
        if(nuudd.ne.0) then
         ppi0=float(npich/2)/float(nuudd)
         prho0=float(nrhoch/2)/float(nuudd)
      endif      
c     determine diagonal mesons (pi0,rho0,eta and omega) from uubar/ddbar:
      npi0=0
      DO 1002 ISG = 1, NSG
         if(K2SGS(ISG,1).eq.-K2SGS(ISG,2)
     1        .and.iabs(K2SGS(ISG,1)).le.2.and.NJSGS(ISG).eq.2) then
            if(RANART(NSEED).le.ppi0) npi0=npi0+1
         endif
 1002 CONTINUE
c
      if(nuudd.gt.1) then
         call index1(MAXSTR,nuudd,xmdiag,indx)
      else
         indx(1)=1
      end if
c
      DO 1003 ix=1,nuudd
         iuudd=indx(ix)
         inatt=ndiag(iuudd)            
         if(ix.le.npi0) then
            kf=111
         elseif(RANART(NSEED).le.(prho0/(1-ppi0+0.00001))) then
            kf=113
         else
c     at T=150MeV, thermal weights for eta and omega(spin1) are about the same:
            if(RANART(NSEED).le.0.5) then
               kf=221
            else
               kf=223
            endif
         endif
         ITYPAR(inatt)=kf
         KATT(inatt,1)=kf
 1003 CONTINUE
c  determine hadron formation time, position and momentum:
      inatt=0
      DO 1006 ISG = 1, NSG
           if(NJSGS(ISG).ne.0) then
            inatt=inatt+1
              K1=K2SGS(ISG,1)
              k1abs=iabs(k1)
              PX1=PXSGS(ISG,1)
              PY1=PYSGS(ISG,1)
              PZ1=PZSGS(ISG,1)
              K2=K2SGS(ISG,2)
              k2abs=iabs(k2)
              PX2=PXSGS(ISG,2)
              PY2=PYSGS(ISG,2)
              PZ2=PZSGS(ISG,2)
              e1=PESGS(ISG,1)
              e2=PESGS(ISG,2)
c
              if(NJSGS(ISG).eq.2) then
               PXAR(inatt)=sngl(px1+px2)
               PYAR(inatt)=sngl(py1+py2)
               PZAR(inatt)=sngl(pz1+pz2)
               PATT(inatt,1)=PXAR(inatt)
               PATT(inatt,2)=PYAR(inatt)
               PATT(inatt,3)=PZAR(inatt)
               etot=e1+e2
              elseif((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) 
     1              then
               PX3=PXSGS(ISG,3)
               PY3=PYSGS(ISG,3)
               PZ3=PZSGS(ISG,3)
               e3=PESGS(ISG,3)
               PXAR(inatt)=sngl(px1+px2+px3)
               PYAR(inatt)=sngl(py1+py2+py3)
               PZAR(inatt)=sngl(pz1+pz2+pz3)
               PATT(inatt,1)=PXAR(inatt)
               PATT(inatt,2)=PYAR(inatt)
               PATT(inatt,3)=PZAR(inatt)
               etot=e1+e2+e3
              endif
              XMAR(inatt)=ULMASS(ITYPAR(inatt))
              PEAR(inatt)=sqrt(PXAR(inatt)**2+PYAR(inatt)**2
     1           +PZAR(inatt)**2+XMAR(inatt)**2)
              PATT(inatt,4)=PEAR(inatt)
              EATT=EATT+PEAR(inatt)
            ipartn=NJSGS(ISG)
            DO 1004 i=1,ipartn
               ftp(i)=ftsgs(isg,i)
               gxp(i)=gxsgs(isg,i)
               gyp(i)=gysgs(isg,i)
               gzp(i)=gzsgs(isg,i)
               pxp(i)=pxsgs(isg,i)
               pyp(i)=pysgs(isg,i)
               pzp(i)=pzsgs(isg,i)
               pmp(i)=pmsgs(isg,i)
               pep(i)=pesgs(isg,i)
 1004       CONTINUE
            call locldr(ipartn,drlocl)
c
            tau0=ARPAR1(1)
            ftavg0=ft0fom+dble(tau0)
            gxavg0=0d0
            gyavg0=0d0
            gzavg0=0d0
            DO 1005 i=1,ipartn
               gxavg0=gxavg0+gxp0(i)/ipartn
               gyavg0=gyavg0+gyp0(i)/ipartn
               gzavg0=gzavg0+gzp0(i)/ipartn
 1005       CONTINUE
            bex=dble(PXAR(inatt))/etot
            bey=dble(PYAR(inatt))/etot
            bez=dble(PZAR(inatt))/etot
            beta2 = bex ** 2 + bey ** 2 + bez ** 2
            gam = 1.d0 / dsqrt(1.d0 - beta2)
            if(beta2.ge.0.9999999999999d0) then
               write(6,*) '2',bex,bey,bez,beta2,gam
            endif
c
            call lorenz(ftavg0,gxavg0,gyavg0,gzavg0,-bex,-bey,-bez)
              GXAR(inatt)=sngl(pxnew)
              GYAR(inatt)=sngl(pynew)
              GZAR(inatt)=sngl(pznew)
              FTAR(inatt)=sngl(enenew)
           endif
 1006   CONTINUE
c     number of hadrons formed from partons inside ZPC:
      nattzp=natt
      MSTJ(24)=mstj24
c      
      RETURN
      END

c=======================================================================
      SUBROUTINE coales

      PARAMETER (MAXSTR=150001)
      IMPLICIT DOUBLE PRECISION(D)
      DOUBLE PRECISION  gxp,gyp,gzp,ftp,pxp,pyp,pzp,pep,pmp
      DIMENSION IOVER(MAXSTR),dp1(2:3),dr1(2:3)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      double precision  dpcoal,drcoal,ecritl
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      common /coal/dpcoal,drcoal,ecritl
cc      SAVE /coal/
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
cc      SAVE /loclco/
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
      SAVE   
c      
      do 1001 ISG=1, NSG
         IOVER(ISG)=0
 1001 continue
C1     meson q coalesce with all available qbar:
      do 150 ISG=1,NSG
         if(NJSGS(ISG).ne.2.or.IOVER(ISG).eq.1) goto 150
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
         if(K2SGS(ISG,1).lt.0) then
            write(6,*) 'Antiquark appears in quark loop; stop'
            stop
         endif
c         
         do 1002 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1002    continue
         call locldr(2,drlocl)
         dr0=drlocl
c     dp0^2 defined as (p1+p2)^2-(m1+m2)^2:
         dp0=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 120 JSG=1,NSG
c     skip default or unavailable antiquarks:
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 120
            if(NJSGS(JSG).eq.2) then
               ipmin=2
               ipmax=2
            elseif(NJSGS(JSG).eq.3.and.K2SGS(JSG,1).lt.0) then
               ipmin=1
               ipmax=3
            else
               goto 120
            endif
            do 100 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(1)*pesgs(jsg,ip)
     1              -pxp(1)*pxsgs(jsg,ip)
     2              -pyp(1)*pysgs(jsg,ip)
     3              -pzp(1)*pzsgs(jsg,ip)
     4              -pmp(1)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
               if(dplocl.gt.dpcoal) goto 120
               ftp(2)=ftsgs(jsg,ip)
               gxp(2)=gxsgs(jsg,ip)
               gyp(2)=gysgs(jsg,ip)
               gzp(2)=gzsgs(jsg,ip)
               pxp(2)=pxsgs(jsg,ip)
               pyp(2)=pysgs(jsg,ip)
               pzp(2)=pzsgs(jsg,ip)
               pmp(2)=pmsgs(jsg,ip)
               pep(2)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
               if(drlocl.gt.drcoal) goto 120
c     q_isg coalesces with qbar_jsg:
               if((dp0.gt.dpcoal.or.dr0.gt.drcoal)
     1              .or.(drlocl.lt.dr0)) then
                  dp0=dplocl
                  dr0=drlocl
                  call exchge(isg,2,jsg,ip)
               endif
 100        continue
 120     continue
         if(dp0.le.dpcoal.and.dr0.le.drcoal) IOVER(ISG)=1
 150  continue
c
C2     meson qbar coalesce with all available q:
      do 250 ISG=1,NSG
         if(NJSGS(ISG).ne.2.or.IOVER(ISG).eq.1) goto 250
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
         do 1003 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1003    continue
         call locldr(2,drlocl)
         dr0=drlocl
         dp0=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 220 JSG=1,NSG
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 220
            if(NJSGS(JSG).eq.2) then
               ipmin=1
               ipmax=1
            elseif(NJSGS(JSG).eq.3.and.K2SGS(JSG,1).gt.0) then
               ipmin=1
               ipmax=3
            else
               goto 220
            endif
            do 200 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(2)*pesgs(jsg,ip)
     1              -pxp(2)*pxsgs(jsg,ip)
     2              -pyp(2)*pysgs(jsg,ip)
     3              -pzp(2)*pzsgs(jsg,ip)
     4              -pmp(2)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
               if(dplocl.gt.dpcoal) goto 220
               ftp(1)=ftsgs(jsg,ip)
               gxp(1)=gxsgs(jsg,ip)
               gyp(1)=gysgs(jsg,ip)
               gzp(1)=gzsgs(jsg,ip)
               pxp(1)=pxsgs(jsg,ip)
               pyp(1)=pysgs(jsg,ip)
               pzp(1)=pzsgs(jsg,ip)
               pmp(1)=pmsgs(jsg,ip)
               pep(1)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
               if(drlocl.gt.drcoal) goto 220
c     qbar_isg coalesces with q_jsg:
               if((dp0.gt.dpcoal.or.dr0.gt.drcoal)
     1              .or.(drlocl.lt.dr0)) then
                  dp0=dplocl
                  dr0=drlocl
                  call exchge(isg,1,jsg,ip)
               endif
 200        continue
 220     continue
         if(dp0.le.dpcoal.and.dr0.le.drcoal) IOVER(ISG)=1
 250  continue
c
C3     baryon q (antibaryon qbar) coalesce with all available q (qbar):
      do 350 ISG=1,NSG
         if(NJSGS(ISG).ne.3.or.IOVER(ISG).eq.1) goto 350
         ibaryn=K2SGS(ISG,1)
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
         do 1004 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1004    continue
         call locldr(2,drlocl)
         dr1(2)=drlocl
         dp1(2)=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         ftp(2)=ftsgs(isg,3)
         gxp(2)=gxsgs(isg,3)
         gyp(2)=gysgs(isg,3)
         gzp(2)=gzsgs(isg,3)
         pxp(2)=pxsgs(isg,3)
         pyp(2)=pysgs(isg,3)
         pzp(2)=pzsgs(isg,3)
         pmp(2)=pmsgs(isg,3)
         pep(2)=pesgs(isg,3)
         call locldr(2,drlocl)
         dr1(3)=drlocl
         dp1(3)=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 320 JSG=1,NSG
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 320
            if(NJSGS(JSG).eq.2) then
               if(ibaryn.gt.0) then
                  ipmin=1
               else
                  ipmin=2
               endif
               ipmax=ipmin
            elseif(NJSGS(JSG).eq.3.and.
     1              (ibaryn*K2SGS(JSG,1)).gt.0) then
               ipmin=1
               ipmax=3
            else
               goto 320
            endif
            do 300 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(1)*pesgs(jsg,ip)
     1              -pxp(1)*pxsgs(jsg,ip)
     2              -pyp(1)*pysgs(jsg,ip)
     3              -pzp(1)*pzsgs(jsg,ip)
     4              -pmp(1)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
               if(dplocl.gt.dpcoal) goto 320
               ftp(2)=ftsgs(jsg,ip)
               gxp(2)=gxsgs(jsg,ip)
               gyp(2)=gysgs(jsg,ip)
               gzp(2)=gzsgs(jsg,ip)
               pxp(2)=pxsgs(jsg,ip)
               pyp(2)=pysgs(jsg,ip)
               pzp(2)=pzsgs(jsg,ip)
               pmp(2)=pmsgs(jsg,ip)
               pep(2)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
               if(drlocl.gt.drcoal) goto 320
c     q_isg may coalesce with q_jsg for a baryon:
               ipi=0
               if(dp1(2).gt.dpcoal.or.dr1(2).gt.drcoal) then
                  ipi=2
                  if((dp1(3).gt.dpcoal.or.dr1(3).gt.drcoal)
     1                 .and.dr1(3).gt.dr1(2)) ipi=3
               elseif(dp1(3).gt.dpcoal.or.dr1(3).gt.drcoal) then
                  ipi=3
               elseif(dr1(2).lt.dr1(3)) then
                  if(drlocl.lt.dr1(3)) ipi=3
               elseif(dr1(3).le.dr1(2)) then
                  if(drlocl.lt.dr1(2)) ipi=2
               endif
               if(ipi.ne.0) then
                  dp1(ipi)=dplocl
                  dr1(ipi)=drlocl
                  call exchge(isg,ipi,jsg,ip)
               endif
 300        continue
 320     continue
         if(dp1(2).le.dpcoal.and.dr1(2).le.drcoal
     1        .and.dp1(3).le.dpcoal.and.dr1(3).le.drcoal)
     2        IOVER(ISG)=1
 350  continue
c      
      RETURN
      END

c=======================================================================
      SUBROUTINE exchge(isg,ipi,jsg,ipj)
c
      implicit double precision  (a-h, o-z)
      PARAMETER (MAXSTR=150001)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      SAVE   
c
      k1=K1SGS(isg,ipi)
      k2=K2SGS(isg,ipi)
      px=PXSGS(isg,ipi)
      py=PYSGS(isg,ipi)
      pz=PZSGS(isg,ipi)
      pe=PESGS(isg,ipi)
      pm=PMSGS(isg,ipi)
      gx=GXSGS(isg,ipi)
      gy=GYSGS(isg,ipi)
      gz=GZSGS(isg,ipi)
      ft=FTSGS(isg,ipi)
      K1SGS(isg,ipi)=K1SGS(jsg,ipj)
      K2SGS(isg,ipi)=K2SGS(jsg,ipj)
      PXSGS(isg,ipi)=PXSGS(jsg,ipj)
      PYSGS(isg,ipi)=PYSGS(jsg,ipj)
      PZSGS(isg,ipi)=PZSGS(jsg,ipj)
      PESGS(isg,ipi)=PESGS(jsg,ipj)
      PMSGS(isg,ipi)=PMSGS(jsg,ipj)
      GXSGS(isg,ipi)=GXSGS(jsg,ipj)
      GYSGS(isg,ipi)=GYSGS(jsg,ipj)
      GZSGS(isg,ipi)=GZSGS(jsg,ipj)
      FTSGS(isg,ipi)=FTSGS(jsg,ipj)
      K1SGS(jsg,ipj)=k1
      K2SGS(jsg,ipj)=k2
      PXSGS(jsg,ipj)=px
      PYSGS(jsg,ipj)=py
      PZSGS(jsg,ipj)=pz
      PESGS(jsg,ipj)=pe
      PMSGS(jsg,ipj)=pm
      GXSGS(jsg,ipj)=gx
      GYSGS(jsg,ipj)=gy
      GZSGS(jsg,ipj)=gz
      FTSGS(jsg,ipj)=ft
c
      RETURN
      END

c=======================================================================
      SUBROUTINE locldr(icall,drlocl)
c
      implicit double precision (a-h, o-z)
      dimension ftp0(3),pxp0(3),pyp0(3),pzp0(3),pep0(3)
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
cc      SAVE /loclco/
      common /prtn23/ gxp0(3),gyp0(3),gzp0(3),ft0fom
cc      SAVE /prtn23/
      common /lor/ enenew, pxnew, pynew, pznew
cc      SAVE /lor/
      SAVE   
c     for 2-body kinematics:
      if(icall.eq.2) then
         etot=pep(1)+pep(2)
         bex=(pxp(1)+pxp(2))/etot
         bey=(pyp(1)+pyp(2))/etot
         bez=(pzp(1)+pzp(2))/etot
c     boost the reference frame down by beta to get to the pair rest frame:
         do 1001 j=1,2
            beta2 = bex ** 2 + bey ** 2 + bez ** 2
            gam = 1.d0 / dsqrt(1.d0 - beta2)
            if(beta2.ge.0.9999999999999d0) then
               write(6,*) '4',pxp(1),pxp(2),pyp(1),pyp(2),
     1              pzp(1),pzp(2),pep(1),pep(2),pmp(1),pmp(2),
     2          dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2+pmp(1)**2)/pep(1),
     3          dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2)/pep(1)
               write(6,*) '4a',pxp(1)+pxp(2),pyp(1)+pyp(2),
     1              pzp(1)+pzp(2),etot
               write(6,*) '4b',bex,bey,bez,beta2,gam
            endif
c
            call lorenz(ftp(j),gxp(j),gyp(j),gzp(j),bex,bey,bez)
            gxp0(j)=pxnew
            gyp0(j)=pynew
            gzp0(j)=pznew
            ftp0(j)=enenew
            call lorenz(pep(j),pxp(j),pyp(j),pzp(j),bex,bey,bez)
            pxp0(j)=pxnew
            pyp0(j)=pynew
            pzp0(j)=pznew
            pep0(j)=enenew
 1001    continue
c     
         if(ftp0(1).ge.ftp0(2)) then
            ilate=1
            iearly=2
         else
            ilate=2
            iearly=1
         endif
         ft0fom=ftp0(ilate)
c     
         dt0=ftp0(ilate)-ftp0(iearly)
         gxp0(iearly)=gxp0(iearly)+pxp0(iearly)/pep0(iearly)*dt0
         gyp0(iearly)=gyp0(iearly)+pyp0(iearly)/pep0(iearly)*dt0
         gzp0(iearly)=gzp0(iearly)+pzp0(iearly)/pep0(iearly)*dt0
         drlocl=dsqrt((gxp0(ilate)-gxp0(iearly))**2
     1        +(gyp0(ilate)-gyp0(iearly))**2
     2        +(gzp0(ilate)-gzp0(iearly))**2)
c     for 3-body kinematics, used for baryons formation:
      elseif(icall.eq.3) then
         etot=pep(1)+pep(2)+pep(3)
         bex=(pxp(1)+pxp(2)+pxp(3))/etot
         bey=(pyp(1)+pyp(2)+pyp(3))/etot
         bez=(pzp(1)+pzp(2)+pzp(3))/etot
         beta2 = bex ** 2 + bey ** 2 + bez ** 2
         gam = 1.d0 / dsqrt(1.d0 - beta2)
         if(beta2.ge.0.9999999999999d0) then
            write(6,*) '5',bex,bey,bez,beta2,gam
         endif
c     boost the reference frame down by beta to get to the 3-parton rest frame:
         do 1002 j=1,3
            call lorenz(ftp(j),gxp(j),gyp(j),gzp(j),bex,bey,bez)
            gxp0(j)=pxnew
            gyp0(j)=pynew
            gzp0(j)=pznew
            ftp0(j)=enenew
            call lorenz(pep(j),pxp(j),pyp(j),pzp(j),bex,bey,bez)
            pxp0(j)=pxnew
            pyp0(j)=pynew
            pzp0(j)=pznew
            pep0(j)=enenew
 1002    continue
c     
         if(ftp0(1).gt.ftp0(2)) then
            ilate=1
            if(ftp0(3).gt.ftp0(1)) ilate=3
         else
            ilate=2
            if(ftp0(3).ge.ftp0(2)) ilate=3
         endif
         ft0fom=ftp0(ilate)
c     
         if(ilate.eq.1) then
            imin=2
            imax=3
            istep=1
         elseif(ilate.eq.2) then
            imin=1
            imax=3
            istep=2
         elseif(ilate.eq.3) then
            imin=1
            imax=2
            istep=1
         endif
c     
         do 1003 iearly=imin,imax,istep
            dt0=ftp0(ilate)-ftp0(iearly)
            gxp0(iearly)=gxp0(iearly)+pxp0(iearly)/pep0(iearly)*dt0
            gyp0(iearly)=gyp0(iearly)+pyp0(iearly)/pep0(iearly)*dt0
            gzp0(iearly)=gzp0(iearly)+pzp0(iearly)/pep0(iearly)*dt0
 1003    continue
      endif
c
      RETURN
      END

c=======================================================================
        subroutine hoscar
c
        parameter (MAXSTR=150001,AMN=0.939457,AMP=0.93828)
        character*8 code, reffra, FRAME
        character*25 amptvn
        common/snn/efrm,npart1,npart2
cc      SAVE /snn/
        common /lastt/itimeh,bimp 
cc      SAVE /lastt/
        COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
cc      SAVE /hbt/
        common/oscar1/iap,izp,iat,izt
cc      SAVE /oscar1/
        common/oscar2/FRAME,amptvn
cc      SAVE /oscar2/
        SAVE   
        data nff/0/
c
c       file header
        if(nff.eq.0) then
           write (19, 101) 'OSCAR1997A'
           write (19, 111) 'final_id_p_x'
           code = 'AMPT'
           if(FRAME.eq.'CMS') then
              reffra = 'nncm'
              xmp=(amp*izp+amn*(iap-izp))/iap
              xmt=(amp*izt+amn*(iat-izt))/iat
              ebeam=(efrm**2-xmp**2-xmt**2)/2./xmt
           elseif(FRAME.eq.'LAB') then
              reffra = 'lab'
              ebeam=efrm
           else
              reffra = 'unknown'
              ebeam=0.
           endif
           ntestp = 1
           write (19, 102) code, amptvn, iap, izp, iat, izt,
     &        reffra, ebeam, ntestp
           nff = 1
           ievent = 1
           phi = 0.
           if(FRAME.eq.'CMS') write(19,112) efrm
        endif
c       comment
c       event header
        write (19, 103) ievent, nlast, bimp, phi
c       particles
        do 99 i = 1, nlast
           ene=sqrt(plast(1,i)**2+plast(2,i)**2+plast(3,i)**2
     1          +plast(4,i)**2)
           write (19, 104) i, INVFLV(lblast(i)), plast(1,i),
     1          plast(2,i),plast(3,i),ene,plast(4,i),
     2          xlast(1,i),xlast(2,i),xlast(3,i),xlast(4,i)
 99     continue
        ievent = ievent + 1
 101        format (a10)
 111        format (a12)
 102        format (a4,1x,a20,1x,'(', i3, ',', i3, ')+(', i3, ',', 
     &           i3, ')', 2x, a4, 2x, e10.4, 2x, i8)
 103        format (i10, 2x, i10, 2x, f8.3, 2x, f8.3)
 104        format (i10, 2x, i10, 2x, 9(e12.6, 2x))
 112        format ('# Center-of-mass energy/nucleon-pair is',
     & f12.3,'GeV')
c
        return
        end

c=======================================================================
        subroutine getnp

        PARAMETER (MAXSTR=150001)
        COMMON /HMAIN1/ NATT, EATT, JATT, NT, NP, N0, N01, N10, N11
cc      SAVE /HMAIN1/
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
        COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
cc      SAVE /HPARNT/
        common/snn/efrm,npart1,npart2
cc      SAVE /snn/
        SAVE   

        if(NATT.eq.0) then
           npart1=0
           npart2=0
           return
        endif
c
        PZPROJ=SQRT(HINT1(6)**2-HINT1(8)**2)
        PZTARG=SQRT(HINT1(7)**2-HINT1(9)**2)
        epsiln=0.01
c
        nspec1=0
        nspec2=0
        DO 1000 I = 1, NATT
           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
     1          .and.PATT(I, 1).eq.0.and.PATT(I, 2).eq.0) then
              if(PATT(I, 3).gt.amax1(0.,PZPROJ-epsiln)) then
                 nspec1=nspec1+1
              elseif(PATT(I, 3).lt.(-PZTARG+epsiln)) then
                 nspec2=nspec2+1
              endif
           endif
 1000    CONTINUE
        npart1=IHNT2(1)-nspec1
        npart2=IHNT2(3)-nspec2

        return
        end

c=======================================================================
c     2/18/03 use PYTHIA to decay eta,rho,omega,k*,phi and Delta
        subroutine resdec(i1,nt,nnn,wid,idecay)

        PARAMETER (hbarc=0.19733)
        PARAMETER (AK0=0.498,APICH=0.140,API0=0.135,AN=0.940,ADDM=0.02)
        PARAMETER (MAXSTR=150001, MAXR=1)
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
cc      SAVE /INPUT2/
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
cc      SAVE /LUJETS/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
cc      SAVE /LUDAT2/
        COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
cc      SAVE /LUDAT3/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
cc      SAVE /PA/
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
cc      SAVE /PB/
        COMMON   /PC/EPION(MAXSTR,MAXR)
cc      SAVE /PC/
        COMMON   /PD/LPION(MAXSTR,MAXR)
cc      SAVE /PD/
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
        common/resdcy/NSAV,iksdcy
cc      SAVE /resdcy/
        common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1       px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
        EXTERNAL IARFLV, INVFLV
        COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
cc      SAVE /tdecay/
        COMMON/RNDF77/NSEED
        COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1       dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2       dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
cc      SAVE /RNDF77/
        SAVE   
        irun=idecay
        if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &       .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &       .or.lb1.eq.24.or.(iabs(lb1).ge.6.and.iabs(lb1).le.9) 
     &       .or.iabs(lb1).eq.16) then
           kf=INVFLV(lb1)
        else
           return
        endif
c
        IP=1
c     label as undecayed and the only particle in the record:
        N=1
        K(IP,1)=1
        K(IP,3)=0
        K(IP,4)=0
        K(IP,5)=0
c
        K(IP,2)=kf
        P(IP,1)=px1
        P(IP,2)=py1
        P(IP,3)=pz1
        em1a=em1
c     eta or omega in ART may be below or too close to (pi+pi-pi0) mass, 
c     causing LUDECY error,thus increase their mass ADDM above this thresh,
c     noting that rho (m=0.281) too close to 2pi thrshold fails to decay:
        if((lb1.eq.0.or.lb1.eq.28).and.em1.lt.(2*APICH+API0+ADDM)) then
           em1=2*APICH+API0+ADDM
c     rho
        elseif(lb1.ge.25.and.lb1.le.27.and.em1.lt.(2*APICH+ADDM)) then
           em1=2*APICH+ADDM
c     K*
        elseif(iabs(lb1).eq.30.and.em1.lt.(APICH+AK0+ADDM)) then
           em1=APICH+AK0+ADDM
c     Delta created in ART may be below (n+pich) mass, causing LUDECY error:
        elseif(iabs(lb1).ge.6.and.iabs(lb1).le.9
     1          .and.em1.lt.(APICH+AN+ADDM)) then
           em1=APICH+AN+ADDM
        endif
        if(em1.ge.(em1a+0.01)) write (6,*) 
     1       'Mass increase in resdec():',nt,em1-em1a,lb1
        e1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
        P(IP,4)=e1
        P(IP,5)=em1
clin-5/2008:
        dpdecp=dpertp(i1)
        call ludecy(IP)
c     add decay time to daughter's formation time at the last timestep:
        if(nt.eq.ntmax) then
           tau0=hbarc/wid
           taudcy=tau0*(-1.)*alog(1.-RANART(NSEED))
           ndaut=n-nsav
           if(ndaut.le.1) then
              write(10,*) 'note: ndaut(<1)=',ndaut
              write(89,*) 'note: ndaut(<1)=',ndaut
              call lulist(2)
              stop
            endif
c     lorentz boost:
           taudcy=taudcy*e1/em1
           tfnl=tfnl+taudcy
           xfnl=xfnl+px1/e1*taudcy
           yfnl=yfnl+py1/e1*taudcy
           zfnl=zfnl+pz1/e1*taudcy
c     at the last timestep, assign rho, K0S or eta (decay daughter) 
c     to lb(i1) only (not to lpion) in order to decay them again:
           if(n.ge.(nsav+2)) then
              do 1001 idau=nsav+2,n
                 kdaut=K(idau,2)
                 if(kdaut.eq.221.or.kdaut.eq.113
     1                .or.kdaut.eq.213.or.kdaut.eq.-213
     2                .or.kdaut.eq.310) then
c     switch idau and i1(nsav+1):
                    ksave=kdaut
                    pxsave=p(idau,1)
                    pysave=p(idau,2)
                    pzsave=p(idau,3)
                    esave=p(idau,4)
                    xmsave=p(idau,5)
                    K(idau,2)=K(nsav+1,2)
                    p(idau,1)=p(nsav+1,1)
                    p(idau,2)=p(nsav+1,2)
                    p(idau,3)=p(nsav+1,3)
                    p(idau,4)=p(nsav+1,4)
                    p(idau,5)=p(nsav+1,5)
                    K(nsav+1,2)=ksave
                    p(nsav+1,1)=pxsave
                    p(nsav+1,2)=pysave
                    p(nsav+1,3)=pzsave
                    p(nsav+1,4)=esave
                    p(nsav+1,5)=xmsave
c     note: phi decay may produce rho, K0s or eta, N*(1535) decay may produce 
c     eta, but only one daughter may be rho, K0s or eta:
                    goto 111
                 endif
 1001         continue
           endif
 111       continue
c     
           enet=0.
           do 1002 idau=nsav+1,n
              enet=enet+p(idau,4)
 1002      continue
c           if(abs(enet-e1).gt.0.02) 
c     1          write(93,*) 'resdec(): nt=',nt,enet-e1,lb1
        endif

 200    format(a20,3(1x,i6))
 210    format(i6,5(1x,f8.3))
 220    format(a2,i5,5(1x,f8.3))

        do 1003 idau=nsav+1,n
           kdaut=K(idau,2)
           lbdaut=IARFLV(kdaut)
c     K0S and K0L are named K+/K- during hadron cascade, and only 
c     at the last timestep they keep their real LB # before output;
c     K0/K0bar (from K* decay) converted to K0S and K0L at the last timestep:
           if(nt.eq.ntmax.and.(kdaut.eq.130.or.kdaut.eq.310
     1          .or.iabs(kdaut).eq.311)) then
              if(kdaut.eq.130) then
                 lbdaut=22
              elseif(kdaut.eq.310) then
                 lbdaut=24
              elseif(iabs(kdaut).eq.311) then
                 if(RANART(NSEED).lt.0.5) then
                    lbdaut=22
                 else
                    lbdaut=24
                 endif
              endif
           endif
c
           if(idau.eq.(nsav+1)) then
              LB(i1)=lbdaut
              E(i1)=p(idau,5)
              px1n=p(idau,1)
              py1n=p(idau,2)
              pz1n=p(idau,3)
clin-5/2008:
              dp1n=dpdecp
           else
              nnn=nnn+1
              LPION(NNN,IRUN)=lbdaut
              EPION(NNN,IRUN)=p(idau,5)
              PPION(1,NNN,IRUN)=p(idau,1)
              PPION(2,NNN,IRUN)=p(idau,2)
              PPION(3,NNN,IRUN)=p(idau,3)
              RPION(1,NNN,IRUN)=xfnl
              RPION(2,NNN,IRUN)=yfnl
              RPION(3,NNN,IRUN)=zfnl
              tfdpi(NNN,IRUN)=tfnl
clin-5/2008:
              dppion(NNN,IRUN)=dpdecp
           endif
 1003   continue
 230    format(a2,i5,5(1x,e8.2))
        return
        end

c=======================================================================
        subroutine inidcy

        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
cc      SAVE /LUJETS/
        common/resdcy/NSAV,iksdcy
cc      SAVE /resdcy/
        SAVE   
        N=1
        NSAV=N
        return
        end

c=======================================================================
clin-6/06/02 local parton freezeout motivated from critical density:
        subroutine local(t)
c
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        PARAMETER (r0=1d0)
        COMMON /para1/ mul
cc      SAVE /para1/
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
cc      SAVE /prec2/
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
cc      SAVE /prec4/
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
cc      SAVE /prec5/
        common /coal/dpcoal,drcoal,ecritl
cc      SAVE /coal/
        SAVE   
c
      do 1001 it=1,301
         if(t.ge.tfrz(it).and.t.lt.tfrz(it+1)) then
            if(it.eq.itlast) then
               return
            else
               itlast=it
               goto 50
            endif
         endif
 1001 continue
      write(1,*) 'local time out of range in LOCAL, stop',t,it
      stop
 50   continue
c
      do 200 ip=1,mul
c     skip partons which have frozen out:
         if(ifrz(ip).eq.1) goto 200
         if(it.eq.301) then
c     freezeout all the left partons beyond the time of 3000 fm:
            etcrit=1d6
            goto 150
         else
c     freezeout when transverse energy density < etcrit:
            etcrit=(ecritl*2d0/3d0)
         endif
c     skip partons which have not yet formed:
         if(t.lt.FT5(ip)) goto 200
         rap0=rap(ip)
         eta0=eta(ip)
         x0=GX5(ip)+vx(ip)*(t-FT5(ip))
         y0=GY5(ip)+vy(ip)*(t-FT5(ip))
         detdy=0d0
         do 100 itest=1,mul
c     skip self and partons which have not yet formed:
            if(itest.eq.ip.or.t.lt.FT5(itest)) goto 100
            ettest=eta(itest)
            xtest=GX5(itest)+vx(itest)*(t-FT5(itest))
            ytest=GY5(itest)+vy(itest)*(t-FT5(itest))
            drt=sqrt((xtest-x0)**2+(ytest-y0)**2)
c     count partons within drt<1 and -1<(eta-eta0)<1:
            if(dabs(ettest-eta0).le.1d0.and.drt.le.r0) 
     1           detdy=detdy+dsqrt(PX5(itest)**2+PY5(itest)**2
     2           +XMASS5(itest)**2)*0.5d0
 100     continue
         detdy=detdy*(dcosh(eta0)**2)/(t*3.1416d0*r0**2*dcosh(rap0))
c     when density is below critical density for phase transition, freeze out:
 150     if(detdy.le.etcrit) then
            ifrz(ip)=1
            idfrz(ip)=ITYP5(ip)
            pxfrz(ip)=PX5(ip)
            pyfrz(ip)=PY5(ip)
            pzfrz(ip)=PZ5(ip)
            efrz(ip)=E5(ip)
            xmfrz(ip)=XMASS5(ip)
            if(t.gt.FT5(ip)) then
               gxfrz(ip)=x0
               gyfrz(ip)=y0
               gzfrz(ip)=GZ5(ip)+vz(ip)*(t-FT5(ip))
               ftfrz(ip)=t
            else
c     if this freezeout time < formation time, use formation time & positions.
c     This ensures the recovery of default hadron when e_crit=infty:
               gxfrz(ip)=GX5(ip)
               gyfrz(ip)=GY5(ip)
               gzfrz(ip)=GZ5(ip)
               ftfrz(ip)=FT5(ip)
            endif
         endif
 200  continue
c
        return
        end

c=======================================================================
clin-6/06/02 initialization for local parton freezeout
        subroutine inifrz
c
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
        SAVE   
c
c     for freezeout time 0-10fm, use interval of 0.1fm; 
c     for 10-100fm, use interval of 1fm; 
c     for 100-1000fm, use interval of 10fm; 
c     for 1000-3000fm, use interval of 100fm: 
        step1=0.1d0
        step2=1d0
        step3=10d0
        step4=100d0
c     
        do 1001 it=1,101
           tfrz(it)=0d0+dble(it-1)*step1
 1001 continue
        do 1002 it=102,191
           tfrz(it)=10d0+dble(it-101)*step2
 1002   continue
        do 1003 it=192,281
           tfrz(it)=100d0+dble(it-191)*step3
 1003   continue
        do 1004 it=282,301
           tfrz(it)=1000d0+dble(it-281)*step4
 1004   continue
        tfrz(302)=tlarge
c
        return
        end
