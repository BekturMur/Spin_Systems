        subroutine chebsdPDDG(ioutstat, ! out status
     +                    nfiles,
     +      aspT,bspT,JxT,JyT,JzT,nconnT,HxT,HyT,HzT,
     +      LsysT,LbathT,LtotT,iRandT,nstepsT,tauT,init_stateT,
     +      autonormT,EpsFT,imflagT,ispecflagT, 
     +      iNtypeT,nooutputT, imarkT,
     +      vecResOut,iptVecResOut,mrkResOut)

        implicit none
        include "chsdpar.h"

c in/out variables

        integer ioutstat 
        integer nfiles
        double precision JxT(maxNconn,maxLen), JyT(maxNconn,maxLen)
        double precision JzT(maxNconn,maxLen)
        double precision HxT(maxLtot,maxLen), HyT(maxLtot,maxLen)
        double precision HzT(maxLtot,maxLen), tauT(maxLen)
        integer LsysT(maxLen),LbathT(maxLen),LtotT(maxLen)
        integer nconnT(maxLen)
        integer aspT(maxNconn,maxLen),bspT(maxNconn,maxLen)
        integer nstepsT(maxLen),iRandT(maxLen)
        integer init_stateT(maxLen),autonormT(maxLen),imflagT(maxLen)
        double precision EpsFT(maxLen)
        integer ispecflagT(maxLen)
        integer imarkT(2,maxLen), nooutputT(maxLen)

        double precision vecResOut(maxLenVecResOut)
        integer iptVecResOut
        integer mrkResOut(maxLenVecResOut)

        integer iNtypeT(maxLtot,maxLen)

c internal variables

        double precision psiR(0:maxStat-1), psiI(0:maxStat-1)
        double precision Jx(maxNconn),Jy(maxNconn), Jz(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot), tau
        integer nconn, asp(maxNconn),bsp(maxNconn), Lsys, Lbath, Ltot
        integer nsteps,iRand,init_state,autonorm,imflag
        double precision timeGl, alp,emax, EpsF
        integer firstime, ispecflag, nooutput
        integer imark(2,maxLen)
        integer iNtype(maxLtot)
        integer ifile, i,k, ios
        integer irnd
        double precision xrnd

        double precision phiR(0:maxStat-1), phiI(0:maxStat-1)
        double precision dpsiR(0:maxStat-1), dpsiI(0:maxStat-1)

c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
        double precision timeGlz, pzR(0:maxStat-1),pzI(0:maxStat-1),xx
        double precision xnR, xnI, echoXr,echoXi,echoYr,echoYi, ddot
        double precision echoZr, echoZi
        double precision nrmz,nrmpsi
        integer nstates, lenrec
        external ddot,addHSx
c END

c executable statements

c        integer lr
c        open(chlog,file=logfnam,access='append')
c        write(chlog,'(a,i4)') 'in chebsd, lr=',lr
c        close(chlog)

        do i=0,maxStat-1
          psiR(i) = 0.0D0
          psiI(i) = 0.0D0
        enddo
        firstime = 1
        timeGl = 0.0
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
        Lsys=LsysT(1)
        Lbath=LbathT(1)
        Ltot=LtotT(1)
        do k=1,Ltot
          Hx(k) = HxT(k,1)
          Hy(k) = HyT(k,1)
          Hz(k) = HzT(k,1)
        enddo
        nconn = nconnT(1)
        do k=1,nconn
          asp(k) = aspT(k,1)
          bsp(k) = bspT(k,1)
          Jx(k) = JxT(k,1)
          Jy(k) = JyT(k,1)
          Jz(k) = JzT(k,1)
        enddo
        iRand = iRandT(1)
        autonorm = autonormT(1)
        nooutput = 1
        init_state = -3 ! RANDOM

        call chinitPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot,iRand,init_state,ioutstat,
     +    psiR,psiI,dpsiR,dpsiI)
        if (ioutstat.ne.0) return

        nstates=2**Ltot
        do i=0,nstates-1
          pzR(i)=0.0D0
          pzI(i)=0.0D0
        enddo
        xx=1.0D0
        do i=1,Ltot
          call addHSx(xx,psiR,psiI,i,nstates,pzR,pzI)
        enddo
        timeGlz=0.0D0
c END
        do ifile=1,maxLen
          imark(1,ifile) = imarkT(1,ifile)
          imark(2,ifile) = imarkT(2,ifile)
        enddo

        iptVecResOut = 0
        ifile = 1
        do while (ifile.le.nfiles)
c          open(chlog,file=logfnam,access='append')
c          write(chlog,'(a)') ' '
c          write(chlog,'(a,i5)')'---> Executing data file #', ifile
c          close(chlog)

c          execute file

          Lsys = LsysT(ifile)
          Lbath = LbathT(ifile)
          Ltot = LtotT(ifile)
          do k=1,Ltot
            Hx(k) = HxT(k,ifile)
            Hy(k) = HyT(k,ifile)
            Hz(k) = HzT(k,ifile)
          enddo
          nconn = nconnT(ifile)
          do k=1,nconn
            asp(k) = aspT(k,ifile)
            bsp(k) = bspT(k,ifile)
            Jx(k) = JxT(k,ifile)
            Jy(k) = JyT(k,ifile)
            Jz(k) = JzT(k,ifile)
          enddo
          iRand = iRandT(ifile)
          nsteps = nstepsT(ifile)
          tau = tauT(ifile)
          init_state = init_stateT(ifile)
          autonorm = autonormT(ifile)
          EpsF = EpsFT(ifile)
          imflag = imflagT(ifile)
          ispecflag = ispecflagT(ifile)
          nooutput = nooutputT(ifile)
c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
          nooutput=1
c END
          do k=1,Ltot
            iNtype(k) = iNtypeT(k,ifile)
          enddo

          do irnd=0,iRand
c            call random_number(xrnd)
            xrnd=rand()
          enddo

c THIS IS MOD FOR 2D NMR MODELING
c BEGIN
          if (imflag.eq.1) then
            call chimstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGl,ioutstat,
     +        psiR,psiI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else if (imflag.eq.0) then
            call chstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGl,ioutstat,
     +        psiR,psiI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else 
           call chspecstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGl,ioutstat,
     +        psiR,psiI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        ispecflag,nooutput, iNtype,
     +        vecResOut,iptVecResOut,mrkResOut)
          endif
          if (ioutstat.ne.0) then
c            open(chlog,file=logfnam,access='append')
c            write(chlog,*) 'Runtime error in the file #',ifile
c            close(chlog)
            ioutstat = ioutstat - 10000*ifile
            return
          endif

          if (imflag.eq.1) then
            call chimstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGlz,ioutstat,
     +        pzR,pzI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else if (imflag.eq.0) then
            call chstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGlz,ioutstat,
     +        pzR,pzI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else 
           call chspecstepsPDDGnmr(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +        Lsys,Lbath,Ltot,nsteps,tau,timeGlz,ioutstat,
     +        pzR,pzI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        ispecflag,nooutput, iNtype,
     +        vecResOut,iptVecResOut,mrkResOut)
          endif
          if (ioutstat.ne.0) then
c            open(chlog,file=logfnam,access='append')
c            write(chlog,*) 'Runtime error in the file #',ifile
c            close(chlog)
            ioutstat = ioutstat - 50000*ifile
            return
          endif

          lenrec=7
          if (lenrec+iptVecResOut.gt.maxLenVecResOut) then
            ioutstat=-456 - 70000*ifile
            return
          endif
          mrkResOut(iptVecResOut+1)=lenrec
          vecResOut(iptVecResOut+1)=timeGl

          do i=0,nstates-1
            phiR(i)=0.0D0
            phiI(i)=0.0D0
          enddo
          xx=1.0D0
          do i=1,Ltot
            call addHSx(xx,psiR,psiI,i,nstates,phiR,phiI)
          enddo
          xnR = ddot(nstates,phiR,1,pzR,1)
          xnI = ddot(nstates,phiI,1,pzI,1)
          echoXr=xnR + xnI
          xnR = ddot(nstates,phiR,1,pzI,1)
          xnI = ddot(nstates,phiI,1,pzR,1)
          echoXi=xnI - xnR
          vecResOut(iptVecResOut+2)=echoXr
          vecResOut(iptVecResOut+3)=echoXi

          do i=0,nstates-1
            phiR(i)=0.0D0
            phiI(i)=0.0D0
          enddo
          xx=1.0D0
          do i=1,Ltot
            call addHSy(xx,psiR,psiI,i,nstates,phiR,phiI)
          enddo

          xnR = ddot(nstates,phiR,1,pzR,1)
          xnI = ddot(nstates,phiI,1,pzI,1)
          echoYr=xnR + xnI
          xnR = ddot(nstates,phiR,1,pzI,1)
          xnI = ddot(nstates,phiI,1,pzR,1)
          echoYi=xnI - xnR
c          vecResOut(iptVecResOut+2)=echoYr
c          vecResOut(iptVecResOut+3)=echoYi








c now calculate norms, use same notations
          xnR = ddot(nstates,pzR,1,pzR,1)
          xnI = ddot(nstates,pzI,1,pzI,1)
          nrmz=(xnR + xnI)/LtotT(1)
          xnR = ddot(nstates,psiR,1,psiR,1)
          xnI = ddot(nstates,psiI,1,psiI,1)
          nrmpsi=xnI + xnR
          vecResOut(iptVecResOut+4)=nrmz
          vecResOut(iptVecResOut+5)=nrmpsi
c except for nrmz,nrmpsi


c Calculation mean value of S_{z} (S_{z}S_{z})

          do i=0,nstates-1
            phiR(i)=0.0D0
            phiI(i)=0.0D0
          enddo
          xx=1.0D0
          do i=1,Ltot
            call addHSz(xx,psiR,psiI,i,nstates,phiR,phiI)
          enddo



          xnR = ddot(nstates,phiR,1,pzR,1)
          xnI = ddot(nstates,phiI,1,pzI,1)
          echoZr=xnR + xnI
          xnR = ddot(nstates,phiR,1,pzI,1)
          xnI = ddot(nstates,phiI,1,pzR,1)
          echoZi=xnI - xnR
c          vecResOut(iptVecResOut+2)=echoZr
c          vecResOut(iptVecResOut+3)=echoZi





          iptVecResOut = iptVecResOut+lenrec
c END
          firstime = 0 ! not the first time already

          if (imark(2,ifile).gt.1) then ! mark of the endcycle
            imark(2,ifile) = imark(2,ifile) - 1
            ifile = imark(1,ifile)
          else 
            ifile = ifile + 1
          endif
        enddo ! the ifile loop

c            open(chlog,file=logfnam,access='append')
c            write(chlog,'(a)')
c     +        '----------------------------------------------'
c            write(chlog,'(2a)') 'Exhausted input file ',
c     +        infilnam(:len_trim(infilnam))
c            write(chlog,'(a)')
c     +        'END OF Chebyshev-based quantum spin dynamics'
c            write(chlog,'(a)')
c     +        '----------------------------------------------'
c            close(chlog)

        return
        end
