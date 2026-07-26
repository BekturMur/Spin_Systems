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
        logical exOut

        double precision phiR(0:maxStat-1), phiI(0:maxStat-1)
        double precision dpsiR(0:maxStat-1), dpsiI(0:maxStat-1)

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

          do k=1,Ltot
            iNtype(k) = iNtypeT(k,ifile)
          enddo

          do irnd=0,iRand
c            call random_number(xrnd)
            xrnd=rand()
          enddo

          if (imflag.eq.1) then
            call chimstepsPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,Lsys,
     +        Lbath,Ltot,iRand,nsteps,tau,init_state,timeGl,ioutstat,
     +        psiR,psiI,autonorm,firstime,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else if (imflag.eq.0) then
            call chstepsPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,Lsys,
     +        Lbath,Ltot,iRand,nsteps,tau,init_state,timeGl,ioutstat,
     +        psiR,psiI,autonorm,firstime,EpsF,phiR,phiI,dpsiR,dpsiI,
     +        nooutput,
     +        vecResOut,iptVecResOut,mrkResOut)
          else 
            call chspecstepsPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,Lsys,
     +        Lbath,Ltot,iRand,nsteps,tau,init_state,timeGl,ioutstat,
     +        psiR,psiI,autonorm,firstime,EpsF,phiR,phiI,dpsiR,dpsiI,
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
