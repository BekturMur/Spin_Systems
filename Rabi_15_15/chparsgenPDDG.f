        subroutine chparsgenPDDG(ioutstat, ! out status
     +                    infilnam, ! name of the input file
     +                    nfiles,
     +      aspT,bspT,JxT,JyT,JzT,nconnT,HxT,HyT,HzT,
     +      LsysT,LbathT,LtotT,iRandT,nstepsT,tauT,init_stateT,
     +      autonormT,EpsFT,imflagT,ispecflagT,
     +      iNtypeT,nooutputT, imarkT)

        implicit none
        include "chsdpar.h"

c in/out variables

        integer ioutstat, nfiles 
        character*32 infilnam

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
        integer iNtypeT(maxLtot,maxLen)

c internal variables

        double precision Jx(maxNconn), Jy(maxNconn)
        double precision Jz(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot)
        double precision Hz(maxLtot), tau
        integer Lsys,Lbath,Ltot
        integer nconn
        integer asp(maxNconn),bsp(maxNconn)
        integer nsteps,iRand
        integer init_state,autonorm,imflag
        double precision EpsF
        integer ispecflag
        integer goOn, ifile, k, firstime
        character*32 curfilnam
        character*128 infline

        integer ios
        logical exOut

        integer nooutput, iNtype(maxLtot)
        integer kBegCycle,nRepCycle,incycle

c executable statements

        inquire(file=infilnam(:len_trim(infilnam)), exist=exOut)
        if (.not.exOut) then
          ioutstat = -101
          open(chlog,file=logfnam,access='append')
          write(chlog,*) 'No input file ', 
     +      infilnam(:len_trim(infilnam))
          close(chlog)
          return
        endif
        open(chin,file=infilnam(:len_trim(infilnam)),iostat=ios)
        if (ios.ne.0) then
          ioutstat=-101
          open(chlog,file=logfnam,access='append')
          write(chlog,*) 'Cannot open input file ',
     +      infilnam(:len_trim(infilnam))
          close(chlog)
          return
        endif

        goOn = -1 ! first time
        firstime = 1
        ifile=0

        incycle=0

        do while (goOn.ne.0) 
          read(chin,'(a)',iostat=ios) infline
          if (ios.ne.0) then
            ioutstat = 0
            goOn = 0

            nfiles = ifile

            close(chin)
            open(chlog,file=logfnam,access='append')
            write(chlog,'(a)')
     +        '----------------------------------------------'
            write(chlog,'(2a)') 'Exhausted input file ',
     +        infilnam(:len_trim(infilnam))
            write(chlog,'(a,i5)') 'Total number of files:', nfiles
            write(chlog,'(a)')
     +        '----------------------------------------------'
            close(chlog)
            return
          endif

          if (index(infline(:len_trim(infline)),'@CYCLE').ne.0) then
            if (incycle.ne.0) then
              ioutstat=-987
              if (goOn.ne.0) close(chin)
              open(chlog,file=logfnam,access='append')
              write(chlog,'(a)') 'ERROR: nested cycle'
              close(chlog)
              return
            else
              open(chlog,file=logfnam,access='append')
              write(chlog,*) ' '
              write(chlog,*) '============================'
              write(chlog,'(a$)') 'START CYCLE '
              close(chlog)
              k=index(infline,'CYCLE')+len('CYCLE')
              read(infline(k:),*,iostat=ios) nRepCycle
              if (ios.ne.0) then
                nRepCycle = 1
                open(chlog,file=logfnam,access='append')
                write(chlog,*) ' '
                write(chlog,*) '!!!'
                write(chlog,*)
     +            'WARNING: NUMBER OF CYCLE REPETITIONS IS NOT SHOWN'
                write(chlog,*) 'ASSUME IT IS 1'
                write(chlog,*) '!!!'
                close(chlog)
              endif
              open(chlog,file=logfnam,access='append')
              write(chlog,*) nRepCycle, 'REPETITIONS'
              close(chlog)
              incycle=1
              kBegCycle=ifile+1
              cycle
            endif
          else if (index(infline(:len_trim(infline)),'@ENDCYCLE').ne.0) 
     +      then
            if (incycle.eq.0) then
              ioutstat=-789
              if (goOn.ne.0) close(chin)
              open(chlog,file=logfnam,access='append')
              write(chlog,'(a)')'ERROR: ENDCYCLE without matching CYCLE'
              close(chlog)
              return
            else
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'END CYCLE'
              write(chlog,*) '============================'
              close(chlog)
              incycle=0
              if (ifile.ge.kBegCycle) then ! if cycle body non-empty
                imarkT(1,ifile)=kBegCycle
                imarkT(2,ifile)=nRepCycle
              endif
              cycle
            endif
          endif

          goOn = index(infline(:len_trim(infline)),'@')
          if (goOn.eq.0) then
            close(chin)
            curfilnam = infilnam
          else
            goOn = 1
            curfilnam = infline((goOn+1):LENfnam)
          endif

          open(chlog,file=logfnam,access='append')
          write(chlog,'(a)') ' '
          write(chlog,'(2a)')'---> Reading data file ',
     +      curfilnam(:len_trim(curfilnam))
          write(chlog,'(a,i5)') 'File number ',ifile+1
          close(chlog)

          call chparsFPDDG(curfilnam,asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +      Lsys,Lbath,Ltot,iRand,nsteps,tau,init_state,firstime,
     +      autonorm,ioutstat,EpsF,imflag,ispecflag,nooutput,iNtype)
          if (ioutstat.ne.0) then
            if (goOn.ne.0) close(chin)
            return
          endif

          ifile = ifile+1
          if (ifile.gt.maxLen) then
            if (goOn.ne.0) close(chin)
            open(chlog,file=logfnam,access='append')
            write(chlog,'(a)') ' '
            write(chlog,'(a)') 'ERROR:'
            write(chlog,'(a)') 'Number of data files larger than maxLen'
            close(chlog)
            ioutstat = -1701
            return
          endif

          LsysT(ifile) = Lsys
          LbathT(ifile) = Lbath
          LtotT(ifile) = Ltot
          do k=1,Ltot
            HxT(k,ifile) = Hx(k)
            HyT(k,ifile) = Hy(k)
            HzT(k,ifile) = Hz(k)
          enddo
          nconnT(ifile) = nconn
          do k=1,nconn
            aspT(k,ifile) = asp(k)
            bspT(k,ifile) = bsp(k)
            JxT(k,ifile) = Jx(k)
            JyT(k,ifile) = Jy(k)
            JzT(k,ifile) = Jz(k)
          enddo
          iRandT(ifile) = iRand
          nstepsT(ifile) = nsteps
          tauT(ifile) = tau
          init_stateT(ifile) = init_state
          autonormT(ifile) = autonorm
          EpsFT(ifile) = EpsF
          imflagT(ifile) = imflag
          ispecflagT(ifile) = ispecflag
          nooutputT(ifile) = nooutput

          imarkT(1,ifile) = ifile+1 ! normal setting
          imarkT(2,ifile) = 1       ! if ENDCYCLE encountered, will be modified

          do k=1,Ltot
            iNtypeT(k,ifile)=iNtype(k)
          enddo

          firstime = 0 ! not the first time already

          call chlogparPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,Lsys,
     +      Lbath,Ltot,iRand,nsteps,tau,init_state,autonorm,imflag,
     +      ispecflag,nooutput,iNtype)

        enddo ! the while(goOn) loop

        if (incycle.gt.0) then
          open(chlog,file=logfnam,access='append')
          write(chlog,'(a)')'ERROR: CYCLE without matching ENDCYCLE'
          close(chlog)
          ioutstat = -879 
        endif

        end
