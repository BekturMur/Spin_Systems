       subroutine chparsFPDDG(curfilnam,asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot,iRand,nsteps,tau,init_state,firstime,
     +    autonorm,ioutstat,EpsF,imflag,ispecflag,nooutput,iNtype)

        implicit none
        include "chsdpar.h"

c in/out variables

        integer ioutstat,firstime, autonorm, imflag
        character*32 curfilnam
        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn), EpsF
        integer asp(maxNconn), bsp(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot), tau
        integer Lsys,Lbath,Ltot,nconn, nsteps, iRand, init_state
        integer ispecflag, nooutput
        integer iNtype(maxLtot)

c internal variables

        integer i,j,ki, iconn, nlin, ntauint
        integer errsteps,errtau
        double precision x,y,z
        integer ios, istline
        character*128 curline
        logical eX

        double precision hx0,hx1,hy0,hy1,hz0,hz1,ox1,oy1,oz1
        double precision phix,phiy,phiz

c executable statements

        inquire(file=curfilnam(:len_trim(curfilnam)), exist=eX)
        if (.not.eX) then
          ioutstat = -103
          open(chlog,file=logfnam,access='append')
          write(chlog,*) 'No data file ',curfilnam(:len_trim(curfilnam))
          close(chlog)
          return
        endif
        open(chcur, file=curfilnam(:len_trim(curfilnam)), iostat=ios)
        if (ios.ne.0) then
          ioutstat = -103
          open(chlog,file=logfnam,access='append')
          write(chlog,*) 'Error opening data file ',
     +      curfilnam(:len_trim(curfilnam))
          close(chlog)
          return
        endif

c first, re-initialize all parameters

        do i=1,maxLtot
          Hx(i)=0.0
          Hy(i)=0.0
          Hz(i)=0.0
        enddo
        do j=1,maxNconn
          Jx(j) = 0.0
          Jy(j) = 0.0
          Jz(j) = 0.0
          asp(j) = 0
          bsp(j) = 0
        enddo
        nconn = 0
        iRand = -1 ! default value
        nsteps = 0
        errsteps = 1 ! error a priori
        tau = 0.0
        errtau = 1 ! error a priori
        imflag = 0 ! a priori, real-time propagation
        init_state = -10 ! not first time, use previous result
        if (firstime.ne.0) then
          init_state = -1 ! first time, default: SPIN UP
        endif
        EpsF = Eps ! default: EpsF = Eps
        Lsys = 0
        Lbath = 0
        Ltot = 0
        autonorm = 0 ! default: no autonorm
        ispecflag = 0
        nooutput = 0 ! default: print output

        do j=1,maxLtot
          iNtype(j) = -1
        enddo

c reading actual values

        read(chcur,'(a)',iostat=ios) curline
        do while (index(curline,'END').eq.0.and.ios.eq.0) 

          if (index(curline,'HEISENBERG INTERACTION').ne.0) then

            iconn = 0
            read(chcur,'(a)',iostat=ios) curline
            do while (index(curline,'END').eq.0.and.ios.eq.0) 
              read(curline,*,iostat=istline) i,j,x,y,z
              if (istline.ne.0) then
                ioutstat = -11
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'i,j,x,y,z expected in data file ',
     +            curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              if (i.eq.j) then
                ioutstat = -12
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'i=j in data file ',
     +            curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              if (i.gt.maxLtot.or.j.gt.maxLtot.or.i.le.0.or.j.le.0) then
                ioutstat = -13
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'i or j >maxLtot or <0 in data file ',
     +            curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              if (nconn+1.gt.maxNconn) then ! is it safe to make nconn=nconn+1?
                ioutstat = -14
                nconn=nconn+1
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'nconn >maxNconn in data file ',
     +            curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              if (i.gt.j) then
                ki = i
                i = j
                j = ki
              endif
              ki=1 ! check if J_{i,j} already exists
              do while ((asp(ki).ne.i.or.bsp(ki).ne.j).and.ki.le.iconn)
                ki=ki+1
              enddo
              if (ki.le.iconn) then ! J_{ij} exists, update it
                Jx(ki) = x
                Jy(ki) = y
                Jz(ki) = z
              else
                iconn=iconn+1
                nconn=max(nconn,iconn) 
                asp(iconn) = i
                bsp(iconn) = j
                Jx(iconn) = x
                Jy(iconn) = y
                Jz(iconn) = z
              endif
              Ltot = max(Ltot,i,j)
              read(chcur,'(a)',iostat=ios) curline
            enddo ! reading Heisenberg interactions

          else if (index(curline,'SPIN :').ne.0) then

            nlin = index(curline,'SPIN :')+len('SPIN :')
            read(curline(nlin:),*,iostat=istline) i
            if (istline.ne.0) then
              ioutstat = -21
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'spin index expected in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            if (i.le.0.or.i.gt.maxLtot) then
              ioutstat = -22
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'spin index <0 or >maxLtot in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0) then
              ioutstat = -23
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Error reading SPIN ',i,' in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(curline,*,iostat=istline) hx0,hx1,ox1,hy0,hy1,oy1,
     +        hz0,hz1,oz1,phix,phiy,phiz
            if (istline.ne.0) then
              ioutstat = -24
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'hx,y,z etc. expected in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            Hx(i) = hx0
            Hy(i) = hy0
            Hz(i) = hz0
            Ltot = max(Ltot,i)
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -25
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END for spin ',i,' in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif

          else if (index(curline,'SUB-SYSTEM SIZE').ne.0) then

            read(chcur,*,iostat=ios) Lsys
            if (ios.ne.0) then
              ioutstat = -31
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing L_system in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            Ltot = max(Ltot,Lsys+Lbath)
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -32
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of L_system in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
 
          else if (index(curline,'SPINS TO BE TRACED OUT').ne.0) then

            read(chcur,*,iostat=ios) Lbath
            if (ios.ne.0) then
              ioutstat = -41
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing L_bath in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            Ltot = max(Ltot,Lsys+Lbath)
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -42
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of L_bath in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
 
          else if (index(curline,'START RANFX FROM').ne.0) then

            read(chcur,*,iostat=ios) iRand
            if (ios.ne.0) then
              ioutstat = -51
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing iRand in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -52
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of iRand in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
 
          else if (index(curline,'NUMBER OF FIELD STEPS').ne.0) then

            errsteps = 0
            read(chcur,*,iostat=ios) nsteps
            if (ios.ne.0) then
              ioutstat = -61
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing Nsteps in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -62
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of Nsteps in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
 
          else if (index(curline,'TIME STEP ;').ne.0) then

            errtau = 0
            read(chcur,*,iostat=ios) tau 
            if (ios.ne.0) then
              ioutstat = -71
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing tau in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -72
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of tau in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif

          else if (index(curline,'IMAGINARY TIME PROPAGATION').ne.0)then
            imflag=1

          else if (index(curline,'SPECIAL PROPAGATION').ne.0)then 
            imflag=2
            if (index(curline,':').ne.0) then
              nlin=index(curline,':')+len(':')
              read(curline(nlin:),*,iostat=istline) ispecflag
              if (istline.ne.0) then
                ioutstat = -9002
                open(chlog,file=logfnam,access='append')
                write(chlog,*)'ispecflag expected after colon',
     +            'in data file', curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
            endif

          else if (index(curline,'NUMBER OF INTERMEDIATE TIME STEPS')
     +      .ne.0) then

            read(chcur,*,iostat=ios) ntauint
            if (ios.ne.0) then
              ioutstat = -73
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing ntauinterm in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            tau = tau * ntauint
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -74
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of ntauinterm in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif

          else if (index(curline,'MINIMUM EXPANSION COEFFICIENT')
     +      .ne.0) then

            read(chcur,*,iostat=ios) EpsF
            if (ios.ne.0) then
              ioutstat = -75
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing EpsF in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -76
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of EpsF in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif

          else if (index(curline,'HNORM AUTO').ne.0) then 

            read(chcur,*,iostat=ios) autonorm
            if (ios.ne.0) then
              ioutstat = -91
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing autonorm in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
            if (autonorm.eq.-1.or.autonorm.eq.0) then
              autonorm= 0
            else
              autonorm= 1
            endif 
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -92
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Missing END of autonorm in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif

          else if (index(curline,'NO OUTPUT').ne.0) then
            nooutput= 1

          else if (index(curline,'SPIN TYPES').ne.0) then
            read(chcur,'(a)') curline
c            read(curline,*,iostat=ios) (iNtype(j),j=1,maxLtot)
            read(curline,*,end=11) iNtype
11          continue
            if (ios.ne.0) then
              ioutstat = -345
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Error reading SPIN TYPES in data file ',
     +          curfilnam(:len_trim(curfilnam))
              write(chlog,*) 'ios=',ios
              close(chlog)
              close(chcur)
              return
            endif
            read(chcur,'(a)',iostat=ios) curline
            if (ios.ne.0.or.index(curline,'END').eq.0) then
              ioutstat = -345
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'No END for SPIN TYPES in data 
     +           file ', curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif
 
          else if (index(curline,'INITIAL STATE').ne.0) then

            if (index(curline,'SPIN UP').ne.0) then
              init_state = -1
            else if (index(curline,'SPIN DOWN').ne.0) then
              init_state = -2
            else if (index(curline,'RANDOM').ne.0) then
              init_state = -3
            else if (index(curline,'GROUND STATE').ne.0) then
              init_state = -4
            else if (index(curline,'DIAGONAL').ne.0) then
              read(chcur,*,iostat=ios) init_state
              if (ios.ne.0) then
                ioutstat = -82
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'Error reading DIAGONAL init_state in 
     +            data file ',curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              if (init_state.lt.0.or.init_state.ge.(2**Ltot)) then
                ioutstat = -83
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'Bad value of DIAGONAL init_state in data
     +            file ',curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
              read(chcur,'(a)',iostat=ios) curline
              if (ios.ne.0.or.index(curline,'END').eq.0) then
                ioutstat = -84
                open(chlog,file=logfnam,access='append')
                write(chlog,*) 'No END for DIAGONAL init_state in data 
     +             file ', curfilnam(:len_trim(curfilnam))
                close(chlog)
                close(chcur)
                return
              endif
            else
              ioutstat = -81
              open(chlog,file=logfnam,access='append')
              write(chlog,*) 'Bad value of init_state in data file ',
     +          curfilnam(:len_trim(curfilnam))
              close(chlog)
              close(chcur)
              return
            endif ! choice of init_state
          else if (index(curline,'END').eq.0) then
            open(chlog,file=logfnam,access='append')
            write(chlog,*) '!!!=========>'
            write(chlog,*) 'The line ',curline(:len_trim(curline)),
     +        ' ignored in data file ', curfilnam(:len_trim(curfilnam))
            close(chlog)
          endif
          read(chcur,'(a)',iostat=ios) curline
 
        enddo ! reading all lines from the data file

c reading is finished, EOF or error

        if (ios.eq.0) then 
c if END of data file: check if all right
          if (errsteps.ne.0) then
            ioutstat = -106
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Bad or no value of Nsteps in data file',
     +        curfilnam(:len_trim(curfilnam))
            close(chlog)
          else if (errtau.ne.0.or.dabs(tau).lt.1.0D-15) then
            ioutstat = -107
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Bad or no value of Tau in data file',
     +        curfilnam(:len_trim(curfilnam))
            close(chlog)
          endif     
          if (Lsys.lt.0.or.Lsys.gt.maxLsys) then
            ioutstat = -108
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Lsys is negative or exceeds maxLsys',
     +        curfilnam(:len_trim(curfilnam))
            close(chlog)
          endif
          if (Ltot.lt.0.or.Ltot.gt.maxLtot) then
            ioutstat = -109
            open(chlog,file=logfnam,access='append')
            write(chlog,*) 'Ltot is negative or exceeds maxLtot',
     +        curfilnam(:len_trim(curfilnam))
            close(chlog)
          endif
        else  ! if (ios.ne.0), i.e. not the end of data file
          ioutstat = -103 
          open(chlog,file=logfnam,access='append')
          write(chlog,*) 'Error reading data file',
     +      curfilnam(:len_trim(curfilnam)),ios
          close(chlog)
        endif ! end of if (ios.eq.-1.or.ios.eq.0)
        close(chcur)
        return
        end

c        subroutine erhand(ioutst,ermess,chcur)
c        integer ioutst
c        character*256 ermess
c        
c        open(chlog,file=logfnam,access='append')
c        write(chlog,*) trim(ermess)
c        close(chlog)

