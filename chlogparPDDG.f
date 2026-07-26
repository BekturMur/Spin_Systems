        subroutine chlogparPDDG(asp,bsp,Jx,Jy,Jz,nconn,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot,iRand,nsteps,tau,init_state,autonorm,imflag,
     +    ispecflag,nooutput,iNtype)

        implicit none
        include "chsdpar.h"

c in/out variables

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        integer asp(maxNconn), bsp(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot), tau
        integer Lsys,Lbath,Ltot,nconn, nsteps, iRand, init_state
        integer autonorm, imflag
        integer ispecflag, nooutput
        integer iNtype(maxLtot)

c internal variables

        integer i

c executable statements

        open(chlog, file=logfnam, access='append')
        write(chlog,*) 'Total number of spins: ', Ltot
        write(chlog,*) 'Number of system spins: ', Lsys
        write(chlog,*) 'Number of bath spins: ', Lbath
        write(chlog,'(a$)') ' Initial state: '
        if (init_state.eq.-1) then
          write(chlog,*) 'spin up'
        else if (init_state.eq.-2) then
          write(chlog,*) 'spin down'
        else if (init_state.eq.-3) then
          write(chlog,*) 'random ; ranfx starts from', iRand
        else if (init_state.eq.-4) then
          write(chlog,*) 'ground state'
        else if (init_state.eq.-10) then
          write(chlog,*) 'previous'
        else if (init_state.ge.0) then
          write(chlog,*) 'diagonal #',init_state
        else
          write(chlog,*) 'SOME PROBLEM!!!'
          close(chlog)
          return
        endif
        write(chlog,*) 'Integration time step: ', tau
        write(chlog,*) 'Number of such time steps: ', nsteps
        write(chlog,*) 'Automatic normalization: ',autonorm
        if (imflag.eq.1) then
          write(chlog,*) 'Imaginary time propagation'
        else if (imflag.eq.0) then
          write(chlog,*) 'Real time propagation'
        else 
          write(chlog,*) 'Special propagation',ispecflag
        endif
        write(chlog,*) 'No output: ',nooutput

        write(chlog,*) ' '
        write(chlog,*) 'Spin types'
        do i=1,Ltot
          write(chlog,'(i5,a$)') iNtype(i),' , '
        enddo
        write(chlog,*) ' '

        write(chlog,*) ' '
        do i=1,Ltot
          write(chlog,'(a,i4,2a,3f10.5)') 'Spin #',i,' local fields', 
     +     ' (x,y,z): ', Hx(i), Hy(i), Hz(i)
        enddo

        write(chlog,*) ' '
        write(chlog,*) 'Exchange interactions: connection, spin pair,', 
     +   ' Jxx, Jyy, Jzz'
        do i=1,nconn
          write(chlog,'(i4,a,i4,a,i4,a,3f10.5)') i,' (',asp(i),',',
     +      bsp(i),') :', Jx(i),Jy(i),Jz(i)
        enddo

        close(chlog)

        return
        end


