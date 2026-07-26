        subroutine chimstepsPDDG(asp0,bsp0,Jx,Jy,Jz,nconn0,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot0,iRand,nsteps,tau,init_state,timeGl,ioutstat,
     +    psiR,psiI,autonorm,firstime,EpsF,phiR,phiI,dpsiR,dpsiI,
     +    nooutput,
     +    vecResOut,iptVecResOut,mrkResOut)

        implicit none
        include "chsdpar.h"

c in/out variables

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        integer asp0(maxNconn), bsp0(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot)
        double precision tau, timeGl, EpsF
        integer autonorm, firstime
        integer Lsys,Lbath,Ltot0,nconn0,nsteps,iRand,init_state,ioutstat
        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)

        double precision phiI(0:maxStat-1),phiR(0:maxStat-1)
        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)

        integer nooutput

        double precision vecResOut(maxLenVecResOut)
        integer iptVecResOut
        integer mrkResOut(maxLenVecResOut)

c internal variables

        integer nconn, Ltot, asp(maxNconn),bsp(maxNconn) 
        double precision jefx(maxNconn),jefy(maxNconn),jefz(maxNconn)
        double precision hefx(maxLtot),hefy(maxLtot),hefz(maxLtot)
        integer nstates, nstat1,nstat1max,i,k,iconn,j,jp
        integer istep,ia,ib, iord 
c        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)
c        double precision phiI(0:maxStat-1),phiR(0:maxStat-1)
c        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)
c        double precision pbufI(0:maxStat-1),pbufR(0:maxStat-1)
c        double precision dpbI(0:maxStat-1),dpbR(0:maxStat-1)
        double precision xnorm,xnR,xnI, x,y,z
        double precision HmR(2**maxLsys,2**maxLsys)
        double precision HmI(2**maxLsys,2**maxLsys)

        double complex HmT(2**maxLsys,2**maxLsys)
        double complex work((2**maxLsys)*((2**maxLsys)+1))
        double precision w(2**maxLsys),rwork(3*(2**maxLsys)-2)
        integer lwork

        double precision vecR(2**maxLsys,2**maxLsys)
        double precision vecI(2**maxLsys,2**maxLsys)
        double precision tPre, emax,alp,atab(0:nordMax)
        integer nR,nI, nord,ncalc
        logical coeflarge

        double precision xrnd

        double precision dnrm2,ddot
        external dnrm2,ddot
        external daxpy,dscal,dcopy,dswap,ribesl
        external steprnPDDG,addJSxyz,addHSx,addHSy,addHSz
        external actSx,actSy,actSz


c executable statements

        nstates = 2**Ltot0
        Ltot = Ltot0
        nconn = nconn0
c        do i=0,maxStat-1 
c          psiI(i)=psi0I(i)
c          psiR(i)=psi0R(i)
c        enddo
        do i=1,maxNconn
          asp(i) = asp0(i)
          bsp(i) = bsp0(i)
        enddo
        do i=1,2**maxLsys
          do k=1,2**maxLsys
            HmR(k,i) = 0
            HmI(k,i) = 0
          enddo
        enddo

        if (init_state.ne.-10) then ! if initial state is going to be (re)set
          do i=0,nstates-1
            psiR(i) = 0.0D0
            psiI(i) = 0.0D0
          enddo
        endif

        if (init_state.eq.-1) then ! SPIN UP: default value
          psiR(nstates-1) = 1
        else if (init_state.eq.-2) then ! SPIN DOWN
          psiR(0) = 1
        else if (init_state.eq.-3) then ! RANDOM
          do i=0,nstates-1
            xrnd = rand()
            psiR(i) = 2*xrnd-1
          enddo
          do i=0,nstates-1
            xrnd = rand()
            psiI(i) = 2*xrnd-1
          enddo
          xnR = dnrm2(nstates,psiR(0),1)
          xnI = dnrm2(nstates,psiI(0),1)
          xnorm = 1.0D0/sqrt(xnR*xnR + xnI*xnI)
          call dscal(nstates,xnorm,psiR(0),1)
          call dscal(nstates,xnorm,psiI(0),1)
        else if (init_state.eq.-4) then ! GROUND STATE
          nstat1 = 2**Lsys
          nstat1max = 2**maxLsys
          do iconn=1,nconn
            if (asp(iconn).le.Lsys.and.bsp(iconn).le.Lsys) then
              x=-0.25D0*Jx(iconn)
              y=-0.25D0*Jy(iconn)
              z=-0.25D0*Jz(iconn)
              do k=1,nstat1
                psiR(k-1)=1
                call actSx(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSx(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,x,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,x,phiI(0),1,HmI(1,k),1)
                call actSy(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSy(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,y,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,y,phiI(0),1,HmI(1,k),1)
                call actSz(psiR,psiI,asp(iconn),nstat1,dpsiR,dpsiI)
                call actSz(dpsiR,dpsiI,bsp(iconn),nstat1,phiR,phiI)
                call daxpy(nstat1,z,phiR(0),1,HmR(1,k),1)
                call daxpy(nstat1,z,phiI(0),1,HmI(1,k),1)
                psiR(k-1)=0
              enddo
            endif
          enddo
          do i=1,Lsys 
            x=-0.5D0*Hx(i)
            y=-0.5D0*Hy(i)
            z=-0.5D0*Hz(i)
            do k=1,nstat1
              psiR(k-1)=1
              call actSx(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,x,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,x,phiI(0),1,HmI(1,k),1)
              call actSy(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,y,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,y,phiI(0),1,HmI(1,k),1)
              call actSz(psiR,psiI,i,nstat1,phiR,phiI)
              call daxpy(nstat1,z,phiR(0),1,HmR(1,k),1)
              call daxpy(nstat1,z,phiI(0),1,HmI(1,k),1)
              psiR(k-1)=0
            enddo
          enddo

          lwork=((2**maxLsys)*((2**maxLsys)+1))
          do k=1,nstat1
          do i=1,nstat1
            HmT(i,k)=dcmplx(HmR(i,k),HmI(i,k))
          enddo
          enddo
          call zheev('V','U',nstat1,HmT,nstat1max,w,work,lwork,rwork,
     +      ioutstat)
          if (ioutstat.ne.0) then
c            open(chlog,file=logfnam,access='append')
c            write(chlog,*) 'ZHEEV failed, error #', ioutstat
c            close(chlog)
            ioutstat= -201
            return
          endif
          do k=1,nstat1
            vecR(k,1) = dble(HmT(k,1))
            vecI(k,1) = dimag(HmT(k,1))
          enddo

c          jp=max(nstates-1-2**Lsys,0)
          jp=2**Lbath-1
          do k=0,jp
            xrnd = rand()
            dpsiR(k)=2*xrnd-1
          enddo
          do k=0,jp
            xrnd = rand()
            dpsiI(k)=2*xrnd-1
          enddo
          jp = 2**Lbath
          xnR = dnrm2(jp,dpsiR(0),1)
          xnI = dnrm2(jp,dpsiI(0),1)
          xnorm = 1.0D0/sqrt(xnR*xnR + xnI*xnI)
          call dscal(jp,xnorm,dpsiR(0),1)
          call dscal(jp,xnorm,dpsiI(0),1)
          j=1
          jp=0
          do k=0,nstates-1
            psiR(k)= vecR(j,1)*dpsiR(jp) - vecI(j,1)*dpsiI(jp)
            psiI(k)= vecR(j,1)*dpsiI(jp) + vecI(j,1)*dpsiR(jp)
            j=j+1
            if (j.gt.nstat1) then  ! nstat1= 2**Lsys
              j=1
              jp=jp+1
            endif
          enddo
        else if (init_state.ge.0) then
          psiR(init_state)=1
        endif
        if (firstime.ne.0.and.nooutput.eq.0) then
          call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +      phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +      Hx,Hy,Hz,nconn,
     +      vecResOut,iptVecResOut,mrkResOut,ioutstat)
          if (ioutstat.ne.0) then
          return
          endif
        endif

c
c psi(t=0) initialized, let's go on
c

        emax= 0
        do i=1,Ltot
          emax= emax + 0.5D0*(dabs(Hx(i))+dabs(Hy(i))+dabs(Hz(i)))
        enddo

CDIR$ NOVECTOR

        do iconn=1,nconn
          emax= emax + 0.25D0*(dabs(Jx(iconn))+dabs(Jy(iconn))+
     +      dabs(Jz(iconn)))
        enddo
        alp= emax*dabs(tau)
        nord= max(3,int(1.1*alp)) ! preliminary estimate
        if (alp.lt.2*EpsF) then
          atab(0)= 1.0D0
          emax= 1.0D0 ! precaution for too small emax
          if (extraPrec.ne.0) then
            nord= 1
            atab(1)= 0.5D0*tau 
          else
            nord= 0
          endif
        else
          coeflarge= .true.
          do while(coeflarge.and.(nord.le.nordMax))
            call ribesl(alp,0.0D0,nord+1,1,atab(0),ncalc) 
            if (ncalc.le.0) then
c              open(chlog,file=logfnam,access='append')
c              write(chlog,*) 'error in chsteps.f: NCALC=',ncalc
c              close(chlog)
              ioutstat= -202
              return
            endif
            i= 0
            do while(coeflarge.and.(i+2).lt.nord.and.(i+2).lt.ncalc-1)
              i=i+1
              x= abs(atab(i))
              y= abs(atab(i+1))
              z= abs(atab(i+2))
              coeflarge= (x.gt.EpsF).or.(y.gt.EpsF).or.(z.gt.EpsF)
     +          .or.(y.ge.x).or.(z.ge.y)
            enddo 
            if (coeflarge) then
              if (ncalc.lt.(nord+1)) then
c                open(chlog,file=logfnam,access='append')
c                write(chlog,*) 'RIBESL is not precise enough:',
c     +            ' increase EpsF in CHSDPAR.H'
c                close(chlog)
                ioutstat= -203
                return
              else ! nord is not sufficient
                nord= nord+1
              endif
            endif
          enddo
          if (coeflarge) then
c            open(chlog,file=logfnam,access='append')
c            write(chlog,*) 'increase the value of EpsF',
c     +        'or the value of nordMax in CHSDPAR.H'
c            close(chlog)

            vecResOut(1) = EpsF
            vecResOut(2) = alp
            vecResOut(3) = emax
            vecResOut(4) = dble(nord)
            vecResOut(5) = tau

            ioutstat= -204
            return
          endif
          if (extraPrec.ne.0) then
            nord= max(2,i-1)
          else
            nord= i-1
          endif
        endif ! end if (alp.lt.2*EpsF)
c        open(chlog,file=logfnam,access='append')
c        write(chlog,*) 'EPSILON=',EpsF
c        write(chlog,*) 'ALPHA=',alp
c        write(chlog,*) 'EMAX=',emax
c        write(chlog,*) 'final NORD=',nord
c        close(chlog)

        if (tau.lt.0) then
          do iord=1,nord,2
            atab(iord)=-atab(iord)
          enddo
        endif

        do i=1,Ltot
          hefx(i)= -Hx(i)/emax ! \
          hefy(i)= -Hy(i)/emax ! - = -0.5*H*2 (since 2*Ham calculated)
          hefz(i)= -Hz(i)/emax ! /
        enddo
        do iconn=1,nconn
          jefx(iconn)= -0.5D0*Jx(iconn)/emax ! \
          jefy(iconn)= -0.5D0*Jy(iconn)/emax ! - = -0.25*J*2 (since 2*Ham calculated)
          jefz(iconn)= -0.5D0*Jz(iconn)/emax ! /
        enddo
        do istep=1,nsteps
          tPre= timeGl
          timeGl= timeGl+tau
          call dcopy(nstates,psiR(0),1,dpsiR(0),1)
          call dcopy(nstates,psiI(0),1,dpsiI(0),1)
          x= atab(0)
          call dscal(nstates,x,psiR(0),1)
          call dscal(nstates,x,psiI(0),1)
          nR= -1
          nI= -1
          do i=0,nstates-1
            phiR(i)= 0
            phiI(i)= 0
          enddo
          do iord=1, nord
            call dscal(nstates,-1.0D0,phiR(0),1)
            call dscal(nstates,-1.0D0,phiI(0),1)
            do iconn=1,nconn
              x=jefx(iconn)
              y=jefy(iconn)
              z=jefz(iconn)
              ia= asp(iconn)
              ib= bsp(iconn)
c              call addJSxxSyy(x,y,dpsiR,dpsiI,ib,ia,nstates,phiR,phiI)
c              call addJSzz(z,dpsiR,dpsiI,ib,ia,nstates,phiR,phiI)
              call addJSxyz(x,y,z,dpsiR,dpsiI,ib,ia,nstates,phiR,phiI)
            enddo
            do i=1,Ltot 
              x=hefx(i)
              y=hefy(i)
              z=hefz(i)
              if (dabs(x).gt.1.0D-14) then
              call addHSx(x,dpsiR,dpsiI,i,nstates,phiR,phiI)
              endif
              if (dabs(y).gt.1.0D-14) then
              call addHSy(y,dpsiR,dpsiI,i,nstates,phiR,phiI)
              endif
              if (dabs(z).gt.1.0D-14) then
              call addHSz(z,dpsiR,dpsiI,i,nstates,phiR,phiI)
              endif
            enddo
            if (iord.eq.1) then
              call dscal(nstates,0.5D0,phiR(0),1)
              call dscal(nstates,0.5D0,phiI(0),1)
            endif
            call dswap(nstates,phiR(0),1,dpsiR(0),1)
            call dswap(nstates,phiI(0),1,dpsiI(0),1)
            x= atab(iord)
            call daxpy(nstates,2*nR*x,dpsiR(0),1,psiR(0),1)
            call daxpy(nstates,2*nI*x,dpsiI(0),1,psiI(0),1)
            nR= -nR
            nI= -nI
          enddo ! expansion series, do iord=1, nord
c          xnR = dnrm2(nstates,psiR(0),1)
c          xnI = dnrm2(nstates,psiI(0),1)
c          xnorm = 1.0D0/sqrt(xnR*xnR + xnI*xnI)
c          call dscal(nstates,xnorm,psiR(0),1)
c          call dscal(nstates,xnorm,psiI(0),1)

          if (nooutput.eq.0) then
            call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +        phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +        Hx,Hy,Hz,nconn,
     +        vecResOut,iptVecResOut,mrkResOut,ioutstat)
            if (ioutstat.ne.0) then
            return
            endif
          endif

          if (autonorm.eq.1) then
            xnR=dnrm2(nstates,psiR(0),1)
            xnI=dnrm2(nstates,psiI(0),1)
            xnorm= 1.0D0/sqrt(xnR*xnR + xnI*xnI)
            call dscal(nstates,xnorm,psiR(0),1)
            call dscal(nstates,xnorm,psiI(0),1)
          endif
        enddo ! loop over time steps

c        do i=0,maxStat-1 
c          psi0I(i)=psiI(i)
c          psi0R(i)=psiR(i)
c        enddo

        return
        end

