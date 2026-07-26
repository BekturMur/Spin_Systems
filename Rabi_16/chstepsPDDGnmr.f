        subroutine chstepsPDDGnmr(asp0,bsp0,Jx,Jy,Jz,nconn0,Hx,Hy,Hz,
     +    Lsys,Lbath,Ltot0,nsteps,tau,timeGl,ioutstat,
     +    psiR,psiI,autonorm,EpsF,phiR,phiI,dpsiR,dpsiI,
     +    nooutput,
     +    vecResOut,iptVecResOut,mrkResOut)

        implicit none
        include "chsdpar.h"

c in/out variables

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        integer asp0(maxNconn), bsp0(maxNconn)
        double precision Hx(maxLtot), Hy(maxLtot), Hz(maxLtot)
        double precision tau, timeGl, EpsF
        integer autonorm
        integer Lsys,Lbath,Ltot0,nconn0,nsteps,ioutstat
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
        integer nstates, i,k,iconn,j,jp
        integer istep,ia,ib, iord 
c        double precision psiI(0:maxStat-1),psiR(0:maxStat-1)
c        double precision phiI(0:maxStat-1),phiR(0:maxStat-1)
c        double precision dpsiI(0:maxStat-1),dpsiR(0:maxStat-1)
c        double precision pbufI(0:maxStat-1),pbufR(0:maxStat-1)
c        double precision dpbI(0:maxStat-1),dpbR(0:maxStat-1)
        double precision xnorm,xnR,xnI, x,y,z

        double precision vecR(2**maxLsys,2**maxLsys)
        double precision vecI(2**maxLsys,2**maxLsys)
        double precision tPre, emax,alp,atab(0:nordMax)
        integer nR,nI, nord,ncalc
        logical coeflarge

        double precision xrnd

        double precision dnrm2,ddot
        external dnrm2,ddot
        external daxpy,dscal,dcopy,dswap,rjbesl

c executable statements

        nstates = 2**Ltot0
        Ltot = Ltot0
        nconn = nconn0
        do i=1,maxNconn
          asp(i) = asp0(i)
          bsp(i) = bsp0(i)
        enddo

        emax= 0.0D0
        do i=1,Ltot
          emax= emax + 0.5D0*(dabs(Hx(i))+dabs(Hy(i))+dabs(Hz(i)))
          
c              open(chlog,file=logfnam,access='append')
c              write(chlog,*) i,emax
c              close(chlog)

        enddo
!DIR$ NOVECTOR
        do iconn=1,nconn
          emax= emax + 0.25D0*(dabs(Jx(iconn))+dabs(Jy(iconn))+
     +      dabs(Jz(iconn)))
        enddo
          
c              open(chlog,file=logfnam,access='append')
c              write(chlog,*) nconn,emax
c              close(chlog)

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
            call rjbesl(alp,0.0D0,nord+1,atab(0),ncalc) 
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
c                write(chlog,*) 'RJBESL is not precise enough:',
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



c        open(77,file='chkord.out',access='append')
c        write(77,*) 'EPSILON=',EpsF
c        write(77,*) 'ALPHA=',alp
c        write(77,*) 'EMAX=',emax
c        write(77,*) 'final NORD=',nord
c        close(77)




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


c              open(chlog,file=logfnam,access='append')
c              do jp=1,nconn
c              write(chlog,*) jp, jefx(jp),jefy(jp),jefz(jp)
c              enddo
c              close(chlog)


        do istep=1,nsteps
          tPre= timeGl
          timeGl= timeGl+tau
          call dcopy(nstates,psiR(0),1,dpsiR(0),1)
          call dcopy(nstates,psiI(0),1,dpsiI(0),1)
          x= atab(0)
          call dscal(nstates,x,psiR(0),1)
          call dscal(nstates,x,psiI(0),1)
          nR= 0
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
            call daxpy(nstates,-2*nI*x,dpsiI(0),1,psiR(0),1)
            call daxpy(nstates,2*nR*x,dpsiI(0),1,psiI(0),1)
            call daxpy(nstates,2*nI*x,dpsiR(0),1,psiI(0),1)
            jp= nR
            nR= nI
            nI= -jp
          enddo ! expansion series, do iord=1, nord

          if (autonorm.eq.1) then
            xnR=dnrm2(nstates,psiR(0),1)
            xnI=dnrm2(nstates,psiI(0),1)
            xnorm= 1.0D0/sqrt(xnR*xnR + xnI*xnI)
            call dscal(nstates,xnorm,psiR(0),1)
            call dscal(nstates,xnorm,psiI(0),1)
          endif

          if (nooutput.eq.0) then
            call steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +        phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +        Hx,Hy,Hz,nconn,
     +        vecResOut,iptVecResOut,mrkResOut,ioutstat)
            if (ioutstat.ne.0) then
            return
            endif
          endif

        enddo ! loop over time steps

c        do i=0,maxStat-1 
c          psi0I(i)=psiI(i)
c          psi0R(i)=psiR(i)
c        enddo

        return
        end

c
c
c
        subroutine steprnPDDG(timeGl,psiR,psiI,Lsys,Ltot,nstates,
     +    phiR,phiI,dpsiR,dpsiI,vecR,vecI,asp,bsp,Jx,Jy,Jz,
     +    Hx,Hy,Hz,nconn,
     +    vecResOut,iptVecResOut,mrkResOut,ioutstat)
        implicit none
        include "chsdpar.h"

        double precision timeGl
        integer Lsys,Ltot,nstates
        double precision psiR(0:maxStat-1), psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1), phiI(0:maxStat-1)
        double precision dpsiR(0:maxStat-1), dpsiI(0:maxStat-1)
        double precision vecR(2**maxLsys,2**maxLsys)
        double precision vecI(2**maxLsys,2**maxLsys)

        double precision Jx(maxNconn),Jy(maxNconn),Jz(maxNconn)
        double precision Hx(maxLtot),Hy(maxLtot),Hz(maxLtot)
        integer asp(maxNconn), bsp(maxNconn), nconn

        double precision vecResOut(maxLenVecResOut)
        integer iptVecResOut,ioutstat
        integer mrkResOut(maxLenVecResOut)

        double precision dnrm2,ddot
        external dnrm2,ddot

        integer i,k,l,ik,il, ip1, lenrec
        double precision x,y,z,sxy12,sxy21,sxz12,sxz21,syz12,syz21
        double precision xnR,xnI, rR1,rI1,rR2,rI2

        ip1 = 2**Lsys
        do l=1,ip1
          do k=1,ip1
            vecR(k,l) = 0
            vecI(k,l) = 0
          enddo
        enddo
        do i=0,nstates-1,ip1
          do l=1,ip1
            il= i+l-1
            x= psiR(il)
            y= psiI(il)
            do k=1,ip1
              ik= i+k-1
cc              vecR(k,l)= vecR(k,l)+psiR(ik)*psiR(il)+psiI(ik)*psiI(il) - density matrix
cc              vecI(k,l)= vecI(k,l)+psiI(ik)*psiR(il)-psiR(ik)*psiI(il) /
              vecR(k,l)= vecR(k,l) + psiR(ik)*x + psiI(ik)*y
              vecI(k,l)= vecI(k,l) + psiI(ik)*x - psiR(ik)*y
            enddo
          enddo
        enddo
        xnR=dnrm2(nstates,psiR(0),1)
        xnI=dnrm2(nstates,psiI(0),1)

        lenrec = ip1*ip1 + 2

        if (iptVecResOut+lenrec.gt.maxLenVecResOut) then
          ioutstat = -456
          return
        endif

        mrkResOut(iptVecResOut+1) = lenrec
        VecResOut(iptVecResOut+1) = timeGl
        VecResOut(iptVecResOut+2) = dsqrt(xnR*xnR + xnI*xnI)
        iptVecResOut = iptVecResOut + 2
        do i=1,ip1
          vecResOut(iptVecResOut+1) = vecR(i,i)
          iptVecResOut = iptVecResOut + 1
        enddo
        do i=2,ip1
          do k=1,i-1
            vecResOut(iptVecResOut+1) = vecR(i,k)
            vecResOut(iptVecResOut+2) = vecI(i,k)
            iptVecResOut = iptVecResOut + 2
          enddo
        enddo

        return
        end

c
c |0> <-> |down>, |1> <-> |up>
c
c Try as an alternative: (i) BLAS and/or (ii) phi(i1:i2)=psi(j1:j2)
c

        subroutine addJSxyz(xJ,yJ,zJ,psiR,psiI,i0,j0,nstates,phiR,phiI)
c i0 > j0 !
        include "chsdpar.h"
        integer i0,j0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xJ,yJ,zJ

        integer ip0,ip1,j,k,l,jp0,jp1,j1,j2,k1,l1
        double precision xpy,xmy
        xpy = xJ + yJ
        xmy = xJ - yJ
        ip0 = 2**(i0-1)
        ip1 = 2*ip0  ! 2**i0
        jp0 = 2**(j0-1)
        jp1 = 2*jp0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1,jp1
            do k1=k,k+jp0-1 
              j = k1+ip0
              j1 = k1+jp0
              j2 = k1+ip0+jp0
              phiR(j2) = phiR(j2) + xmy*psiR(k1) + zJ*psiR(j2)
              phiR(k1) = phiR(k1) + xmy*psiR(j2) + zJ*psiR(k1)
              phiI(j2) = phiI(j2) + xmy*psiI(k1) + zJ*psiI(j2)
              phiI(k1) = phiI(k1) + xmy*psiI(j2) + zJ*psiI(k1)

              phiR(j) = phiR(j) + xpy*psiR(j1) - zJ*psiR(j)
              phiR(j1) = phiR(j1) + xpy*psiR(j) - zJ*psiR(j1)
              phiI(j) = phiI(j) + xpy*psiI(j1) - zJ*psiI(j)
              phiI(j1) = phiI(j1) + xpy*psiI(j) - zJ*psiI(j1)
            enddo
          enddo
        enddo
        return
        end


        subroutine addJSxxSyy(xJ,yJ,psiR,psiI,i0,j0,nstates,phiR,phiI)
c i0 > j0 !
        include "chsdpar.h"
        integer i0,j0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xJ,yJ

        integer ip0,ip1,j,k,l,jp0,jp1,j1,j2,k1,l1
        double precision xpy,xmy
        xpy = xJ + yJ
        xmy = xJ - yJ
        ip0 = 2**(i0-1)
        ip1 = 2*ip0  ! 2**i0
        jp0 = 2**(j0-1)
        jp1 = 2*jp0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1,jp1
            do k1=k,k+jp0-1 
              j = k1+ip0
              j1 = k1+jp0
              j2 = k1+ip0+jp0
              phiR(j2) = phiR(j2) + xmy*psiR(k1)
              phiR(k1) = phiR(k1) + xmy*psiR(j2)
              phiI(j2) = phiI(j2) + xmy*psiI(k1)
              phiI(k1) = phiI(k1) + xmy*psiI(j2)

              phiR(j) = phiR(j) + xpy*psiR(j1)
              phiR(j1) = phiR(j1) + xpy*psiR(j)
              phiI(j) = phiI(j) + xpy*psiI(j1)
              phiI(j1) = phiI(j1) + xpy*psiI(j)
            enddo
          enddo
        enddo
        return
        end

        subroutine addJSzz(xJ,psiR,psiI,i0,j0,nstates,phiR,phiI)
c i0 > j0 !
        include "chsdpar.h"
        integer i0,j0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xJ

        integer ip0,ip1,j,k,l,jp0,jp1,j1,j2,k1,l1
        ip0 = 2**(i0-1)
        ip1 = 2*ip0  ! 2**i0
        jp0 = 2**(j0-1)
        jp1 = 2*jp0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1,jp1
            do k1=k,k+jp0-1 
              j = k1+ip0
              j1 = k1+jp0
              j2 = k1+ip0+jp0
              phiR(j2) = phiR(j2) + xJ*psiR(j2)
              phiR(k1) = phiR(k1) + xJ*psiR(k1)
              phiI(j2) = phiI(j2) + xJ*psiI(j2)
              phiI(k1) = phiI(k1) + xJ*psiI(k1)

              phiR(j) = phiR(j) - xJ*psiR(j)
              phiR(j1) = phiR(j1) - xJ*psiR(j1)
              phiI(j) = phiI(j) - xJ*psiI(j)
              phiI(j1) = phiI(j1) - xJ*psiI(j1)
            enddo
          enddo
        enddo
        return
        end

        subroutine addHSx(xH,psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xH

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0  ! 2**i0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = phiR(j) + xH*psiR(k)
            phiR(k) = phiR(k) + xH*psiR(j)
            phiI(j) = phiI(j) + xH*psiI(k)
            phiI(k) = phiI(k) + xH*psiI(j)
          enddo
        enddo
        return
        end

        subroutine addHSy(xH,psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xH

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = phiR(j) + xH*psiI(k)
            phiR(k) = phiR(k) - xH*psiI(j)
            phiI(j) = phiI(j) - xH*psiR(k)
            phiI(k) = phiI(k) + xH*psiR(j)
          enddo
        enddo
        return
        end

        subroutine addHSz(xH,psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)
        double precision xH 

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = phiR(j) + xH*psiR(j)
            phiR(k) = phiR(k) - xH*psiR(k)
            phiI(j) = phiI(j) + xH*psiI(j)
            phiI(k) = phiI(k) - xH*psiI(k)
          enddo
        enddo
        return
        end

        subroutine actSx(psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0  ! 2**i0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = psiR(k)
            phiR(k) = psiR(j)
            phiI(j) = psiI(k)
            phiI(k) = psiI(j)
          enddo
        enddo
        return
        end

        subroutine actSy(psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = psiI(k)
            phiR(k) = -psiI(j)
            phiI(j) = -psiR(k)
            phiI(k) = psiR(j)
          enddo
        enddo
        return
        end

        subroutine actSz(psiR,psiI,i0,nstates,phiR,phiI)
        include "chsdpar.h"
        integer i0, nstates
        double precision psiR(0:maxStat-1)
        double precision psiI(0:maxStat-1)
        double precision phiR(0:maxStat-1)
        double precision phiI(0:maxStat-1)

        integer ip0,ip1,j,k,l
        ip0 = 2**(i0-1)
        ip1 = 2*ip0
        do l=0,nstates-1,ip1
          do k=l,l+ip0-1
            j = k+ip0
            phiR(j) = psiR(j)
            phiR(k) = -psiR(k)
            phiI(j) = psiI(j)
            phiI(k) = -psiI(k)
          enddo
        enddo
        return
        end

