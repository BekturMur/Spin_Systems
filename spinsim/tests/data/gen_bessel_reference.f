c Reference-data generator for the Bessel tests.
c
c Calls the legacy Cody routines rjbesl.f / ribesl.f exactly as
c chstepsPDDGnmr.f:98 and chimstepsPDDGnmr.f:91 do, and dumps the resulting
c coefficient tables. The C++ implementation is checked against this output, so
c the two are compared through the same interface the physics code uses.
c
c Build:  gfortran -O2 gen_bessel_reference.f ../../legacy/rjbesl.f \
c                      ../../legacy/ribesl.f -o gen_bessel_reference
c Usage:  ./gen_bessel_reference > bessel_reference.txt
c
c Output format, one record per (kind, alpha):
c   kind alpha nb ncalc
c   b(1) b(2) ... b(nb)     (one value per line, %.17e)

        program genbessel
        implicit none

        integer nbMax
        parameter (nbMax=400)
        double precision b(nbMax)
        integer ncalc, nb, i, ia
        double precision alp

c Arguments spanning the regimes the propagator actually meets: alp = emax*tau
c ranges from far below 1 (tiny time step) to a few hundred (long step or
c strong coupling). Includes values where the sequence underflows to zero.
        integer nalp
        parameter (nalp=14)
        double precision alptab(nalp)
        data alptab /1.0D-6, 1.0D-3, 0.01D0, 0.1D0, 0.5D0, 1.0D0,
     +    2.0D0, 5.0D0, 10.0D0, 25.0D0, 50.0D0, 100.0D0, 200.0D0,
     +    350.0D0/

        nb = nbMax

        do ia=1,nalp
          alp = alptab(ia)
          call rjbesl(alp,0.0D0,nb,b,ncalc)
          write(*,'(a,1x,e27.17e3,1x,i6,1x,i6)') 'J', alp, nb, ncalc
          do i=1,nb
            write(*,'(e27.17e3)') b(i)
          enddo
        enddo

        do ia=1,nalp
          alp = alptab(ia)
          call ribesl(alp,0.0D0,nb,1,b,ncalc)
          write(*,'(a,1x,e27.17e3,1x,i6,1x,i6)') 'I', alp, nb, ncalc
          do i=1,nb
            write(*,'(e27.17e3)') b(i)
          enddo
c The exponent width is forced to three digits throughout. Under a plain
c E26.17 the compiler drops the 'E' once the exponent needs three digits,
c emitting "0.233...-107"; a C++ float parser then silently reads 0.233.
        enddo

        end
