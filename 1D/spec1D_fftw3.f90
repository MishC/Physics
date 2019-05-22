      PROGRAM PSIrthet
!
! Split-Operator TDSE in 1D. 
! coordinates. r-propagation by
! 
! To compile gfortran -O2 -o spec1D_fftw3 spec1D_fftw3.f90 -lfftw3
! Test example: ./spec1D_fftw3 < in_stability/e0 (zero interaction, 1D couloum like potential)
!               ./spec1D_fftw3 < in_stability/e1 (non-zero interaction)
!
      implicit none

      integer                   :: nsteps           ! Number of timestep      
      integer                   :: nr ! Problem size

      integer i,j,k,kk,itype,ii, isign
      double precision t,dr,dt,Tprop, Tstart, dk     

! Input variables:
      integer outlen,  nut , ninit, powerof2
      integer			:: status
      character*80 outf, psi_infile
      double precision tprintnext, wL, e0, rmax ,  p, pi, x
      
! wavefunction:
      double complex,   dimension(:),allocatable 	:: psi, cxy, cabsorb
      double complex, dimension(:),allocatable	        :: matarg
! new version: matark -> exp(-i matarg*dt);
      double complex, dimension(:),allocatable	        :: arg, pkexp
      double precision, dimension(:),allocatable	:: absorb, inpsi
      double precision, dimension(:),allocatable	:: zk, b
      double precision, dimension(:),allocatable	:: pk
      double precision rl,rlsq,coeff,one,zero
      double precision zmax, zlim, absorber
! Rydberg ting. Puls er e0*sin^2(alfa t)*cos(wLt)
      double precision alfa,ft
     real start,finish
     double complex cmi
     integer*8                                          :: planfwd       ! FFTW plans
     integer*8                                          :: planbwd       ! 
!
!! Read input, broadcast it to all processors, compute local 
!! sizes and allocate arrays
!    
      cmi=dcmplx(0,1); pi=acos(-1.0d0) 

      read (*,*) powerof2, rmax, ninit
      read (*,'(A)') psi_infile
      read (*,*) p	! b is coulomb softening factor and p is position of coulomb pot from matlab
      read (*,*) Tprop, nut, Tstart, dt
      read (*,*) e0
      read (*,*) wL
      read (*,'(A)')outf
! Find real lenght (LL) of outf-string:		
      do i = 1,79
       if (outf(i+1:i+1) .eq. ' ') then
         outlen = i
         goto 68
       endif
      enddo       
      write(*,*)' Something wrong with outf - I stop'
      stop
 68   continue

     
      nr=2**(powerof2); dr=rmax/(nr-1)
      allocate (inpsi(nr), STAT=status); allocate (psi(nr), STAT=status); 
      allocate (matarg(nr), STAT=status); allocate (arg(nr), STAT=status)
      allocate (absorb(nr), STAT=status); allocate (zk(nr), STAT=status); allocate (pk(nr), STAT=status)
      allocate (pkexp(nr), STAT=status); allocate(b(nr), STAT=status);
      allocate (cabsorb(nr), STAT=status);
      do j=1,nr
        zk(j) = (j-1)*dr     
      enddo
      !b=0.1+abs(10*(zk-p)) !b=0.1+7*(abs(zk-p))**(1.75)
      b=2.0
! stable finite oneD potential:
      matarg=(-1.0d0/sqrt((zk-p)**2+b))*(-cmi*dt)
! Coulomb potential:
      b=0; matarg=(-1.0d0/sqrt((zk-p)**2+b))*(-cmi*dt)
      write(81,*)zk-p
      write(83,*)-1.0d0/sqrt((zk-p)**2+b)
      dk=2*pi/(nr*dr)

! cft ala Hicham:
      do j=0,nr/2-1
        pk(j+1) =  j*dk 
       enddo
       do j=nr/2,nr-1
	 pk(j+1) = (j-nr)*dk
       enddo

      call fft_init(nr, psi, planfwd, planbwd)  
      pkexp=cdexp(-cmi*dt/2*pk*pk)

      open(10,file=psi_infile)   
        do i = 1,ninit-1            
          read(10,*)x
        enddo
        read(10,*)inpsi
        
      close(10)
      psi=inpsi
! 1s in hydrogen:
  !    psi=(zk-p)*dexp(-abs(zk-p)); psi=psi/(sum(abs(psi)**2)*dr)
 !     write(82,*)abs(psi)**2

! New absorbator:
      zmax=rmax/2; zlim=0.95*zmax;
      do ii = 1,nr
       if (abs(zk(ii)-p) .ge. zlim) then
          absorb(ii)=absorber(dabs(zk(ii)-p),zlim,zmax)
        else
          absorb(ii)=1.0d0   
        endif
      enddo    
     
!  NO ABSORBATOR: absorb=1
!
! Initiate:
!
     nsteps=Tprop/dt                
     write(*,*)'Starting: nr =    ',nr, ' nsteps = ',nsteps, sum(abs(psi)**2)*dr
!
! Set up initial fl:
! 
     t=Tstart; nsteps=nsteps-1; tprintnext=Tstart+Tprop/nut;

      open(9,file=outf(1:outlen)//'_z.dat'); write(9,*) zk(1:nr); close(9)
      open(10,file=outf(1:outlen)//'_T.dat')
      open(11,file=outf(1:outlen)//'_psiSQ.dat')
      open(12,file=outf(1:outlen)//'_psiR.dat')
      open(13,file=outf(1:outlen)//'_psiI.dat')
      print *, 'Opened'
      call print(t,nr,psi,outf,outlen)  
!
! 1. Iterate first step dt/2 in r-direction:
!                                           
! isign=1 is flagging "Forward" transform, -1 the inverse:
!     isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*exp(cmplx(0.d0,-dt*pk*pk/4.0d0))
!     isign = -1; call fft(nr,isign,cxy); psi=cxy;
!
     call dfftw_execute(planfwd); psi=psi*cdexp(-cmi*dt*pk*pk/4.0d0)
     call dfftw_execute(planbwd); psi=psi/real(nr)
!
!--------------------------------------------------------------
! main loop: Timepropagation from t=0 to t=(n_steps-1)*dt:
!--------------------------------------------------------------
!
      call cpu_time(start)     
      do kk=1,nsteps
!do kk=1,1000
!        if (mod(kk,10000) .eq. 0) write(*,*)kk, sum(abs(psi)**2)*dr
        t=t+dt/2
! Length gauge:
        if (t .lt. 0) then
           ft = -e0*sin(wL*t)
        else
           ft = -1.5*e0*sin(1.5*wL*t)
        endif
        arg=-cmi*dt*(zk-p)*ft
 
! KH-frame - Single cycle THz pulse. 
!        if  (t .lt. 0)  then
!            ft=e0/wL**2*sin(wL*t)
!        else
!            ft=e0/(1.5*wL**2)*sin(1.5*wL*t)
!        endif
!        arg=-cmi*dt*(mat1-1/sqrt(mat2+ft**2-2.0d0*mat3*ft))

        psi=psi*exp(matarg+arg)*absorb
               
!        isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*pkexp
!        isign = -1; call fft(nr,isign,cxy); psi=cxy;
        call dfftw_execute(planfwd); psi=psi*pkexp
        call dfftw_execute(planbwd); psi=psi/real(nr)

        t=t+dt/2    
        if  (t .gt. tprintnext) then
          write(*,'(2F18.6,F9.5,A)')t,sum(abs(psi)**2)*(zk(2)-zk(1)), real(kk)/nsteps*100, ' % finnished'
          open(95,file='monitor_time.txt')
          write(95,'(2F13.3,F8.4,A)')t,sum(abs(psi)**2)*dr, real(kk)/nsteps*100, ' % finnished'
          close(95)
          tprintnext=tprintnext+Tprop/nut
          call print(t,nr,psi,outf,outlen)
        endif

     enddo
     call cpu_time(finish)
     write(*,*)'Time: ',finish-start

!--------------------------------------------------------------
! END Main loop. What is left is the  final steps & printout:
!--------------------------------------------------------------
     
!
! new version:   matarg -> exp(cmplx(0.0d0,-idt*matarg(j,i)))
!   
     t=t+dt/2
     ft = -1.5*e0*sin(1.5*wL*t)
     arg=-cmi*dt*(zk-p)*ft
      
! KH-frame: 
!    ft=e0/(1.5*wL**2)*sin(1.5*wL*t)
!     arg=-cmi*dt*(mat1-1/sqrt(mat2+ft**2-2.0d0*mat3*ft))

     psi=psi*cdexp(matarg+arg)

! isign=1 is flagging "Forward" transform, -1 the inverse:
!     isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*exp(cmplx(0.d0,-dt*pk*pk/4.0d0))
!     isign = -1; call fft(nr,isign,cxy); psi=cxy;

     call dfftw_execute(planfwd); psi=psi*cdexp(-cmi*dt*pk*pk/4.0d0)
     call dfftw_execute(planbwd); psi=psi/real(nr)
     
     call print(t,nr,psi,outf,outlen)
              write(*,*) outf
        close(10)
        close(11)
	close(12); close(13)

     stop
     end


      subroutine print(t, nr, psi, outf,outlen)          
      implicit none

! wavefunction:

      integer                           :: nr
      integer                           :: nrp
      double complex, dimension(nr) :: psi
      double precision t
      character*80 outf     
      integer outlen
!
! Local variables
!
      integer i,j,kk,kkk,jjj,ierror,nsize


        
                   write(10,'(5000000F16.9)')t
          write(11,'(5000000F16.9)')(abs(psi(j))**2,j=1,nr)
          write(12,'(5000000F16.9)')(real(psi(j)),j=1,nr)
          write(13,'(5000000F16.9)')(aimag(psi(j)),j=1,nr)

          !do j=1,nr
            !write(11,'(F16.9)')(abs(psi(j))**2)
            !write(12,'(F16.9)')(real(psi(j)))
           ! write(13,'(F16.9)')(aimag(psi(j)))
          !enddo
	!write(*,*)'leaving print:'

      return
      end

subroutine fft(n,dir,xy)
  ! http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/dft/
  ! http://paulbourke.net/miscellaneous/dft/
  !   computes an in-place complex-to-complex FFT 
  !   x and y are the real and imaginary arrays of 2^m points.
  !   dir =  1 gives forward transform
  !   dir = -1 gives reverse transform
  !    do i=0,n/2-1
  !      k(i) = 2*pi * (i/(n*dt))
  !    end do
  !    do i=n/2,n-1
  !      k(i) = 2*pi * ((i-n)/(n*dt))
  !    end do
  implicit none
  integer, intent(in) :: n,dir
  double complex, intent(inout) :: xy(1:n)
  
  integer :: m,i,i1,j,k,i2,l,l1,l2
  double precision ::  c1,c2,tx,ty,t1,t2,u1,u2,z
  double precision, allocatable ::  x(:),y(:)  

  m = nint(log(real(n))/log(2d0))
  if (2**m.ne.n) then
    write(0,'("fft: n is not a power of 2")')
    stop
  end if

  allocate(x(0:n-1),y(0:n-1))
  x(0:n-1) = real(xy(1:n))
  y(0:n-1) = imag(xy(1:n))

  ! do the bit reversal 
  i2 = n/2
  j = 0
  do i=0,n-2
    if (i < j) then
      tx = x(i)
      ty = y(i)
      x(i) = x(j)
      y(i) = y(j)
      x(j) = tx
      y(j) = ty
    end if
    k = i2
    do while (k <= j)
      j = j-k
      k = k/2
    end do
    j = j+k
  end do

  ! compute the FFT 
  c1 = -1 
  c2 = 0
  l2 = 1
  do l=0,m-1
    l1 = l2
    l2 = l2*2
    u1 = 1 
    u2 = 0
    do j=0,l1-1
      do i=j,n-1,l2
        i1 = i + l1
        t1 = u1 * x(i1) - u2 * y(i1)
        t2 = u1 * y(i1) + u2 * x(i1)
        x(i1) = x(i) - t1 
        y(i1) = y(i) - t2
        x(i) = x(i) + t1
        y(i) = y(i) + t2
      end do
      z =  u1 * c1 - u2 * c2
      u2 = u1 * c2 + u2 * c1
      u1 = z
    end do
    c2 = sqrt((1 - c1) / 2)
    if (dir.eq.1)  c2 = -c2
    c1 = sqrt((1 + c1) / 2)
  end do

  ! scaling for forward transform 
  if (dir.eq.1) then
    xy = dcmplx(x/n,y/n)
  else
    xy = dcmplx(x,y)
  end if

  return

!  write(*,'(i9,2es24.16)') (i-n, x(i),y(i), i=n/2,n-1)
!  write(*,'(i9,2es24.16)') (i,   x(i),y(i), i=0,n/2-1)

end subroutine fft

      subroutine fft_init(nr, psi, planfwd,planbwd)
!
! Purpose:
! Set up FFTW plans with MEASURE parameter to find optimal algorithm 
!
! Date        By                  Comment
! --------    -----------------   ---------------------
! 14.11.06    Raymond Nepstad     First version
!
      implicit none
      include "/usr/include/fftw3.f"

      integer, parameter                                   :: dbl = selected_real_kind(p=13)
      integer, intent(in)                                  :: nr
      integer :: nzp
      double complex,  dimension(nr), intent(inout)        :: psi
      integer*8,intent(out)                                :: planfwd
      integer*8,intent(out)                                :: planbwd
      !integer*8 :: FFTW_FORWARD, FFTW_MEASURE, FFTW_BACKWARD
!
!      call dfftw_plan_many_r2r(planfwd,&         ! 0. Transform pointer (integer*8)
      nzp=1
      call dfftw_plan_many_dft(planfwd,&         ! 0. Transform pointer (integer*8)
                         1,&               ! 1. Rank  
                         nr,&              ! 2. Size of transform 
                         nzp,&             ! 3. # of transforms to compute
                         psi,&              ! 4. Input array ref
                         nr,&              ! 5. Size of input array
                         1,&               ! 6. Transform stride
                         nr,&              ! 7. Distance to next transform
                         psi,&              ! 8. Output array
                         nr,&              ! 9. Size of output array
                         1,&               ! 10. Output array stride
                         nr,&              ! 11. Distance to next output array
                         FFTW_BACKWARD,&    ! 12. Sine-transform flag
                         FFTW_MEASURE)     ! 13. Transform flag
						 
!call dfftw_plan_many_r2r(planbwd,&         ! 0. Transform pointer (integer*8)
        call dfftw_plan_many_dft(planbwd,&         ! 0. Transform pointer (integer*8)
                         1,&               ! 1. Rank  
                         nr,&              ! 2. Size of transform 
                         nzp,&             ! 3. # of transforms to compute
                         psi,&              ! 4. Input array ref
                         nr,&              ! 5. Size of input array
                         1,&               ! 6. Transform stride
                         nr,&              ! 7. Distance to next transform
                         psi,&              ! 8. Output array
                         nr,&              ! 9. Size of output array
                         1,&               ! 10. Output array stride
                         nr,&              ! 11. Distance to next output array
                         FFTW_FORWARD,&   ! 12. Sine-transform flag
                         FFTW_MEASURE)     ! 13. Transform flag

        end subroutine fft_init

double precision function absorber(x,x0,xmax)
!      
! Computes the value of a logistic absorber - in matlab:
!    L=1; x0=0; xmax=10; k=12/(2*xmax); np=64; x=linspace(-xmax,xmax,np);
!    y=L./(1+exp(-k*(x-x0))); yy=1-y; plot(x,yy,'r-')       
!
     implicit none      
     double precision                     :: x, x0, xmax     
     double precision                     :: k, rmidpoint,y
 
     rmidpoint=(xmax+x0)/2.0d0; k=15.0d0/(xmax-x0)
     y=1/(1+exp(-k*(x-rmidpoint)))

     absorber=1-y
     return
     end
