      PROGRAM PSIrthet
!
! Split-Operator TDSE in 1D. 
! coordinates. r-propagation by
! 
! Psi(z_i, r_j) (Wavefunction multiplied with r) i=1..nr; j=1..=L_max+1 
! fl(l,j) = f_l(r_j)
! z_i are zeros of P_{L_max} + 1
! w_i are corresponding weights
! Array of r values: rvs(r_i), i=1..nr
!   For fft r(1) = -rmax, r(nr)=rmax-dr, dr=2*rmax/nr
!   For krankn r(1)=0, r(nr) = rmax-dr, dr=rmax/nr
! Array of z values: zvs(z_i), i=1,
!
! Zeros and weights for the Gauss-Legendre rules  are found from 
! by a public domain routine developed by W. Gautschi and available
! from Netlib.  
!
! Any single m value can be treated.
!
      implicit none

      integer                   :: nsteps           ! Number of timestep      
      integer                   :: nr ! Problem size

      integer i,j,k,kk,itype,ii, isign
      double precision t,dr,dt,Tprop, Tstart, dk     

! Input variables:
      integer outlen,  nut , ninit
      integer			:: status
      character*80 outf, psi_infile
      double precision tprintnext, wL, e0, rmax ,  p, pi, powerof2, x
      
! wavefunction:
      double complex,   dimension(:),allocatable 	:: psi, cxy
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
     complex cmi
     integer*8                                          :: planfwd       ! FFTW plans
     integer*8                                          :: planbwd       ! 
!
!! Read input, broadcast it to all processors, compute local 
!! sizes and allocate arrays
!    
      cmi=cmplx(0,1); pi=acos(-1.0d0) 

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

      do j=1,nr
        zk(j) = (j-1)*dr     
      enddo
      !b=0.1+abs(10*(zk-p))
      b=0.1+7*(abs(zk-p))**(1.75)
      matarg=(-1.0d0/sqrt((zk-p)**2+b))*(-cmi*dt)
      write(81,*)zk
      write(81,*)-1.0d0/sqrt((zk-p)**2+b)
      dk=2*pi/(nr*dr)
!! cft from tor (2003)
!      do j=1,nr/2+1
!       pk(j)=(j-1)*dk
!      enddo
!      do j=nr/2+2,nr
!        pk(j)=(-nr+j-1)*dk
!      enddo
! cft ala Hicham:
      do j=0,nr/2-1
        pk(j+1) =  j*dk 
       enddo
       do j=nr/2,nr-1
	 pk(j+1) = (j-nr)*dk
       enddo

      call fft_init(nr, psi, planfwd, planbwd)  

      pkexp=exp(cmplx(0.d0,-dt/2*pk*pk))

      open(10,file=psi_infile)   
        do i = 1,ninit-1            
          read(10,*)x
        enddo
        read(10,*)inpsi
      close(10)
      psi=inpsi

! OLD absorbator close to both edges (0 and rmax):
!      do ii = 1,nr
!        absorb(ii)=1.0d0
! rmax:
!        if (abs(zk(ii)) .ge. 0.9d0*rmax) then
!	    absorb(ii)=(dcos(abs(abs(zk(ii))-0.9d0*rmax)/(0.1d0*rmax)*pi/2.d0))**(1.d0/8.d0)
!        endif
! 0:
!         if (abs(zk(ii)) .le. 0.1d0*rmax) then
!	    absorb(ii)=(dcos(abs(abs(zk(ii))-0.1d0*rmax)/(0.1d0*rmax)*pi/2.d0))**(1.d0/8.d0)
!        endif       
!      enddo

! New absorbator:
      zmax=rmax/2; zlim=0.95*zmax;
      do ii = 1,nr
        if (abs(zk(ii)-p) .ge. zlim) then
          absorb(ii)=absorber(dabs(zk(ii)-p),zlim,zmax)
        else
          absorb(ii)=1.0d0   
        endif
      enddo    
!
! Initiate:
!
     nsteps=Tprop/dt                
     write(*,*)'Starting: nr =    ',nr, ' nsteps = ',nsteps
!
! Set up initial fl:
! 
     t=Tstart; nsteps=nsteps-1; tprintnext=Tstart+Tprop/nut;

      open(9,file=outf(1:outlen)//'_z.dat'); write(9,*) zk(1:nr); close(9)
      open(10,file=outf(1:outlen)//'_T.dat')
      open(11,file=outf(1:outlen)//'_psiSQ.dat')
      open(12,file=outf(1:outlen)//'_psiR.dat')
      open(13,file=outf(1:outlen)//'_psiI.dat')
     
      call print(t,nr,psi,outf,outlen)  
!
! 1. Iterate first step dt/2 in r-direction:
!                                           
! isign=1 is flagging "Forward" transform, -1 the inverse:
!     isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*exp(cmplx(0.d0,-dt*pk*pk/2.0d0))
!     isign = -1; call fft(nr,isign,cxy); psi=cxy;
!
     call dfftw_execute(planfwd); psi=psi*exp(cmplx(0.d0,-dt*pk*pk/2.0d0))
     call dfftw_execute(planbwd); psi=psi/real(nr)

!
!--------------------------------------------------------------
! main loop: Timepropagation from t=0 to t=(n_steps-1)*dt:
!--------------------------------------------------------------
!
      call cpu_time(start)     
      do kk=1,nsteps
 !write(*,*)kk, t, tprintnext

        t=t+dt/2
! Velocity gauge:
        if (t .lt. 0) then
           ft = -e0/wL*(cos(wL*t)+1)
        else
           ft = -e0/wL*(cos(1.5*wL*t)+1) 
        endif


        psi=psi*exp(matarg)*absorb
       
!        isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*pkexp
!        isign = -1; call fft(nr,isign,cxy); psi=cxy;
        call dfftw_execute(planfwd); psi=psi* pkexp * exp(cmplx(0.d0,dt*ft*pk))
        call dfftw_execute(planbwd); psi=psi/real(nr)

        t=t+dt/2    
        if  (t .gt. tprintnext) then
          write(*,'(2F13.3,F8.2,A)')t,sum(abs(psi)**2)*(zk(2)-zk(1)), real(kk)/nsteps*100, ' % finnished'
          open(95,file='monitor_time.txt')
          write(95,'(2F13.3,F8.2,A)')t,sum(abs(psi)**2)*(zk(2)-zk(1)), real(kk)/nsteps*100, ' % finnished'
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
     ft = -e0/wL*(cos(1.5*wL*t)+1) 
     
     psi=psi*exp(matarg)

! isign=1 is flagging "Forward" transform, -1 the inverse:
!     isign = 1; cxy=psi; call fft(nr,isign,cxy); cxy=cxy*exp(cmplx(0.d0,-dt*pk*pk/2.0d0))
!     isign = -1; call fft(nr,isign,cxy); psi=cxy;

     call dfftw_execute(planfwd); psi=psi*exp(cmplx(0.d0,-dt*pk*pk/2.0d0))*exp(cmplx(0.d0,dt*ft*pk/2.0d0))
     call dfftw_execute(planbwd); psi=psi/real(nr)
     
     call print(t,nr,psi,outf,outlen)

        close(10)
        close(11)
	close(12); close(13)

     stop
     end


      subroutine print (t, nr, psi, outf,outlen)          
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

!	write(*,*)'leaving print:'

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
!
      integer, parameter                                   :: dbl = selected_real_kind(p=13)
      integer, intent(in)                                  :: nr
      integer :: nzp
      double complex,  dimension(nr), intent(inout)        :: psi
      integer*8,intent(out)                                :: planfwd
      integer*8,intent(out)                                :: planbwd
!integer*8 :: FFTW_FORWARD, FFTW_MEASURE, FFTW_BACKWARD
!
!call dfftw_plan_many_r2r(planfwd,&         ! 0. Transform pointer (integer*8)
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
!						 
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
!
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
