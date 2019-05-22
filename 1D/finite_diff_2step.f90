      PROGRAM PSIrthet
!
! Finite difference 1D TDSE solver, (PRA 77 014701)
! 
! Psi(t+dt) = Psi(t-dt) - 2 i H dt Psi(t)
! Variables: Psi_new=Psi(t+dt); Psi_old=Psi(t-dt); Psi=Psi(t)
! 
      implicit none

      integer                   :: nsteps,nsteps2           ! Number of timestep      
      integer                   :: nr ! Problem size

      integer i,j,k,kk,itype,ii, isign
      double precision t,dr,dt,Tprop, Tstart, dk     

! Input variables:
      integer outlen,  nut , ninit
      integer			:: status
      character*80 outf, psi_infile
      double precision tprintnext, wL, e0, rmax , b, p, pi, powerof2, x
      
! wavefunction:
      double complex,   dimension(:),allocatable 	:: psi_Even, psi_Odd, Dl, Du, Dgnal
      double precision, dimension(:),allocatable	:: absorb, inpsi, Vt
      double precision, dimension(:),allocatable	:: zk
      
! Rydberg ting. Puls er e0*sin^2(alfa t)*cos(wLt)
      double precision ft, dzinv2, start, finish
      complex cmi, cmi2dt
!
!! Read input, broadcast it to all processors, compute local 
!! sizes and allocate arrays
!    
      cmi=cmplx(0,1); pi=acos(-1.0d0) 

      read (*,*) powerof2, rmax, ninit
      read (*,'(A)') psi_infile
      read (*,*) b, p	! b is coulomb softening factor and p is position of coulomb pot from matlab
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
      allocate (inpsi(nr), STAT=status); allocate (Vt(nr), STAT=status);
      allocate (psi_Odd(nr), STAT=status); allocate (psi_Even(nr), STAT=status)
      allocate (absorb(nr), STAT=status); allocate (zk(nr), STAT=status); 
      allocate (Dl(nr), STAT=status); allocate (Du(nr), STAT=status); allocate (Dgnal(nr), STAT=status); 

! z in -rmax/2,rmax/2:
      do j=1,nr
        zk(j) = (j-1)*dr-p     
      enddo
      Vt=(-1.0d0/sqrt(zk**2+b))
    
      open(10,file=psi_infile)   
        do i = 1,ninit-1            
          read(10,*)x
        enddo
        read(10,*)inpsi
      close(10)
      psi_Even=inpsi

! absorbator close to both edges (0 and rmax):
      do ii = 1,nr
        absorb(ii)=1.0d0
! rmax:
        if (abs(zk(ii)) .ge. 0.8d0*(rmax-p)) then
	    absorb(ii)=(dcos(abs(abs(zk(ii))-0.8d0*(rmax-p))/(0.2d0*(rmax-p))*pi/2.d0))**(1.d0/8.d0)
        endif
! 0:
!         if (abs(zk(ii)) .le. 0.1d0*rmax) then
!	    absorb(ii)=(dcos(abs(abs(zk(ii))-0.1d0*rmax)/(0.1d0*rmax)*pi/2.d0))**(1.d0/8.d0)
!        endif       
      enddo    

!
! Initiate:
!
      nsteps2=Tprop/dt/2; nsteps=2*nsteps2                
      write(*,*)'Starting: nr =    ',nr, ' nsteps = ',nsteps
        
!
! Set up initial:
! 
      t=Tstart; nsteps=nsteps-1; tprintnext=Tstart+Tprop/nut;
      psi_Odd=Psi_Even; Dl=0; Du=0
      dzinv2=-0.5/dr**2; cmi2dt=cmi*dt*2.0d0; 

      write(67,*)zk
      write(67,*)Vt
      write(67,*)absorb
     

      open(9,file=outf(1:outlen)//'_z.dat'); write(9,*) zk(1:nr); close(9)
      open(10,file=outf(1:outlen)//'_T.dat')
      open(11,file=outf(1:outlen)//'_psiSQ.dat')
      open(12,file=outf(1:outlen)//'_psiR.dat')
      open(13,file=outf(1:outlen)//'_psiI.dat')
     
      call print(t,nr,psi_Even,outf,outlen)  
!
!--------------------------------------------------------------
! main loop: Timepropagation from t=0 to t=(n_steps-1)*dt:
!--------------------------------------------------------------
!
     call cpu_time(start)     
     do kk=1,nsteps2
        t=t+dt
! Length gauge:
        if (t .lt. 0) then
           ft = -e0*sin(wL*t)
        else
           ft = -1.5*e0*sin(1.5*wL*t)
        endif       
 
! Two step formula: Psi(t+dt) = Psi(t-dt) - 2*i * H * dt * Psi(t)
           
        Dl(2:nr)=psi_Odd(1:nr-1)*dzinv2*cmi2dt 
        Du(1:nr-1)=psi_Odd(2:nr)*dzinv2*cmi2dt      
        Dgnal=(-2.0d0*dzinv2+Vt+ft*zk)*psi_Odd*cmi2dt
        psi_Even = (psi_Even -(Dl+Du+Dgnal)) *absorb

        t=t+dt
! Length gauge:
        if (t .lt. 0) then
           ft = -e0*sin(wL*t)
        else
           ft = -1.5*e0*sin(1.5*wL*t)
        endif
        
        Dl(2:nr)=psi_Even(1:nr-1)*dzinv2*cmi2dt 
        Du(1:nr-1)=psi_Even(2:nr)*dzinv2*cmi2dt      
        Dgnal=(-2.0d0*dzinv2+Vt+ft*zk)*psi_Even*cmi2dt
        psi_Odd = (psi_Odd - (Dl+Du+Dgnal)) *absorb
  
    
        if  (t .gt. tprintnext) then
          write(*,'(F13.3,F8.2,A)')t, real(2*kk)/nsteps*100, ' % finnished'
          open(95,file='monitor_time.txt')
          write(95,'(F13.3,F8.2,A)')t, real(2*kk)/nsteps*100, ' % finnished'
          close(95)
          tprintnext=tprintnext+Tprop/nut
          call print(t,nr,psi_Odd,outf,outlen)
        endif

     enddo
     call cpu_time(finish)
     write(*,*)'Time: ',finish-start
!--------------------------------------------------------------
! END Main loop. What is left is the  final steps & printout:
!--------------------------------------------------------------
          
     call print(t,nr,psi_Odd,outf,outlen)

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


        
          write(10,'(5000000F24.9)')t
          write(11,'(5000000F16.9)')(abs(psi(j))**2,j=1,nr)
          write(12,'(5000000F16.9)')(real(psi(j)),j=1,nr)
          write(13,'(5000000F16.9)')(aimag(psi(j)),j=1,nr)

!	write(*,*)'leaving print:'

      return
      end


