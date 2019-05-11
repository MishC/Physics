program  Hmatrix_Bsplines
!i -angular number (number of channels), N- number of states per channel, [a,b] -box size
!N_b-number of B-splines per channel, n -order of B-splines, bp -how many breakpoints; 
            use parameters
			use sorting_module
			use eigreal
            implicit none
					integer ::  i, j, k, kk,hh, ind1,ind, ind2,bp,nx, s,numb !indexes, size, etc.
                    integer, parameter :: n=5, lmin=0, lmax=3, m=10, alfa=0, omega=100, Nb=250
                    character(string_length)::format 
                    !!!!!!!!!!!!!!Space parametres inside box
					real(kind = dp) 	  	            :: dX !number of points break points
                    real(kind = dp), dimension(:), allocatable :: x_break !<alfa:dX:omega> 
                    real(kind=dp)                              ::lfac                    
                              
                    !!!!!!!!!!!!!!For Gaus-LEgendre integration  :::Golub Welsch algorithm    
                    real(kind = dp),            DIMENSION(n)                :: B !Helping number for G-L
                    REAL(kind = dp),            DIMENSION(n,n)              :: JJ !helping matrix in G-L integration
                    real(kind = dp),            DIMENSION(n)                :: xx, WW ! 
                    REAL(kind = dp),            DIMENSION(:),   ALLOCATABLE :: x ! Gauss Legendre integration points for B-splines
                    REAL(kind = dp),            DIMENSION(:),   ALLOCATABLE :: W ! Weights inside integration points for B-splines
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Parameters for dsyev!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					INTEGER								:: Ma, Na
					INTEGER								:: LDA
					INTEGER								:: LWORK
                    REAL(kind = dp),            DIMENSION(:),   ALLocatable :: WORK                 
                    INTEGER :: info
 
					REAL(kind = dp),			DIMENSION(:),		ALLOCATABLE		:: E, E1,E2
					REAL(kind = dp),			DIMENSION(n,n)				:: V
                    LOGICAL                                                  :: dovecs
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!B-splines parametres!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                    real(kind = dp),			dimension(:,:,:),	ALLOCATABLE:: Bsp
		     		real(kind = dp),			dimension(:,:,:),	ALLOCATABLE	:: diffBsp			
					
				    real(kind = dp)             :: AA,BB, CC, DD, EE !PArametres for Jordan-Gordon integration    

                    REAL(kind = dp),            DIMENSION(:),   ALLOCATABLE :: t    
                    REAL(kind = dp),            DIMENSION(:),   ALLOCATABLE :: yy1, yy2, yy3 

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! matrices and parametres!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    REAL(kind = dp),            DIMENSION(:,:), ALLOCATABLE :: H, overlap,overlap1
                    INTEGER :: sh, s_x,l
                    REAL(kind=dp) ::  int1, int2
                    REAL(kind = dp),            DIMENSION(:,:), ALLOCATABLE :: psi1, psi2
                    REAL(kind = dp),            DIMENSION(:,:), ALLOCATABLE :: BIGPSI
 
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EIGENVALUE:::PROBLEM:::!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                     REAL(kind = dp),            DIMENSION(:), ALLOCATABLE :: En,Erg,En2,Erg2
                     REAL(kind = dp),            DIMENSION(:,:), ALLOCATABLE :: Vec,Vv
                     ! REAL(kind = dp),            DIMENSION(:,:), ALLOCATABLE ::  Vector
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Coupling::Elements!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      REAL(kind = dp),            DIMENSION(:),   ALLocatable :: kopdip1, kopdip2 
                      REAL(dp), DIMENSION(:), Allocatable :: kopdip11
                     ! REAL(kind = dp),            DIMENSION(:),   ALLocatable :: ind_L1, ind_L2,ind_R1,ind_R2
                      INTEGER :: cup 
                      integer, dimension(:), allocatable :: icod, indl1,indr1, indl2, indr2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!						
						kk=0
                        hh=0
                                           
                        !dX- break point interval, x_break- number of break points 
						dX=(omega-alfa)/(Nb-1.)
						bp=floor(omega/dX)+2
						!print*, bp
						allocate(x_break(bp))
						x_break(1)=alfa;
						do i=2,bp
						  x_break(i)=x_break(i-1)+dX
                        end do
                         !print *, xbreak(bp)

                        open(unit=1,file='space_parametres.txt', status='replace')	
                        write(1,*) "No. of break points:", bp
                        write(1,*) "Size of space", size(x_break(:))
                        write(1,*) "first three points", x_break(1), x_break(2), x_break(3)
                        write(1,*) "Last point", x_break(bp)
                        close(1)
						
                        
                        do i=1,n
					
					        B(i)=0.5/sqrt(1.d0-((2.d0*i)**(-2.d0)))
				  		
		                end do
		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gauss-legendre integration, Golub-Welsch algorithm!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						
	            !!!!!!JACOBI MATRIX!!!!!
                
      			do j=1,n
	    			do i=1,n
	  			
						JJ(i,j)=0
							
		  			end do
			   end do
			
					do i=1,n-1
							JJ(i,i+1)=B(i)
							JJ(i+1,i)=B(i)				
					end do
					

				!DIAGONALIZATION START HERE


				Ma = size(JJ(:,1))
				Na = size(JJ(1,:))
				

			  !	ALLOCATE(V(Ma,Na))
				ALLOCATE(E(Ma))
				
                
				V(:,:)= JJ(:,:)
					
				
				LDA = Ma

				!Calculating efficient worksize by passing -1 to LAPACK algorithm
				LWORK =  -1
				ALLOCATE(WORK(1))
				!WRITE(*,*) 'CALL ONE!', Ma, LDA, WORK, LWORK
				call dsyev('V', 'L', Ma,  V , LDA, E, WORK, LWORK, info)
				!WRITE(*,*) 'info on call One: ', info
				
				!Reallocate workspace to optimal size
			  LWORK = INT(WORK(1))

 			  DEALLOCATE(WORK)
 			  ALLOCATE(WORK(LWORK))

				!WRITE(*,*) 'CALL TWO!', Ma, LDA, SIZE(WORK), LWORK
				!Call the Actual Eigenvalue algorithm
				V(:,:) = JJ(:, :)
				call dsyev('V', 'L', Ma,  V, LDA, E,  WORK, LWORK, info)
				!WRITE(*,*) 'info on call One: ', info
				!DIAGONALIZATION ENDS HERE
					
			   DEALLOCATE(WORK)		
										
					do j=1,n	
					  WW(j)= 2.*V(1,j)**2.
					  xx(j)=E(j)				
					end do	


			ALLOCATE(x(Nb*n-n))
			ALLOCATE(W(Nb*n-n))
				
			
			do i=2,Nb

				x((i-2)*n+1:(i-1)*n)=(x_break(i)-x_break(i-1))/2.d0*xx(1:n)+(x_break(i)+x_break(i-1))/2.d0
				W((i-2)*n+1:(i-1)*n)=WW(1:n);

			end do
							
		    open(unit=2,file='firstcheck',status='replace')   
			
                write(2,*) "Recurrence coefficient", B(:)
                write(2,*) "Jacobi matrix size", Ma, Na
			    write(2,*) "Number of nods", size(x) 
                write(2,*) "First 4 points", x(1),x(2),x(3),x(4)
                write(2,*) "Number of weights", size(W)
                write(2,*) "Max limit for space", max(x,1.)
                write(2,*) "Max of weights", max(W,1.)
                
			
             close(2) 
			

					nx=n
   			    	s=size(x)
			
   
		Allocate(Bsp(1:s,1:(Nb-1+2*n),1:n+1))
         Bsp(:,:,:)=0.d0

        Allocate(t(1:Nb+2*n))
		t(1:Nb+2*n)=0
 		do i=1,n
 			t(i)=0
 		end do	
        
 		do i=(n+1),(Nb+n-1)
             Bsp((i-n-1)*nx+1:(i-n)*nx,i,1) = 1.d0
 			t(i)=(i-n-1)*dX
 		end do
          
 	do i= (Nb+n),(Nb+2*n)
 			t(i)=omega
 		end do
 		
 		Allocate(yy1(1:s))	
 		Allocate(yy2(1:s))	
 		!	yy1=0
 		!	yy2=0	
        
!            print *, 'Allocated'

 		do k=1,n
 				do i=1, (Nb-1+2*n-k)
 						
 						yy1=Bsp(1:s,i,k)
 						yy2=Bsp(1:s,i+1,k)
 						!print *,yy1	
 			AA=t(i+k)-t(i)
            
 			BB=t(i+k+1)-t(i+1)
 			
 				if (AA==0) then
     				AA=1
 				end if    
 				
 				if (BB==0) then
     				BB=1
 				end if 
 		!!	ALLOCATE(x(Nb*n-n))
 			Bsp(1:s,i,k+1)= ((x(:)-t(i)) / AA ) * yy1(:) + ( (t(i+k+1)-x(:) ) /(BB) ) *yy2(:)
 			
            end do ! (i)
            end do !(k)
!format=(f3.6)"
            
		   open(unit=3, file="Bsplines.dat", status='replace')
				
           !do k=1,n+1
              !do j=1,Nb-1+2*n
              
                do i=1,5
			
                  write(3,*)  Bsp(i,1:5,n+1) 
			  
                 end do
            
            close(3)
!            print *,'toto jemaximum',maxloc(Bsp)
!			gnuplot  -> plot 'Bsplines1.dat'
			
			Allocate(diffBsp(1:s,1:Nb-1+2*n,1:n+1))
			diffBsp(1:s,1:Nb-1+2*n,1:n+1)=0.d0;
				

				deallocate(yy1)
				deallocate(yy2)
				
        		allocate(yy1(1:s))
				allocate(yy2(1:s))
				allocate(yy3(1:s))
			     
                
				do k=2,n
        				 do i=1,(Nb-1+2*n-k) 
			
				
				yy1=Bsp(1:s,i,k-1)
                
				yy2=Bsp(1:s,i+1,k-1)

                !if (k==5 .and. i==250) then
               !do ind2=1220,1245
                !  print*, '\n', yy2(ind2)
                !end do
              !end if
				yy3=Bsp(1:s,i+2,k-1)
					
				AA=t(i+k)-t(i)
					if (abs(AA) <tol) then
              			AA=1.d0
					end if
   		
    			BB=t(i+k+1)-t(i+1)    
      				if (abs(BB)<tol) then
                        BB=1.d0
  					end if
  
                CC=t(i+k-1)-t(i)  
					if (abs(CC) < tol) then
    			      CC=1.d0
                    end if

                DD=t(i+k)-t(i+1)    
					  if (abs(DD) < tol) then
              			DD=1.d0
                      end if

                EE=t(i+k+1)-t(i+2)    
                  if (abs(EE) < tol) then
                       EE=1.d0
            		end if
			
                diffBsp(1:s,i,k+1)=(k*(k-1)/AA)*(yy1(1:s)/CC-yy2(1:s)/DD)-(k*(k-1)/BB)*(yy2(1:s)/DD-yy3(1:s)/EE)

              end do
            end do
			
			open(4, file="Diffbsp.dat", status="replace")
            
            !do k=1,n
                !do j=1,Nb-1+2*n-k
                  do ind=1,5
					write(4,*)  diffBsp(ind,1:5,n+1)
				end do
             ! end do
            !end do

            close(4)
			

  !!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%Calculation of matrix elements of Hamiltonian: TDSE Hydrogen potential%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(H(1:Nb-3+n,1:Nb-3+n))
  ALLOCATE(overlap(1:Nb-3+n,1:Nb-3+n))
   H(1:Nb-3+n,1:Nb-3+n)=0.d0
   overlap(1:Nb-3+n,1:Nb-3+n)=0.d0
   sh=size(H,1)
    Allocate(psi1(1:s,1:m))
     psi1(1:s,1:m)=0.d0
     ! Allocate(indl1())
    !  Allocate(indl2
     !do i = 1, N
     !  psi = - a 
   Allocate(psi2(1:s,1:m))
      psi2(1:s,1:m)=0.d0
  cup=lmax*m*m
  numb=(lmax+1)*m
  allocate(Bigpsi(s, m*lmax))
  bigpsi(:,:)=0.d0
  Allocate(kopdip1(cup))
  kopdip1(:)=0.d0
  Allocate(kopdip2(cup))
  !kopdip2(:)=0.d0
  Allocate(kopdip11(numb))
  !kopdip11=0.d0
    Allocate(indl2(numb))
    Allocate(indr2(numb))
    Allocate(indl1(cup))
    Allocate(indr1(cup))
    Allocate(icod(cup))

    do l=lmin, lmax-1
        print *, l
          lfac=(l+1)**2.d0
          lfac=lfac/(4*(l+1)**(2.d0)-1)
          lfac=sqrt(lfac)         
        print *, 'lfak',lfac
          if (l==0) then
          
          do i=1,sh
            do j=i,min(i+n,sh)
            !!Kinetic Energy!!
            int1= -.25d0*dX*sum(W(1:s)*Bsp(1:s,i+1,n+1)*diffBsp(1:s,j+1,n+1))
            H(i,j)=H(i,j)+int1

            !!Potential Energy!!
            int1=.5d0*dX*sum(W(1:s)*Bsp(1:s,i+1,n+1)*Bsp(1:s,j+1,n+1)*(-1.d0/x(1:s) + (l*(l+1.d0)/(2.d0*x(1:s)**2))))
            H(i,j)=H(i,j)+int1
            !!Overlap Matrix!!
            
            int2=dX*.5d0*sum(W(1:s)*Bsp(1:s,i+1,n+1)*Bsp(1:s,j+1,n+1))
            overlap(i,j)=overlap(i,j)+int2  
          
            if (j .ne. i) then
               H(j,i)=H(i,j)
               overlap(j,i)=overlap(i,j)
            end if
          
          end do !j
          end do !i
        
          open (11, file='hamiltonian.dat', status='replace')
          
          do i=1,4
              write(11,*) H(i,1:4)
          end do

          close(11)

		  		Ma = size(H(:,1))
				Na = size(H(1,:))
				
                
                   open(10, file='overlap.dat', status='replace')
                 do ind2=1,5 
                   write(10,*) overlap(ind2,1:5)
                 end do
                   close(10)
           
                   
                 ALLOCATE(Vec(Ma,Ma))
                 ALLOCATE(En(Ma))
				LDA = Ma
                 dovecs=.true.
                 call realgeneig(Ma,dovecs,H,overlap,En,Vec,info)

				!DIAGONALIZATION ENDS HERE
         		 ! DEALLOCATE(WORK)
                !    write(*,*), En(1:5),'\n'        
                 Allocate(Erg(Ma))
               
                 !!call geneig(Ma,dovecs,H,S,En,Vec,info)
                 Erg=En    
                call sorting(Ma,Erg)
                open(12,file='lzerostavy', status='replace')
                do i=1,Ma
                write(12,*) Erg(i),'\n'
              end do
            close(12)
            open(14,file='lzerovektory',status='replace')
               
                do i=1,5
                 write(14,*) Vec(i,1:5)
               end do
           
           close(14)
        
           write(*,*) 'Size of space', s
           write(*,*) 'No. of points inside l', m

      !!ALLOCATE(psi1(s,m))
        !!          psi1(:,:)=0.0d0;
            
           do i=1,m
              do j=1,sh
                    !!PRINT *, SHAPE(psi1), SHAPE(B) 
                    psi1(:,i)=psi1(:,i)+Vec(j,i)*Bsp(:,j+1,n+1)
                 !         print *, psi1(:,1) 
                end do
              
                if (psi1(2,i)<0) then
                  psi1(:,i)=-psi1(:,i)
                end if
              
              end do !!1...N
 
             
          ALLOCATE(E1(ma))
              E1=Erg(1:ma)       
          !Allocate(VV(Ma,Ma)) 
          !Allocate(E2(ma))
       end if  !l==lmin!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
          
       if (Allocated(H) .eqv. .false.) then
      Allocate(H(1:Nb-3+n,1:Nb-3+n))
       end if

       if (Allocated(overlap) .eqv. .false.) then
         Allocate (overlap(1:Nb-3+n,1:Nb-3+n))
       end if
       
       H(1:Nb-3+n,1:Nb-3+n)=0.d0
       overlap(1:Nb-3+n,1:Nb-3+n)=0.d0
       sh=size(H,1)
                 !if (Allocated(VV) .eqv. .true.) then
                  !    deAllocate(VV(Ma,Ma))
                 !end if
                !Allocate vv(Ma,Ma) 
       do i=1,sh
         do j=i,min(i+n,sh)
         
    !!Kinetic Energy!!
            int1= -.25d0*dX*sum(W(1:s)*Bsp(1:s,i+1,n+1)*diffBsp(1:s,j+1,n+1))
            H(i,j)=H(i,j)+int1

            !!Potential Energy!!
            int1=.5d0*dX*sum(W(1:s)*Bsp(1:s,i+1,n+1)*Bsp(1:s,j+1,n+1)*(-1.d0/x(1:s) + ((l+1)*(l+2.d0)/(2.d0*x(1:s)**2))))
            H(i,j)=H(i,j)+int1
            !!Overlap Matrix!!
            
            int2=dX*.5d0*sum(W(1:s)*Bsp(1:s,i+1,n+1)*Bsp(1:s,j+1,n+1))
            overlap(i,j)=overlap(i,j)+int2  
          
            if (j .ne. i) then
               H(j,i)=H(i,j)
               overlap(j,i)=overlap(i,j)
            end if
          
         end do !j
     end do !i
        
          open (15, file='h2.dat', status='replace')
          
          do i=1,4
              write(15,*) H(i,1:4)
          end do

          close(15)

		  		Ma = size(H(:,1))
				Na = size(H(1,:))
			        	
                
                   open(16, file='S2.dat', status='replace')
                 do ind2=1,5 
                   write(16,*) overlap(ind2,1:5)
                 end do
                   close(16)
              
                   !!!FORTRAN DO NOT LIKE THIS 
                      if (allocated(En2)  .eqv. .false.)then
                     allocate(En2(MA))
                 end if
                
                 if (Allocated(VV) .eqv. .false.) then
                   Allocate(VV(Ma,Ma))
                 end if


                 LDA = Ma
                 dovecs=.true.
                 call realgeneig(Ma,dovecs,H,overlap,En2,Vv,info)
              
           do i=1,m
                do j=1,sh
              
                    psi2(1:s,i)=psi2(1:s,i)+Vv(j,i)*Bsp(1:s,j+1,n+1)
  
                end do
              
                if (psi2(2,i)<0) then
                  psi2(:,i)=-psi2(:,i)
                end if
              
           end do !!1...N


           if (ALLOCATED(E2) .eqv. .false.) then
          Allocate(E2(Ma))
         end if

             
              E2(:)=En2(:)
        
              !!!!!!!Couplings with EMF!!!!  
         
!            if (allocated(kopdip1) .eqv. .false.)
!          Allocate(kopdip1(cup))
!        end if 
!
!        if (allocated(kopdip2) .eqv. .false.)
!          Allocate(kopdip2(cup))
!        end if

          do j=1,m
            do i=1,m

              hh=hh+1
              icod(hh)=1
              indl1(hh)=l*m+j
              indr1(hh)=(l+1)*m+i 
              
              !!!Velocity gauge
                 kopdip1(hh)= dX*0.5d0*sum(W(1:s)*psi1(1:s,j)*psi2(1:s,i)*x(1:s)) *(En2(i)-E1(j))*Lfac
                 !PRINT *, KOPDIP1(hh)
              !!!Length gauge
                 kopdip2(hh)=dX*0.5d0*sum(W(1:s)*psi1(1:s,j)*psi2(1:s,i)*x(1:s))*Lfac
                 !print*, kopdip2
            end do
          end do
               
      !  if (Allocated(kopdip11) .eqv . .false.)
      !          allocate(kopdip11(numb))
     ! end if
     
      if (l==lmin) then
          do j=1,m
            kk=kk+1
            indl2(kk)=(l+1)*m+j
            indr2(kk)=(l+1)*m+j
            kopdip11(kk)=E1(j)

          end do
      end if
       
      
      do i=1,m
          kk=kk+1
          indl2(kk)=(l+1)*m+i
          indr2(kk)=(l+1)*m+i
          kopdip11(kk)=E2(i)
      end do

                 E1=E2
                 !!Allocate(Bigpsi(s,l*m)
                 Bigpsi(1:s,1+l*m:(1+l)*m)=psi1
                 psi1=psi2

                 write(*,*)'Size of indexu %d', size(indl1)
      end do !for l=lmin:lmax 
      k=0
      ind=0
!      numb=(lmax+1)*m
!      allocate(kopdip11(numb))
!      kopdip11(:)=0
!      
!   do l=lmin,(lmax-1)
!       if (l==lmin) then
!       do j=1:N
!           
!       ind=ind+1
!       indl2(ind)=l*m+j
!       indr2(ind)=l*m+j
!       kopdip11(ind)=E1(j)
!       end do
!       end if
!      do i=1:m
!        k=k+1
!        indl2(k)=(l+1)*m+i
!        indr2(k)=(l+1)*m+i
!        kopdip11(k)=E2(i)
!    end do
!      end do

 open(unit=5,file="ENERGY",status="replace")
 do i=1,numb
 write(5,*), kopdip11(i)
 end do
 close(5)

 open(unit=6,file="COUPLINGS",status="replace")
 !!format='(i5, i5, e10.5, "\n")'
 do i=1,cup
   !do j=1,sh
   write(6,*), indl1(i),indr1(i), kopdip1(i) 
 !end do
 end do
 close(6)
            
 	  end program Hmatrix_Bsplines


 
 
