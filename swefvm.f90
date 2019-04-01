!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c Finite Volume Method for the 2D Shallow Water Equation
!   c with Sediment Transport on Unstructured Mesh
!   c (the sediment transport part is solved seperately with 
!   c the flow.)
!   c by Xiaofeng Liu, Hydrosystems Lab, UIUC, 10/26/2008
!   C liu19@illinois.edu

    program main
	USE COMMON_MODULE
    implicit none
	real*8 timestart,timestop,timetotal

	
	write(*,*) '-------------------------------------------------------------------'
    write(*,*) '|   Shallow    | HydroSed2D: Shallow Water Equation with Sediment |'
	write(*,*) '|   Water      |                                                  |'
	write(*,*) '|   Equation   |    by     : Xiaofeng Liu, Hydrosystems Lab, UIUC |'
	write(*,*) '|   with       |                                                  |'
	write(*,*) '|   Sediment   |    version: 1.0  (12/5/2008)                     |' 
	write(*,*) '-------------------------------------------------------------------'
    write(*,*) 'Web  : http://vtchl.uiuc.edu/~liu19     '
	write(*,*) 'Email: liu19@illinois.edu'
	write(*,*) '--------------------------------------------------------------------'
    write(*,*) '|Base on HydroSed2D£¬Mingliang Zhang and Hongxing Zhang further     |'
	write(*,*) '|developed the depth-averaged 2D hydrodynamic model by introducing  |'
    write(*,*) '|treatment technology of wet-dry boundary                           |'                                                                                                          
	write(*,*) '|Mingliang Zhang (zhmliang_mail@126.com);                           |'
	write(*,*) '|Hongxing Zhang (zhxing611@163.com)                                 |'
	write(*,*) '|School of Ocean Science and Environment, Dalian Ocean University   |'
	write(*,*) '|--------------------------------------------------------------------'
	write(*,*) ' '
	
	call cpu_time(timestart)

	!initialization
	call init

	call results_output


!-------------read input water level data
	  open(2,FILE='input/h.dat',STATUS='UNKNOWN')
	    do aa=1,nDEMPoints
   	     read(2,*) wse(aa,1),wse(aa,2)
		 enddo

	  close(2)
     
	!calculation
	do while (t <= tstop) 
		write(*,*) 't = ', t, ' s out of ', tstop, 's' 
		nStep=nStep+1

		call swe
	    	

		call results_output
    
	 if(mod(nstep,250).eq.0)write(111,*) t, eta(22554)  
	 if(mod(nstep,250).eq.0)write(112,*) t, eta(3225)  
	 if(mod(nstep,250).eq.0)write(113,*) t, eta(5426)

		t=t+dt

		tscount=tscount+1
		if(tscount==sedInterval)then
			scalc=1	
			tscount=0
		else
			scalc=0
		end if

    end do

	call cpu_time(timestop)

	timetotal = timestop - timestart

	write(*,*) 'Total CPU time = ', timetotal, ' seconds.'

	end
