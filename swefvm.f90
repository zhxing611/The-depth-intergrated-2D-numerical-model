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

	!calculation
	do while (t <= tstop) 
		write(*,*) 't = ', t, ' s out of ', tstop, 's' 
		nStep=nStep+1

		call swe		
		call results_output

        if(mod(nstep,100).eq.0)write(11,*) t, Q1(15473) !P1
        if(mod(nstep,100).eq.0)write(12,*) t, Q1(13018)!P2
		if(mod(nstep,100).eq.0)write(13,*) t, Q1(12942)!P3
        if(mod(nstep,100).eq.0)write(14,*) t, Q1(12867)!P4
		if(mod(nstep,100).eq.0)write(15,*) t, Q1(11071)!P5
        if(mod(nstep,100).eq.0)write(16,*) t, Q1(15247)!P6
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
