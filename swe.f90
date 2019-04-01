!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    HydroSed2D Copyright (C) 2008 Xiaofeng Liu
!
!    License
!
!    This file is part of HydroSed2D.
!
!    HydroSed2D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HydroSed2D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HydroSed2D.  If not, see <http://www.gnu.org/licenses/>.
!
!    Base on HydroSed2D, Mingliang Zhang and Hongxing Zhang further developed the depth-averaged 2D hydrodynamic model 
!    by introducing treatment technology of wet-dry boundary and considering vegetation effects. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   Shallow water equation: one step
	subroutine swe()
	USE COMMON_MODULE
	implicit none
    
	integer i,j,k,mm,o,w
	integer*4 dryindex,wetindex
	real*8  Q1FN,EtaFN
	character*20 temp
	character*14 tem
    integer NoUse	

    if(inctr==3) then 

		 call inputcells_q_calc()               
	end if
    !estimate the wet or dry for each edge
    call wetdry_edge

    !chang eta
   !call change_eta
   ! call cellCenterToNodes(eta,nodeZsurf)
   ! call change_z
   !! call updateMeshData

	!calculate the gradient of Q at cell centers
	call gradientQ

	!calculate the gradient of u and v at cell centers
	call gradientUV
	call edgeGradientUV

	!calculate the gradient of eta at cell centers
	call gradientEta

	!calculate the face limiters
	call calc_faceLimiters

!******************************************
	do i=1,nFaces

	   wetindex=0
	   dryindex=0
       do j=1,faceEdgesNum(i) 	
          curEdge = faceEdges(i, j)	      
 	      wetindex=wetindex+edgedrywet(curEdge)
       enddo
      
	  if(wetindex.gt.3)then
		 calc_flag(i)=1	
	  elseif(wetindex==0)then
	     calc_flag(i)=0
	  else
	     calc_flag(i)=2			 
	 endif
    
   enddo 
!****************************************************
	do  i=1,nFaces

     curFace=i     
	   if(calc_flag(i)==1.or.calc_flag(i)==2)then
   					
			call border_values  
	   endif


		if(Q1(i)<=mindeep)then
				Q1(i)=mindeep
				do j=1,faceEdgesNum(i)		
					Qb1(i,j)=Q1(i)
				end do
		end if
	end do
    

	!calculate the water depth Q1 at each edge

	call calc_edgeQ1

	do i=1,nFaces	
		
		curFace=i

		if(calc_flag(i)==1.or.calc_flag(i)==2)then
			call border_averages
			call lambda_calc
			call R_calc
			call L_calc
			call modA_calc
			call FI_of_Qb_calc
			call HiVi_calc
       endif
	end do


	do i=1,nFaces
		
		curFace=i

		if(calc_flag(i)==1.or.calc_flag(i)==2)then
			call FI_calc
			call FV_calc        
			call FdotN_calc
			call RofQi_calc
		end if
	end do

	call int_coeffs

	do 115 i=1,nFaces

		curFace=i

		if(calc_flag(i)==1.or.calc_flag(i)==2)then
			call integrate

			if(Q1(i)<=drydeep)then     
				Q2(i)=0
				Q3(i)=0
			end if

			if(Q1(i)<=drydeep)then	
				UM(i)=0
				VN(i)=0
               eta(i)=faceCenters(i,3)
			else
				UM(i)=Q2(i)/Q1(i)
				VN(i)=Q3(i)/Q1(i)
				eta(i)=faceCenters(i,3)+Q1(i)
			end if

	
      else
	     Q2(i)=0
	 	 Q3(i)=0
		 UM(i)=0
		 VN(i)=0
         eta(i)=faceCenters(i,3)
	endif
	      

115	end do

do i=1,nFaces   !calculate tsunami force and maximum force
   
	curFace=i
	

   if(Q1(i).le.drydeep)then 
          Q1(i)=drydeep
          Q1max(i)=Q1(i)
	elseif(Q1(i).gt.drydeep.and.Q1(i).gt.Q1max(i))then
	      Q1max(i)=Q1(i)
   endif
   if(Q1(i).le.drydeep)then
      UM(i)=0.0
      UMmax(i)=UM(i)
   elseif(Q1(i).gt.drydeep.and.UM(i).gt.UMmax(i))then
        UMmax(i)=UM(i)
        VNmax(i)=VN(i)
   endif
   if(Q1(i).le.drydeep)then 
    Tsunamiforce(i)=0.0
	csed(i)=Tsunamiforce(i)
	elseif(Q1(i).gt.drydeep.and.UM(i).gt.0.0)then
	   Tsunamiforce(i)=0.5*1000.0*Q1(i)*(UM(i)**2+VN(i)**2) 
	   if(Tsunamiforce(i).gt.csed(i))then
	   csed(i)=Tsunamiforce(i)
	   endif
   endif
 enddo


	


	!update the values of Q for ghost cells
	call update_ghostcell_value

	!interpolate cell center Q, U and V to nodes
	call cellCenterToNodes(Q1,nodeQ1)
	call cellCenterToNodes(Q2,nodeQ2)
	call cellCenterToNodes(Q3,nodeQ3)
	call cellCenterToNodes(UM,nodeU)
	call cellCenterToNodes(VN,nodeV)
	call cellCenterToNodes(eta,nodeZsurf)
	call cellCenterToNodes(deep,nodedeep)
	call cellCenterToNodes(csed,nodecsed)
    call cellCenterToNodes(Q1max,nodeQ1max)
	call cellCenterToNodes(UMmax,nodeUMmax)
	call cellCenterToNodes(VNmax,nodeVNmax)
	!calculate the Courant number
	call calc_Courant

	calc_flag=0

	return
	end

!c
!c********************************************************************
!c  calculate the inlet flow rate per unit width for each inlet cell   
	subroutine inputcells_q_calc()
	USE COMMON_MODULE,ONLY: inputcellno,inputcelli,inputcellq,edgeMarkers,&
	                        faceEdges,edgeLength,t,faceEdgesNum,faceCenters

	implicit none
	integer i,j
	integer itemp,Zno
	real*8 zp,Ztemp,dl
	real*8 qalpha(1000),qtotal,qtemp(1000)
	real*8 inputq
	integer NoUse
	
!	call get_input_qtotal(t,inputq)

!	qalpha=1.0/inputcellno

		do i=1,inputcellno
		dl=-1
	    !find the length of the inlet side of current cell
		do j=1,faceEdgesNum(i)
			if(edgeMarkers(faceEdges(inputcelli(i),j)).eq.1)then
				dl=edgeLength(faceEdges(inputcelli(i),j))
			end if
		end do
		
		if(dl.lt.0) then
			write(*,*) 'Error on the inlet cells'			
			stop
		endif
        
		if(faceCenters(inputcelli(i),1).gt.406150)then
		    inputcellq(i)=57/3.0/dl

		else
            inputcellq(i)=101/2.0/dl
		endif

	end do

	end subroutine

!c
!c********************************************************************
!c  total inlet flow rate at given time
	subroutine get_input_qtotal(int,inq)
	USE COMMON_MODULE
	implicit none
	integer i
	real*8 inq,int,intt

	call from_t_get_q(intt,inq)

	end subroutine

!c********************************************************************
!c  inquery the total flow rate from the table
	subroutine from_t_get_q(int,inq)
	USE COMMON_MODULE
	integer i
	real*8 inq,int,intt

	!not implemented yet. For now, just uniform distribution
	inq = 101

 403	end subroutine

!c********************************************************************
!c  get inlet q for cell ii
	subroutine get_inputcell_q(ii,temp)
	USE COMMON_MODULE,ONLY: inputcellno,inputcelli,inputcellq
	implicit none
	integer ii,j
	real*8 temp

	do j=1,inputcellno
		if(ii==inputcelli(j))then
			temp=inputcellq(j)
			goto 404
		end if
	end do	
404	end subroutine


!c
!c***************************************************************
!c  This subroutine performs the integration of the governing equations
!c  using the second order Adams-Bashforth time-stepping procedure.

	subroutine integrate
	USE COMMON_MODULE,ONLY: Q1,Q2,Q3,face2DArea,curFace,a,b,d,&
		Rem1,Rem2,Rem3,oldRem1,oldRem2,oldRem3,&
		t,dt,ts,td,nFaces,Sox,Soy,&
		tauwx,tauwy,Swx,Swy,nb,coarse,frctl,g,drydeep,mindeep,&
		sedimentctl,scalc,sedInterval,&
		faceNeighbors,binfo,faceEdgesNum

	implicit none

	integer k,j
	real*8 temp2,temp3,temp1
	real*8 Umax,Uvec,Upar,Vpar,Upar2,Vpar2,Uscale
	real*8 stemp,ddt,wsl,Hmax,Htemp
   
	Q1(curFace)=Q1(curFace)+((dt/face2DArea(curFace))*((a*Rem1(curFace))+(b*oldRem1(curFace))))
	Q2(curFace)=Q2(curFace)+((dt/face2DArea(curFace))*((a*Rem2(curFace))+(b*oldRem2(curFace))))       
	Q3(curFace)=Q3(curFace)+((dt/face2DArea(curFace))*((a*Rem3(curFace))+(b*oldRem3(curFace)))) 

    do j=1,faceEdgesNum(curFace)
	 if(faceNeighbors(curFace,j).gt.0)then
		k=faceNeighbors(curFace,j)
	 else
		k=0
	 endif
		
	 if(k.gt.0.and.Q1(k)>drydeep.and.Q1(curFace)<=0.0)then
	    Q1(k)=Q1(k)+Q1(curFace)*face2DArea(curFace)/face2DArea(k)
	    Q1(curFace)=0.0
     endif
    enddo
    
	Q1(curFace)=max(Q1(curFace),drydeep)
   
    
	return
	end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   calculate the border values of conserved variables
	subroutine border_values
	USE COMMON_MODULE,ONLY: curFace,binfo,faceEdgesNum
	implicit none
	
	integer*4 j

	do j=1,faceEdgesNum(curFace)

		if(binfo(curFace,j)==4)then     

			call bvalue_calc1(curFace,j)

		else
			call bvalue_calc2(curFace,j)
			call gbvalue_calc(curFace,j) 

		end if
	end do


    return
	end


!c
!c***************************************************************
!c
	subroutine bvalue_calc1(i,j)
	USE COMMON_MODULE,ONLY: Qb1,Qb2,Qb3,Q1,Q2,Q3,gQ1,gQ2,gQ3,edgeCenterCoor,&
	                        faceCenters,faceEdges,gradQ1,gradQ2,gradQ3,faceLimiters,&
							eta,gEta,gradEta,Etab,etaLimiter,drydeep,curEdge,edgedrywet
	implicit none

	integer i,j,k
	real*8 phi
	real*8 Q1FN,Q2FN,Q3FN,EtaFN
	real*8 gradTemp(3),rtemp(3)
	real*8 temp,Qedge

    if(Q1(i)>drydeep)then

	call neighbor_value(eta,gEta,i,j,EtaFN)
	call neighbor_value(Q1,gQ1,i,j,Q1FN)
	call neighbor_value(Q2,gQ2,i,j,Q2FN)
	call neighbor_value(Q3,gQ3,i,j,Q3FN)

	!r vector: from face center to edge center
	rtemp(1)=edgeCenterCoor(faceEdges(i,j),1)-faceCenters(i,1)
	rtemp(2)=edgeCenterCoor(faceEdges(i,j),2)-faceCenters(i,2)
	rtemp(3)=0.0

	!gradEta
	gradTemp(1)=gradEta(i,1)
	gradTemp(2)=gradEta(i,2)
	gradTemp(3)=0.0
	call vec_dot(gradTemp,rtemp,temp)
	!unlimited value
	Qedge=Q1(i) + temp
	
	Etab(i,j)=Eta(i)  + etaLimiter(i)*temp
	Qb1(i,j)=Etab(i,j) - edgeCenterCoor(faceEdges(i,j),3)
   

	!gradQ2
	gradTemp(1)=gradQ2(i,1)
	gradTemp(2)=gradQ2(i,2)
	gradTemp(3)=0.0
	call vec_dot(gradTemp,rtemp,temp)
	!unlimited value
	Qedge=Q2(i) + temp

	Qb2(i,j)=Q2(i) + faceLimiters(i,2)*temp


	!gradQ3
	gradTemp(1)=gradQ3(i,1)
	gradTemp(2)=gradQ3(i,2)
	gradTemp(3)=0.0
	call vec_dot(gradTemp,rtemp,temp)
	!unlimited value
	Qedge=Q3(i) + temp

	Qb3(i,j)=Q3(i) + faceLimiters(i,3)*temp

     curEdge = faceEdges(i, j)
			 
    if(edgedrywet(curEdge)==0)then
       Qb1(i,j)=drydeep
       Qb2(i,j)=0.0
	   Qb3(i,j)=0.0
	endif

	else
	 Qb1(i,j)=drydeep
	 Qb2(i,j)=0.0
	 Qb3(i,j)=0.0
	endif

	end subroutine
	
!c***************************************************************
	subroutine bvalue_calc2(i,j)
	USE COMMON_MODULE,ONLY:Qb1,Qb2,Qb3,Q1,Q2,Q3,binfo,mindeep, &
	                       g,inctr,outctr,inletH
	implicit none

	integer i,j,k
	real*8 Q1FN,Q2FN,Q3FN,Q1tmp,Q2tmp,Q3tmp
	real*8 phi1,phi2,phi3
	real*8 uB,uI,hI,hB
	
	!piecewise constant interpolation for boundary direction edge
	Qb1(i,j)=Q1(i) 
	Qb2(i,j)=Q2(i)
	Qb3(i,j)=Q3(i)

    end subroutine

!****************************************************************
!   calculate the boundary value of q for corresponding ghost cells
	subroutine gbvalue_calc(i,j)
	USE COMMON_MODULE,ONLY:Qb1,Qb2,Qb3,Q1,Q2,Q3,binfo,mindeep, &
	                       g,inctr,outctr,gQb1,gQb2,gQb3,faceEdges,&
						   boundaryEdgeGhostCells,faceCenters,csed,gcsed,&
						   inletH,inletUB,inletVB,inletQ,outletH,&
						   outletUB,outletVB,outletQ,gQ1,gQ2,gQ3,faceEdgeNormals,&
	                       aa,nDEMPoints,wse,t
	implicit none
	integer i,j,k,ghostCell
	real*8 Q1FN,Q2FN,Q3FN,Q1tmp,Q2tmp,Q3tmp,csedtmp,temp,temp1
	real*8 phi1,phi2,phi3
	real*8 uB,vB,uI,vI,hI,hB,UQt,h
	real*8 qcondition_inlet
	real*8 nx,ny
	real*8 waterlevel


	
	ghostCell=boundaryEdgeGhostCells(faceEdges(i,j))    

	nx=faceEdgeNormals(i,j,1)
	ny=faceEdgeNormals(i,j,2)

	hI=Qb1(i,j)
	if(hI<=mindeep)then 
		hI=mindeep
	end if

	uI=Qb2(i,j)/hI
	vI=Qb3(i,j)/hI

	UQt=uI*nx+vI*ny

	if(binfo(i,j)==1)then
     
			call get_inputcell_q(i,temp) 
			
			temp1=-1.0*temp       

!			if(Q2tmp>0)then
!				if((Q2tmp/dsqrt(g*hI)/hI).lt.1) then !subcritical					
					hB=qcondition_inlet(UQt,hI,temp1)
!					print*,hb,Qb1(i,j)
!					pause
					Q1tmp=hB
                    Q2tmp=temp*(-nx)
					Q3tmp=temp*(-ny)
!print*,temp,Q2tmp,Q3tmp
!				else                                 !supercritical
!					write(*,*) 'Supercritical inlet'
!					hB=inletH
!					vB=inletVB
!					Q1tmp = hB
!					Q3tmp = vB*hB
!				end if
!			else
!				Q1tmp=Qb1(i,j)
!				Q3tmp=Qb3(i,j) 
!			end if
		

		gQb1(ghostCell)=Q1tmp
		gQb2(ghostCell)=Q2tmp
		gQb3(ghostCell)=Q3tmp


	 elseif(binfo(i,j)==2)then  
       
  	      do aa=1,nDEMPoints-1
            if(t.le.1200.0)then
               waterlevel=3.0*sin(2*3.14/1200*t)
             else
			   waterlevel=0.0
			 endif
		 enddo
         waterlevel=102.0+waterlevel
	     hB= waterlevel 

	     uB=Qb2(i,j)/hI+2*dsqrt(g)*(dsqrt(hI)-dsqrt(hB))
	     vB=Qb3(i,j)/hI+2*dsqrt(g)*(dsqrt(hI)-dsqrt(hB))
!
!         gQb1(ghostCell)= hB  !Qb1(i,j)
!       	 gQb2(ghostCell)= Qb2(i,j)  !gQb1(ghostCell)*uI   !
!		 gQb3(ghostCell)= Qb3(i,j)  !gQb1(ghostCell)*vI   !

	     if(uI.gt.0.0)then
            gQb1(ghostCell)= hB  !Qb1(i,j)
       	    gQb2(ghostCell)= hB*uB  !gQb1(ghostCell)*uI   !
		    gQb3(ghostCell)= 0 !  

	   else
            gQb1(ghostCell)= Qb1(i,j)
       	    gQb2(ghostCell)= Qb2(i,j)  !gQb1(ghostCell)*uI   !
		    gQb3(ghostCell)= 0  !

	     endif
	   

	elseif(binfo(i,j)==3)then   

		gQb1(ghostCell)=Q1(i)     
		gQb2(ghostCell)=0.0D0
		gQb3(ghostCell)=0.0D0

	else
		write(*,*) 'Wrong!'
		stop
	end if

	!update the ghost cells values: there are the same as the border values
	gQ1(ghostCell)=gQb1(ghostCell)
	gQ2(ghostCell)=gQb2(ghostCell)
	gQ3(ghostCell)=gQb3(ghostCell)

	end subroutine
!c***********************************************************************************
!FIofQb1   = Inviscid flux based on Riemann states at cell interface - continuity
!FIofQb2   = Inviscid flux based on Riemann states at cell interface - x-momentum
!FIofQb3   = Inviscid flux based on Riemann states at cell interface - y-momentum
!FI1       = Inviscid flux based on Roe's flux funtion at cell interface - continuity 
!FI2       = Inviscid flux based on Roe's flux funtion at cell interface - x-momentum
!FI3       = Inviscid flux based on Roe's flux funtion at cell interface - y-momentum
	subroutine FI_of_Qb_calc
	USE COMMON_MODULE,ONLY: Qb1,Qb2,Qb3,FI1,FI2,FI3,FIofQb1,FIofQb2,FIofQb3,&
		g,curFace,faceEdgeNormals,Q1,Q2,Q3,drydeep,mindeep,&
		scalc,sedInterval,t,dt,ts,td,sedimentctl,u,v,boundaryEdgeGhostCells,&
		gFIofQb1,gFIofQb2,gFIofQb3,faceNeighbors,faceEdges,gQb1,&
		edgeCenterCoor,faceEdgesNum,gQb2,gQb3
	implicit none
	
	integer j,ghostCell
	real*8 infI1,infI2,infI3
	real*8 gI1,gI2,gI3


	do j=1,faceEdgesNum(curFace)
		if(Qb1(curFace,j)<=drydeep)then	  
			infI1=0.0
			infI2=0.0 
			infI3=0.0
			gI1=0.0 
			gI2=0.0
			gI3=0.0   
		else
			infI1=Qb2(curFace,j)
			infI2=Qb2(curFace,j)**2/Qb1(curFace,j)+0.5*g*Qb1(curFace,j)**2
			infI3=Qb2(curFace,j)*Qb3(curFace,j)/Qb1(curFace,j)
			gI1=Qb3(curFace,j)
			gI2=infI3
			gI3=Qb3(curFace,j)**2/Qb1(curFace,j)+0.5*g*Qb1(curFace,j)**2
	    end if
	
		FIofQb1(curFace,j)=(infI1*(faceEdgeNormals(curFace,j,1)))+(gI1*(faceEdgeNormals(curFace,j,2)))
		FIofQb2(curFace,j)=(infI2*(faceEdgeNormals(curFace,j,1)))+(gI2*(faceEdgeNormals(curFace,j,2)))
		FIofQb3(curFace,j)=(infI3*(faceEdgeNormals(curFace,j,1)))+(gI3*(faceEdgeNormals(curFace,j,2)))

		!for ghost cells
		if(faceNeighbors(curFace,j).lt.0) then
			ghostCell=boundaryEdgeGhostCells(faceEdges(curFace,j))


			if(gQb1(ghostCell)<=drydeep)then	  
				infI1=0.0
				infI2=0.0 
				infI3=0.0
				gI1=0.0 
				gI2=0.0
				gI3=0.0   
			else
			infI1=gQb2(ghostCell)
			infI2=gQb2(ghostCell)**2/gQb1(ghostCell)+0.5*g*gQb1(ghostCell)**2
			infI3=gQb2(ghostCell)*gQb3(ghostCell)/gQb1(ghostCell)
			gI1=gQb3(ghostCell)
			gI2=infI3
			gI3=gQb3(ghostCell)**2/gQb1(ghostCell)+0.5*g*gQb1(ghostCell)**2
			end if

			!here the minus sign in front of the edge normal vector is because the outside normal
			!vector for the ghost cell is the opposite of the inside cell
			gFIofQb1(ghostCell)=(infI1*(-faceEdgeNormals(curFace,j,1)))+(gI1*(-faceEdgeNormals(curFace,j,2)))
			gFIofQb2(ghostCell)=(infI2*(-faceEdgeNormals(curFace,j,1)))+(gI2*(-faceEdgeNormals(curFace,j,2)))
			gFIofQb3(ghostCell)=(infI3*(-faceEdgeNormals(curFace,j,1)))+(gI3*(-faceEdgeNormals(curFace,j,2)))

		end if
	end do

	return
	end

!c**************************************************************************
!c      This subroutine calculates the viscous fluxes at each cell interface      
!c      with the velocity gradients being effectively central difference,
!c      based on the cell centre values of consecutive cells.
	subroutine FV_calc
	USE COMMON_MODULE,ONLY: g,edgeGradU,edgeGradV,curFace,faceEdgeNormals,&
							visc,drydeep,mindeep,QB1,nb,u,v,&
							edgeCenterCoor,faceEdges,Qb1av,FV2,FV3,faceEdgesNum
	implicit none

	integer j,edgeNum
	real*8 ht
	real*8 nx,ny
	
	do j=1,faceEdgesNum(curFace)
      if(Qb1(curFace,j)>drydeep)then
	  
	    visc=0.2*0.4*Qb1(curFace,j)*dsqrt((g*nb(curFace)**2*(u(curFace,j)**2+v(curFace,j)**2))/Qb1(curFace,j)**(1/3))	
		
    	  ht=Qb1av(curFace,j)
		
		  nx=faceEdgeNormals(curFace,j,1)
		  ny=faceEdgeNormals(curFace,j,2)
	
		  edgeNum=faceEdges(curFace,j)

		  FV2(curFace,j)=-visc*(ht*edgeGradU(edgeNum,1)*nx+ht*edgeGradU(edgeNum,2)*ny)
		  FV3(curFace,j)=-visc*(ht*edgeGradV(edgeNum,1)*nx+ht*edgeGradV(edgeNum,2)*ny)
	  end if
	end do

	end subroutine

!c***************************************************************
!c  This subroutine calculates the results of the vector product F.n
	subroutine FdotN_calc
	USE COMMON_MODULE,ONLY: FdotN1,FdotN2,FdotN3,curFace,FI1,FI2,FI3,&
		FIofQb1,FIofQb2,FIofQb3,FV1,FV2,FV3,sedimentctl,t,dt,ts,td,&
		scalc,sedInterval,binfo,edgeLength,faceEdges,ELEDGES,faceEdgesNum
	implicit none

	real*8 length(ELEDGES)
	integer i
	
	do i=1,faceEdgesNum(curFace)
		length(i)=edgeLength(faceEdges(curFace,i))
	end do

	FdotN1(curFace)=0.0
	FdotN2(curFace)=0.0
	FdotN3(curFace)=0.0
	do i=1,ELEDGES
		FdotN1(curFace)=FdotN1(curFace)+FI1(curFace,i)*length(i)
		FdotN2(curFace)=FdotN2(curFace)+(FI2(curFace,i)+FV2(curFace,i))*length(i)    
		FdotN3(curFace)=FdotN3(curFace)+(FI3(curFace,i)+FV3(curFace,i))*length(i)   

	end do

	return       
	end

!c***************************************************************
!c  This subroutine calculates the source term vector components
	subroutine HiVi_calc
	USE COMMON_MODULE,ONLY: maxfaces_,g,Q1,Q2,Q3,face2DArea,HiVi1,HiVi2,HiVi3,curFace,&
		Sox,Soy,tauwx,tauwy,Swx,Swy,t,dt,ts,td,nb,frctl,drydeep,mindeep,sedimentctl,&
		scalc,sedInterval, edgeCenterCoor, faceEdges, faceEdgeNormals,edgeLength,faceEdgesNum,&
		gradEta,edgeQ1,Cd,faceCenters
	implicit none

	integer i
    real*8 Sf,Sfx,Sfy,ssfx,ssfy,vfx,vfy
	real*8 Chezy(maxfaces_),ARR

	real*8 length(10),nx(10),ny(10),h(10)
	real*8 slopeTermx,slopeTermy

!   bottom friction terms
	ssfx=nb(curFace)**2*g*Q2(curFace)*dsqrt(Q2(curFace)**2+Q3(curFace)**2)/&
     		(Q1(curFace)**(7.0/3.0))*dsqrt(1+Sox(curFace)**2+Soy(curFace)**2)
	ssfy=nb(curFace)**2*g*Q3(curFace)*dsqrt(Q2(curFace)**2+Q3(curFace)**2)/&
     		(Q1(curFace)**(7.0/3.0))*dsqrt(1+Sox(curFace)**2+Soy(curFace)**2)

!   vegetation resistance
    vfx=0.0
	vfy=0.0
	if(Q1(curFace).lt.0.5114754)then
	  Cd(curFace)=8.994709
    elseif(Q1(curFace).ge.0.5114754.and.Q1(curFace).lt.0.9836066)then
	  Cd(curFace)=8.994709+(Q1(curFace)-0.5114754)/(0.9836066-0.5114754)*(7.107584-8.994709)
	
	elseif(Q1(curFace).ge.0.9836066.and.Q1(curFace).lt.1.495082)then
	 Cd(curFace)=7.107584+(Q1(curFace)-0.9836066)/(1.495082-0.9836066)*(5.908289-7.107584)
	
	 elseif(Q1(curFace).ge.1.495082.and.Q1(curFace).lt.2.045902)then
	 Cd(curFace)=5.943563+(Q1(curFace)-1.495082)/(2.045902-1.495082)*(5.026455-5.908289)

	elseif(Q1(curFace).ge.2.045902.and.Q1(curFace).lt.2.518033)then
	Cd(curFace)=5.026455+(Q1(curFace)-2.045902)/(2.518033-2.045902)*(4.268078-5.026455)

	elseif(Q1(curFace).ge.2.518033.and.Q1(curFace).lt.3.108197)then
	Cd(curFace)=4.268078+(Q1(curFace)-2.518033)/(3.108197-2.518033)*(3.703704-4.268078)
	
	elseif(Q1(curFace).ge.3.108197.and.Q1(curFace).lt.3.619672)then
	Cd(curFace)=3.703704+(Q1(curFace)-3.108197)/(3.619672-3.108197)*(3.368607-3.703704)

	elseif(Q1(curFace).ge.3.619672.and.Q1(curFace).lt.4.131148)then
	Cd(curFace)=3.368607+(Q1(curFace)-3.619672)/(4.131148-3.619672)*(3.24515-3.368607)

    elseif(Q1(curFace).ge.4.131148.and.Q1(curFace).lt.4.52459)then
	Cd(curFace)=3.24515+(Q1(curFace)-4.131148)/(4.52459-4.131148)*(3.227513-3.24515)

	elseif(Q1(curFace).ge.4.52459.and.Q1(curFace).lt.5.07541)then
	Cd(curFace)=3.227513+(Q1(curFace)-4.52459)/(5.07541-4.52459)*(3.174603-3.227513)

	
	elseif(Q1(curFace).ge.5.07541.and.Q1(curFace).lt.5.62623)then
	Cd(curFace)=3.174603+(Q1(curFace)-5.07541)/(5.62623-5.07541)*(3.104056-3.174603)

	elseif(Q1(curFace).ge.5.62623.and.Q1(curFace).lt.6.098361)then
	Cd(curFace)=3.104056+(Q1(curFace)-5.62623)/(6.098361-5.62623)*(3.015873-3.104056)

    elseif(Q1(curFace).ge.6.098361.and.Q1(curFace).lt.6.531148)then
	Cd(curFace)=3.015873+(Q1(curFace)-6.098361)/(6.531148-6.098361)*(2.945326-3.015873)

   elseif(Q1(curFace).ge.6.531148.and.Q1(curFace).lt.7.081967)then
	Cd(curFace)=2.945326+(Q1(curFace)-6.531148)/(7.081967-6.531148)*(2.857143-2.945326)

    elseif(Q1(curFace).ge.7.081967.and.Q1(curFace).lt.7.632787)then
	Cd(curFace)=2.857143+(Q1(curFace)-7.081967)/(7.632787-7.081967)*(2.680776-2.857143)

    elseif(Q1(curFace).ge.7.632787.and.Q1(curFace).lt.8.02623)then
	Cd(curFace)=2.680776+(Q1(curFace)-7.632787)/(8.02623-7.632787)*(2.557319-2.680776)
	else

	Cd(curFace)=2.557319

	endif
	 
if(faceCenters(curFace,1).ge.5700.0.and.faceCenters(curFace,1).le.5800.0)then
       vfx=0.5*0.1*Cd(curFace)*0.224*min(8.0,Q1(curFace))*Q2(curFace)*dsqrt(Q2(curFace)**2+Q3(curFace)**2)/(Q1(curFace)**2)
       vfy=0.5*0.1*Cd(curFace)*0.224*min(8.0,Q1(curFace))*Q3(curFace)*dsqrt(Q2(curFace)**2+Q3(curFace)**2)/(Q1(curFace)**2)
endif




!   slope terms
	do i=1,faceEdgesNum(curFace)
		length(i)=edgeLength(faceEdges(curFace,i))
		h(i)=edgeQ1(faceEdges(curFace,i))
		nx(i)=faceEdgeNormals(curFace,i,1)
		ny(i)=faceEdgeNormals(curFace,i,2)
	enddo
	
	slopeTermx=0.0
	slopeTermy=0.0
	
	do i=1,faceEdgesNum(curFace)
		slopeTermx=slopeTermx + 0.5*g*h(i)**2*nx(i)*length(i)
		slopeTermy=slopeTermy + 0.5*g*h(i)**2*ny(i)*length(i)
	enddo

	slopeTermx=slopeTermx - g*Q1(curFace)*gradEta(curFace,1)*face2DArea(curFace)
	slopeTermy=slopeTermy - g*Q1(curFace)*gradEta(curFace,2)*face2DArea(curFace)
    

!        ***  new full expression: slope source term is well balanced      **

	HiVi2(curFace)=(Swx-ssfx-vfx)*face2DArea(curFace) + slopeTermx
	HiVi3(curFace)=(Swy-ssfy-vfy)*face2DArea(curFace) + slopeTermy

!c       ***  old full expression: slope source term is not well balanced  **
!	HiVi2(curFace)=(Swx-ssfx+(g*Q1(curFace)*Sox(curFace)))*face2DArea(curFace)
!	HiVi3(curFace)=(Swy-ssfy+(g*Q1(curFace)*Soy(curFace)))*face2DArea(curFace)
	
	return
	end

!c
!c***************************************************************
!c  This subroutine calculates the summation of the F.n term 
!c  integrated round the cell and the source term vector
!c  The subroutine also passes the previous time step value to an
!c  old array for use in the AB2 integration.
	subroutine RofQi_calc
	USE COMMON_MODULE,ONLY: maxfaces_,Rem1,Rem2,Rem3,oldRem1,oldRem2,oldRem3,&
		FdotN1,FdotN2,FdotN3,HiVi1,HiVi2,HiVi3,curFace,t,dt,ts,td,drydeep,mindeep&
		,Q1,Q2,Q3,sedimentctl,scalc,sedInterval,calc_flag,faceCenters
	implicit none

    real*8 Q1FN,etaFN
	integer j

	if(t.gt.0D0)then
	  oldRem1(curFace)=Rem1(curFace)
	  oldRem2(curFace)=Rem2(curFace)
	  oldRem3(curFace)=Rem3(curFace)
	end if

	Rem1(curFace)=-FdotN1(curFace)
	Rem2(curFace)=-FdotN2(curFace)+HiVi2(curFace)
	Rem3(curFace)=-FdotN3(curFace)+HiVi3(curFace)

	
    if(Q1(curFace)<=drydeep)then  
		Rem2(curFace)=0.0
		Rem3(curFace)=0.0

		oldRem2(curFace)=0.0
		oldRem3(curFace)=0.0
	end if


	return
	end

!c***************************************************************
!C  This subroutine calculates the limiter value
!C  Q0: value at cell center
!C  Qedge: unlimited edge value
!C  Qnb:value at neighbor cell center
!C  beta: limiter parameter. (constant value if beta < 0)
	subroutine limiter(Q0,Qedge,Qnb,phi)
	USE COMMON_MODULE,ONLY: beta,VSMALL
	implicit none
	
	real*8 Q0,Qedge,Qnb,phi
	real*8 Q0_min,Q0_max
	real*8 r

	Q0_min=min(Q0,Qnb)
	Q0_max=max(Q0,Qnb)

	if(beta<0) then
		phi=0.0
	else
		if(dabs(Qedge-Q0).lt.VSMALL)then
			r=1.0
		elseif(Qedge>Q0)then
			r=(Q0_max-Q0)/(Qedge-Q0)
		else
			r=(Q0_min-Q0)/(Qedge-Q0)
		endif

!		phi=max(min(beta*r,1.0),min(r,beta))
		phi=min(r,1.0)
	endif

	end subroutine

!c********************************************************************
!c  This subroutine returns the value of j direction neighbor of cell cc
	subroutine neighbor_value(h,gh,cc,j,hFN)
	USE COMMON_MODULE,ONLY: maxfaces_,maxboundaryedges_,nFaces,binfo,&
				edgeMarkers,faceEdges,faceNeighbors,boundaryEdgeGhostCells
	implicit none
	
	integer cc,j,ghostCell
	real*8 h(maxfaces_),gh(maxboundaryedges_),hFN		

	if(edgeMarkers(faceEdges(cc,j))/=4)then	
		ghostCell=boundaryEdgeGhostCells(faceEdges(cc,j))

!		print*,ghostCell,'ghostCell---1'

		hFN=gh(ghostCell)
	else                                    
		hFN=h(faceNeighbors(cc,j))
    end if

    end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!   find the position of edge i in face's edges list
    function edgePositionInFace(face, edge)
	USE COMMON_MODULE,ONLY: faceEdges,faceEdgesNum
    implicit none

	integer i,face,edge    
	integer edgePositionInFace
 
	do i=1,faceEdgesNum(face)
		edgePositionInFace = -1
		if(faceEdges(face,i)==edge) then
			edgePositionInFace = i
			exit
		endif	
    enddo
	          
    return     
    end

!c********************************************************************
!c  This subroutine returns the neighbour values  of
!c  Riemann state values
	subroutine neighbor_value2(h,gh,cc,i,hFN)
	USE COMMON_MODULE,ONLY:maxfaces_,maxboundaryedges_,nFaces,binfo,&
			edgeMarkers,faceEdges,faceNeighbors,boundaryEdgeGhostCells,ELEDGES
	implicit none
	
	integer cc,i,pos,ghostCell
	real*8 h(maxfaces_,ELEDGES),gh(maxboundaryedges_),hFN

	integer edgePositionInFace

	if(edgeMarkers(faceEdges(cc,i))/=4)then	
		ghostCell=boundaryEdgeGhostCells(faceEdges(cc,i))

!		print*,ghostCell,'ghostCell'
!		pause

		hFN=gh(ghostCell)
	else
		pos=edgePositionInFace(faceNeighbors(cc,i),faceEdges(cc,i))
		if(pos==-1)then
			write(*,*) 'Something wrong!'
			stop
		else
			hFN=h(faceNeighbors(cc,i),pos)
		end if
    end if

    end subroutine

!c***************************************************************
!   solve characteristic equation for the inlet   	
!   given unit discharge-->water depth
!   using Newton-Raphson method
	real*8 function qcondition_inlet(uI,hI,Q2tmp)
	USE COMMON_MODULE,ONLY: g,drydeep,mindeep,VLARGE
	implicit none
	REAL*8 uI,hI,Q2tmp
	REAL*8 xk,xkplus1,fxk,fprimexk,err
	
	if(hI<=mindeep)then 	
		hI=mindeep
	end if

	err=VLARGE
	xk=hI
	do while(err>1E-6)  
		fxk=xk-(1/2.0/dsqrt(g)*(uI-Q2tmp/xk)+dsqrt(hI))**2
		fprimexk=1-Q2tmp*(dsqrt(hI)+(-Q2tmp/xk+uI)/(2*dsqrt(g)))/dsqrt(g)/xk**2
		xkplus1=xk-fxk/fprimexk
		err = dabs(xkplus1-xk)
		xk=xkplus1
	end do
	
	qcondition_inlet=xkplus1
	
	end function qcondition_inlet


!c***************************************************************
!c  This subroutine calculates the Roe-type Riemann averages at
!c  each cell interface for h, u and v
	subroutine border_averages
	USE COMMON_MODULE,ONLY:nFaces,Qb1,Qb2,Qb3,Q1,Q2,Q3,Qb1av,&
		Q2in,u,v,alpha,g,curFace,binfo,drydeep,&
		mindeep,t,dt,ts,td,inctr,outctr,inletH,gQb1,gQb2,gQb3,&
		edgeCenterCoor,faceEdges,faceEdgesNum

	implicit none
	integer:: j,status
	real(8):: Q1tmp,Q2tmp,Q3tmp,hs_av,uB,hB,hI,bottom,Qb1FN,QB2FN,Qb3FN
	real(8):: h1,h2
	
	do j=1,faceEdgesNum(curFace)

		if(binfo(curFace,j)==3)then			
		    Qb1av(curFace,j)=Qb1(curFace,j)
 		    u(curFace,j)=0
			v(curFace,j)=0
		else
  		    call neighbor_value2(Qb1,gQb1,curFace,j,Qb1FN)
		    call neighbor_value2(Qb2,gQb2,curFace,j,Qb2FN)
		    call neighbor_value2(Qb3,gQb3,curFace,j,Qb3FN)
	
  		    if(Qb1(curFace,j)<=mindeep)then
		        h1=0+mindeep
	        else
			    h1=dsqrt(Qb1(curFace,j))
	        end if
	
	        if(Qb1FN<=mindeep)then
			    h2=0+mindeep
	        else
			    h2=dsqrt(Qb1FN)
		    end if
	
		    bottom=h1+h2
		    Qb1av(curFace,j)=0.5D0*(Qb1(curFace,j)+Qb1FN)   
									                          
		    if(Qb1(curFace,j)<=mindeep.and.Qb1FN<=mindeep)then
			    u(curFace,j)=0
			    v(curFace,j)=0
		    else if(Qb1(curFace,j)<=mindeep.and.Qb1FN>mindeep)then		 
			    u(curFace,j)=Qb2FN/Qb1FN
			    v(curFace,j)=Qb3FN/Qb1FN
		    else if(Qb1FN<=mindeep.and.Qb1(curFace,j)>mindeep)then		
			    u(curFace,j)=Qb2(curFace,j)/Qb1(curFace,j)
			    v(curFace,j)=Qb3(curFace,j)/Qb1(curFace,j)
		    else
		        u(curFace,j)=(Qb2(curFace,j)/Qb1(curFace,j)*h1+Qb2FN/Qb1FN*h2)/bottom
			    v(curFace,j)=(Qb3(curFace,j)/Qb1(curFace,j)*h1+Qb3FN/Qb1FN*h2)/bottom
		    end if
		end if
	end do

	return
	end


!c***********************************************************************
!c  This subroutine calculates the eigenvalues at each cell interface used 
!c  in the eigenvalue/eigenvector decomposition of the flux jacobian A
	subroutine lambda_calc
	USE COMMON_MODULE,ONLY: Qb1av,Qb1,Qb2,Qb3,g,Q1,Q2,Q3,&
		lambda1,lambda2,lambda3,curFace,modA,R,modLAMBDA,L,faceEdgeNormals,&
		u,v,alpha,drydeep,mindeep,edgeCenterCoor,faceEdges,faceEdgesNum
	implicit none
	
	integer j
	
	do j=1,faceEdgesNum(curFace)
		if(Qb1av(curFace,j)<=mindeep)then	  
		    alpha(curFace,j)=0
	    else
			alpha(curFace,j)=dsqrt(g*Qb1av(curFace,j))        
	    end if
		lambda1(curFace,j)=u(curFace,j)*(faceEdgeNormals(curFace,j,1))+&
				           v(curFace,j)*(faceEdgeNormals(curFace,j,2))
		lambda2(curFace,j)=lambda1(curFace,j) - alpha(curFace,j)
		lambda3(curFace,j)=lambda1(curFace,j) + alpha(curFace,j)
	end do

	return
	end
!c***********************************************************************

	subroutine lambda_change
	USE COMMON_MODULE,ONLY: Qb1av,Qb1,Qb2,Qb3,g,Q1,Q2,Q3,&
		lambda1,lambda2,lambda3,curFace,modA,R,modLAMBDA,L,faceEdgeNormals,&
		u,v,alpha,drydeep,mindeep,edgeCenterCoor,faceEdges,faceEdgesNum,dt,&
        gQb1,gQb2,gQb3,edgeLength,pcoor,edgePoints
	implicit none
	
	integer j
	real*8 tolvalue
    real*8 tempL,tempR
    real*8 lambdaL2,lambdaR2,lambdaL3,lambdaR3
    real(8):: Qb1FN,QB2FN,Qb3FN

	do j=1,faceEdgesNum(curFace)
	tolvalue=0.1*dabs(pcoor(edgePoints(faceEdges(curFace, j),1), 1)-pcoor(edgePoints(faceEdges(curFace, j),2), 1))/dt
  	
	call neighbor_value2(Qb1,gQb1,curFace,j,Qb1FN)
    call neighbor_value2(Qb2,gQb2,curFace,j,Qb2FN)
	call neighbor_value2(Qb3,gQb3,curFace,j,Qb3FN)
    
	tempL=(Qb2(curFace,j)/Qb1(curFace,j))*faceEdgeNormals(curFace,j,1)+(Qb3(curFace,j)/Qb1(curFace,j))*faceEdgeNormals(curFace,j,2)    
	tempR=(Qb2FN/Qb1FN)*faceEdgeNormals(curFace,j,1)+(Qb3FN/Qb1FN)*faceEdgeNormals(curFace,j,2)

    lambdaL2=tempL-Qb1(curFace,j)
    lambdaR2=0.25*(tempL+3*tempR)+0.5*(Qb1(curFace,j)-3.0*Qb1FN)

    lambdaL3=0.25*(tempL+3*tempR)+0.5*(Qb1(curFace,j)-3.0*Qb1FN)
    lambdaR3=tempR+Qb1FN
	
	if(dabs(lambda2(curFace,j))<=tolvalue)then
		lambda2(curFace,j)=lambda2(curFace,j)*((lambdaR2-lambda2(curFace,j))/(lambdaR2-lambdaL2))
	endif
	if(dabs(lambda3(curFace,j))<=tolvalue)then
	    lambda3(curFace,j)=lambda3(curFace,j)*((lambdaR3-lambda3(curFace,j))/(lambdaR3-lambdaL3))
	endif

	end do

	return
	end
!c***********************************************************************

!c  The subroutine calculates the right eigenvector matrix of the 
!c  flux jacobian A
	subroutine R_calc
	USE COMMON_MODULE,ONLY: lambda1,lambda2,lambda3,curFace,modA,R,modLAMBDA,L,&
		faceEdgeNormals,u,v,alpha,faceEdgesNum

	implicit none
	integer j

	do j=1,faceEdgesNum(curFace)
		R(curFace,j,1,1)=0D0
		R(curFace,j,1,2)=1D0
		R(curFace,j,1,3)=1D0
		R(curFace,j,2,1)=(faceEdgeNormals(curFace,j,2))
		R(curFace,j,2,2)=u(curFace,j)-alpha(curFace,j)*faceEdgeNormals(curFace,j,1)
		R(curFace,j,2,3)=u(curFace,j)+alpha(curFace,j)*faceEdgeNormals(curFace,j,1)
		R(curFace,j,3,1)=-faceEdgeNormals(curFace,j,1)
		R(curFace,j,3,2)=v(curFace,j)-alpha(curFace,j)*faceEdgeNormals(curFace,j,2)
		R(curFace,j,3,3)=v(curFace,j)+alpha(curFace,j)*faceEdgeNormals(curFace,j,2)
	end do

	return
	end
		
!c
!c***************************************************************
!c  The subroutine calculates the left eigenvector matrix of the 
!c  flux jacobian A

	subroutine L_calc
	USE COMMON_MODULE,ONLY: lambda1,lambda2,lambda3,curFace,modA,R,modLAMBDA,L,&
		faceEdgeNormals,u,v,alpha,faceEdgesNum

	implicit none
	integer j
	real*8 a,b,e,d

	do j=1,faceEdgesNum(curFace)
		a=2.0D0*alpha(curFace,j)
		if(a/=0)then
			b=lambda1(curFace,j)/a
			e=(faceEdgeNormals(curFace,j,1))/a
			d=(faceEdgeNormals(curFace,j,2))/a
		else 
			b=0D0
			e=0D0
			d=0D0
		end if

		L(curFace,j,1,1)=-((u(curFace,j)*(faceEdgeNormals(curFace,j,2)))-&
		                   (v(curFace,j)*(faceEdgeNormals(curFace,j,1))))
		L(curFace,j,1,2)= faceEdgeNormals(curFace,j,2)
		L(curFace,j,1,3)=-faceEdgeNormals(curFace,j,1)
		L(curFace,j,2,1)=0.5D0 + b 
		L(curFace,j,2,2)=-e
		L(curFace,j,2,3)=-d
		L(curFace,j,3,1)=0.5D0 - b
		L(curFace,j,3,2)=e
		L(curFace,j,3,3)=d

	end do

	return
	end

!c***************************************************************
!c  This subroutine calculates the absolute flux jacobian based on 
!c  absolute eigenvalues.  So modA=R*modLAMBDA*L
!c***************************************************************
	subroutine modA_calc
	USE COMMON_MODULE,ONLY: lambda1,lambda2,lambda3,curFace,modA,R,modLAMBDA,L,&
							u,v,alpha,faceEdgesNum
	implicit none
	integer a,j,k
	real*8 matrix1(3,3),matrix2(3,3),matrix3(3,3)
	real*8 matrix4(3,3),matrix5(3,3)

	do j=1,faceEdgesNum(curFace)
		modLAMBDA(curFace,j,1,1)=dabs(lambda1(curFace,j))
		modLAMBDA(curFace,j,2,2)=dabs(lambda2(curFace,j))
		modLAMBDA(curFace,j,3,3)=dabs(lambda3(curFace,j))

		do k=1,3
		  do a=1,3
		    matrix1(k,a)=R(curFace,j,k,a)
		    matrix2(k,a)=L(curFace,j,k,a)
		    matrix3(k,a)=modLAMBDA(curFace,j,k,a)
		  end do
		end do

		matrix4=matmul(matrix3,matrix2)
		matrix5=matmul(matrix1,matrix4)

		do k=1,3
		  do a=1,3
		     modA(curFace,j,k,a)=matrix5(k,a)
		  end do
		end do
		
	end do

	return
	end

!c********************************************************
!c  This subroutine calculates the intercell flux based on 
!c  Roe's flux funtion
	subroutine FI_calc
	USE COMMON_MODULE,ONLY:Qb1,Qb2,Qb3,nFaces,g,&
		curFace,Q2in,FI1,FI2,FI3,FIofQb1,FIofQb2,FIofQb3,&
		modA,R,modLAMBDA,L,faceEdgeNormals,t,dt,ts,td,&
		Q1,Q2,Q3,binfo,drydeep,mindeep,inctr,outctr,inletH,&
		gFIofQb1,gFIofQb2,gFIofQb3,gQb1,gQb2,gQb3,faceEdgesNum,&
		eta,gEta,faceCenters,gQ1,calc_flag,edgePoints,faceEdges,&
		curEdge,edgedrywet,faceNeighbors
		
	implicit none
	integer a,j,k,pos,edgePositionInFace
	real*8 infI1(3),infI2(3),infI3(3)
	real*8 gI1(3),gI2(3),gI3(3)
	real*8 Qb1FN,Qb2FN,Qb3FN
	real*8 FIofQb1FN,FIofQb2FN,FIofQb3FN
    real*8 infI1FN,infI2FN,infI3FN
	real*8 gI1FN,gI2FN,gI3FN
	real*8 dQ(3),result(3),tempmodA(3,3)
	real*8 hI,hB,uB,uI
	real*8 qcondition_inlet
    real*8 Q1FN,Q2FN,Q3FN,EtaFN


	do j=1,faceEdgesNum(curFace)

		!get the invicid flux of neighbor
	    call neighbor_value2(FIofQb1,gFIofQb1,curFace,j,FIofQb1FN)
		call neighbor_value2(FIofQb2,gFIofQb2,curFace,j,FIofQb2FN)
		call neighbor_value2(FIofQb3,gFIofQb3,curFace,j,FIofQb3FN)

		!get the conservative variables of neighbor
		call neighbor_value2(Qb1,gQb1,curFace,j,Qb1FN)	
		call neighbor_value2(Qb2,gQb2,curFace,j,Qb2FN)
		call neighbor_value2(Qb3,gQb3,curFace,j,Qb3FN)
		
		dQ(1)=Qb1FN-Qb1(curFace,j)     
		dQ(2)=Qb2FN-Qb2(curFace,j)
		dQ(3)=Qb3FN-Qb3(curFace,j)
		
		do k=1,3
			do a=1,3
				tempmodA(k,a)=modA(curFace,j,k,a)  
			end do
		end do

		result=matmul(tempmodA,dQ)

		!the signs infront of FIofQb1FN, FIofQb1FN, FIofQb1FN are negative
		!since the normal vector is the negative for the current edge j with
		!respect to each neighboring face
		FI1(curFace,j)=0.5D0*(-FIofQb1FN+FIofQb1(curFace,j)-result(1))
		FI2(curFace,j)=0.5D0*(-FIofQb2FN+FIofQb2(curFace,j)-result(2))
		FI3(curFace,j)=0.5D0*(-FIofQb3FN+FIofQb3(curFace,j)-result(3))
        


        curEdge = faceEdges(curFace, j)
		pos=edgePositionInFace(curFace,curEdge)

		if(pos==-1)then
			write(*,*) 'Something wrong!'
			stop
		else
            if(Q1(curFace)<=drydeep)then
               if(edgedrywet(curEdge)==1)then  			   
		          FI1(curFace,pos)=-FIofQb1FN
		          FI2(curFace,pos)=0.0   
		          FI3(curFace,pos)=0.0
				else
                  FI1(curFace,pos)=0.0D0
		          FI2(curFace,pos)=0.0D0
		          FI3(curFace,pos)=0.0D0      				
				endif   
              
			else

			   if(edgedrywet(curEdge)==0)then
                  FI1(curFace,pos)=0.0D0
		          FI2(curFace,pos)=0.0D0
		          FI3(curFace,pos)=0.0D0  
			   elseif(edgedrywet(curEdge)==1)then
			      FI1(curFace,pos)=FI1(curFace,pos)
		          FI2(curFace,pos)=0.0D0
		          FI3(curFace,pos)=0.0D0 
			   else
                  FI1(curFace,pos)=FI1(curFace,pos)
		          FI2(curFace,pos)=FI2(curFace,pos)
		          FI3(curFace,pos)=FI3(curFace,pos)            
			   endif	
		   end if
       endif
                      		     
      end do
      
      return
      end

!	Interpolate cell center values to node values
	subroutine cellCenterToNodes(cellValues,nodeValues)
	USE COMMON_MODULE,ONLY: faceCenters,nNodes,pointFaces,maxfaces_,maxnodes_,&
							pointNFaces,pcoor,maxnodefaces_,facePoints,VSMALL,&
							Q1,drydeep,calc_flag
	implicit none

	real(8) cellValues(maxfaces_), nodeValues(maxnodes_)
	real*8 dtemp1,dtemp2,temp
	real*8 v1(3),v2(3),c1(3)
	integer i,j

    do i = 1, nNodes                  
         !temp variables for inverse distance weigh
         dtemp1 = 0.0
         dtemp2 = 0.0

         !current node coordinates
         v1(1) = pcoor(i,1)
         v1(2) = pcoor(i,2)
         v1(3) = pcoor(i,3)
         
         do j =1, pointNFaces(i)         
            if(pointFaces(i, j).ne.-1) then
!			  if(calc_flag(pointFaces(i,j))==1.or.calc_flag(pointFaces(i,j))==2)then
                !current face center coordinates
                c1(1) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),1)+&
                                pcoor(facePoints(pointFaces(i,j),2),1)+ &
                                pcoor(facePoints(pointFaces(i,j),3),1))
                c1(2) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),2)+ &
                                pcoor(facePoints(pointFaces(i,j),2),2)+  &
                                pcoor(facePoints(pointFaces(i,j),3),2))
                c1(3) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),3)+ &
                                pcoor(facePoints(pointFaces(i,j),2),3)+  &
                                pcoor(facePoints(pointFaces(i,j),3),3))

                !distance vector from current node to face center
                v2(1) = v1(1) - c1(1)  
                v2(2) = v1(2) - c1(2)  
                v2(3) = v1(3) - c1(3)  
                
                !distance between current node and face center
                call vec_mag(v2, temp)
                dtemp1 = dtemp1 + cellValues(pointFaces(i,j))/temp
                dtemp2 = dtemp2 + 1.0/temp               
               end if
!		   endif
         end do 
         
		 if(dabs(dtemp2).gt.VSMALL)then
			nodeValues(i) = dtemp1/dtemp2
		 else
		    nodeValues(i) = 0.0
		 end if
    end do

	end subroutine	

!	Interpolate cell center values to node values
	subroutine cellCenterToNodes2(cellValues,nodeValues)
	USE COMMON_MODULE,ONLY: faceCenters,nNodes,pointFaces,maxfaces_,maxnodes_,&
							pointNFaces,pcoor,maxnodefaces_,facePoints,VSMALL,&
							Q1,drydeep,calc_flag
	implicit none

	real(8) cellValues(maxfaces_), nodeValues(maxnodes_)
	real*8 dtemp1,dtemp2,temp
	real*8 v1(3),v2(3),c1(3)
	integer i,j

    do i = 1, nNodes                  
         !temp variables for inverse distance weigh
         dtemp1 = 0.0
         dtemp2 = 0.0

         !current node coordinates
         v1(1) = pcoor(i,1)
         v1(2) = pcoor(i,2)
         v1(3) = pcoor(i,3)
         
         do j =1, pointNFaces(i)
        if(pointFaces(i, j).ne.-1) then

		   if(calc_flag(pointFaces(i,j))==1)then                                
                !current face center coordinates
                c1(1) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),1)+&
                                pcoor(facePoints(pointFaces(i,j),2),1)+ &
                                pcoor(facePoints(pointFaces(i,j),3),1))
                c1(2) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),2)+ &
                                pcoor(facePoints(pointFaces(i,j),2),2)+  &
                                pcoor(facePoints(pointFaces(i,j),3),2))
                c1(3) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),3)+ &
                                pcoor(facePoints(pointFaces(i,j),2),3)+  &
                                pcoor(facePoints(pointFaces(i,j),3),3))

                !distance vector from current node to face center
                v2(1) = v1(1) - c1(1)  
                v2(2) = v1(2) - c1(2)  
                v2(3) = v1(3) - c1(3)  
                
                !distance between current node and face center
                call vec_mag(v2, temp)
                dtemp1 = dtemp1 + cellValues(pointFaces(i,j))/temp
                dtemp2 = dtemp2 + 1.0/temp               
		   else
              if(calc_flag(pointFaces(i,j))==2)then
			     if(Q1(pointFaces(i,j))>drydeep)then
                  c1(1) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),1)+&
                                pcoor(facePoints(pointFaces(i,j),2),1)+ &
                                pcoor(facePoints(pointFaces(i,j),3),1))
                  c1(2) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),2)+ &
                                pcoor(facePoints(pointFaces(i,j),2),2)+  &
                                pcoor(facePoints(pointFaces(i,j),3),2))
                  c1(3) = 1.0/3.0*(pcoor(facePoints(pointFaces(i,j),1),3)+ &
                                pcoor(facePoints(pointFaces(i,j),2),3)+  &
                                pcoor(facePoints(pointFaces(i,j),3),3))

                !distance vector from current node to face center
                   v2(1) = v1(1) - c1(1)  
                   v2(2) = v1(2) - c1(2)  
                   v2(3) = v1(3) - c1(3)  
                
                !distance between current node and face center
                   call vec_mag(v2, temp)
                   dtemp1 = dtemp1 + cellValues(pointFaces(i,j))/temp
                   dtemp2 = dtemp2 + 1.0/temp 
				  else
				   dtemp1 = dtemp1 + cellValues(pointFaces(i,j))
                   dtemp2 = dtemp2 + 1.0             
		          endif
			  endif
		    endif
			endif
         end do 
         
		 if(dabs(dtemp2).gt.VSMALL)then
			nodeValues(i) = dtemp1/dtemp2
		 else
		    nodeValues(i) = 0.0
		 end if
    end do

	end subroutine	


!   calculate the mean value of H, U, V in the domain
	subroutine meanHUV
	USE COMMON_MODULE,ONLY: nNodes,nodeQ1,nodeU,nodeV,nodecsed,nodedeep,nodeQ1max,nodeUMmax,nodeVNmax,meanH,meanU,meanV,meanC,meandeep,meanQ1max,meanUMmax,meanVNmax
	implicit none

	real*8 dHtemp,dUtemp,dVtemp,dCtemp,ddeeptemp,dQ1maxtemp,dUMmaxtemp,dVNmaxtemp
	integer i

    dHtemp = 0.0
    dUtemp = 0.0
    dVtemp = 0.0
    dCtemp = 0.0
	ddeeptemp = 0.0
	dQ1maxtemp = 0.0
	dUMmaxtemp = 0.0
	dVNmaxtemp = 0.0
   

    do i = 1, nNodes
		dHtemp=dHtemp+nodeQ1(i)
		dUtemp=dUtemp+nodeU(i)
		dVtemp=dVtemp+nodeV(i)  

		ddeeptemp = ddeeptemp+nodedeep(i)
	    dQ1maxtemp = dQ1maxtemp+nodeQ1max(i)
		dUMmaxtemp = dUMmaxtemp+nodeUMmax(i)
		dVNmaxtemp = dVNmaxtemp+nodeVNmax(i)
		         
	end do

	meanH = dHtemp/nNodes
	meanU = dUtemp/nNodes
	meanV = dVtemp/nNodes

	meandeep = ddeeptemp/nNodes
    meanQ1max = dQ1maxtemp/nNodes
	meanUMmax = dUMmaxtemp/nNodes
    meanVNmax = dVNmaxtemp/nNodes
	

	end subroutine


!	calculate the gradient of Q at the cell centers using Gauss's theorem
	subroutine gradientQ
	USE COMMON_MODULE,ONLY: nFaces,nNodes,nodeQ1,nodeQ2,nodeQ3,gradQ1,gradQ2,gradQ3,&
							faceEdgeNormals,faceEdges,face2DArea,edgeLength,edgePoints,&
							faceEdgesNum,CALC_FLAG,Q1,drydeep
	implicit none

	integer i,j
	integer node1,node2,node3
	real*8 temp

	gradQ1=0.0
	gradQ2=0.0
	gradQ3=0.0

	do i = 1, nFaces
	 if(Q1(i)>drydeep)then
	   do j = 1, faceEdgesNum(i)	    
		temp=(nodeQ1(edgePoints(faceEdges(i,j),1))+nodeQ1(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		gradQ1(i,1)=gradQ1(i,1)+1.0/face2DArea(i)*temp
	
		temp=(nodeQ1(edgePoints(faceEdges(i,j),1))+nodeQ1(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		gradQ1(i,2)=gradQ1(i,2)+1.0/face2DArea(i)*temp

		temp=(nodeQ2(edgePoints(faceEdges(i,j),1))+nodeQ2(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		gradQ2(i,1)=gradQ2(i,1)+1.0/face2DArea(i)*temp

		temp=(nodeQ2(edgePoints(faceEdges(i,j),1))+nodeQ2(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		gradQ2(i,2)=gradQ2(i,2)+1.0/face2DArea(i)*temp

		temp=(nodeQ3(edgePoints(faceEdges(i,j),1))+nodeQ3(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		gradQ3(i,1)=gradQ3(i,1)+1.0/face2DArea(i)*temp

		temp=(nodeQ3(edgePoints(faceEdges(i,j),1))+nodeQ3(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		gradQ3(i,2)=gradQ3(i,2)+1.0/face2DArea(i)*temp
	   end do
	 endif
	end do

	end subroutine

!	calculate the gradient of u and v at the cell centers using Gauss's theorem
	subroutine gradientUV
	USE COMMON_MODULE,ONLY: nFaces,nNodes,nodeU,nodeV,gradU,gradV,faceEdgesNum,&
							faceEdgeNormals,faceEdges,face2DArea,edgeLength,&
							edgePoints,CALC_FLAG,Q1,drydeep
	implicit none

	integer i,j
	integer node1,node2,node3
	real*8 temp

	gradU=0.0
	gradV=0.0

	do i = 1, nFaces
	 if(Q1(i)>drydeep)then
	   do j = 1, faceEdgesNum(i)	    
		temp=(nodeU(edgePoints(faceEdges(i,j),1))+nodeU(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		gradU(i,1)=gradU(i,1)+1.0/face2DArea(i)*temp
	
		temp=(nodeU(edgePoints(faceEdges(i,j),1))+nodeU(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		gradU(i,2)=gradU(i,2)+1.0/face2DArea(i)*temp

		temp=(nodeV(edgePoints(faceEdges(i,j),1))+nodeV(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		gradV(i,1)=gradV(i,1)+1.0/face2DArea(i)*temp

		temp=(nodeV(edgePoints(faceEdges(i,j),1))+nodeV(edgePoints(faceEdges(i,j),2)))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		gradV(i,2)=gradV(i,2)+1.0/face2DArea(i)*temp
	   end do
	 endif
	end do

	end subroutine

!	calculate the gradient of eta at the cell centers using Gauss's theorem
	subroutine gradientEta
	USE COMMON_MODULE,ONLY: nFaces,nNodes,nodeZsurf,gradEta,faceEdgesNum,&
							faceEdgeNormals,faceEdges,face2DArea,edgeLength,edgePoints,&
							t,BINFO,faceCenters,eta,Q1,nodeQ1,drydeep,&
                            geta,gQ1,calc_flag,pcoor,aa,nDEMPoints,wse
	implicit none

	integer i,j
	integer node1,node2,node3
	integer tindex
	real*8 temp,hB1	
	real*8 waterlevel
	real*8 Q1FN,Q2FN,Q3FN,EtaFN

		gradEta=0.0

do i = 1, nFaces
  do j = 1, faceEdgesNum(i)
	   
	 if(binfo(i,j)==2)then       
   	    
  	     
       do aa=1,nDEMPoints-1
           if(t.le.1200.0)then
               waterlevel=3.0*sin(2*3.14/1200.0*t)
            else
			   waterlevel=0.0
			endif
		 enddo
	     waterlevel=102.0+waterlevel
		 hB1= waterlevel

		 temp=(hB1*faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j)))
		 gradEta(i,1)=gradEta(i,1)+1.0/face2DArea(i)*temp
	
		 temp=(hB1*faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j)))
		 gradEta(i,2)=gradEta(i,2)+1.0/face2DArea(i)*temp
	   
	   else

		 	temp=(nodeZsurf(edgePoints(faceEdges(i,j),1))+nodeZsurf(edgePoints(faceEdges(i,j),2)))/2* &
		           faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		     gradEta(i,1)=gradEta(i,1)+1.0/face2DArea(i)*temp
	
		     temp=(nodeZsurf(edgePoints(faceEdges(i,j),1))+nodeZsurf(edgePoints(faceEdges(i,j),2)))/2* &
		           faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		     gradEta(i,2)=gradEta(i,2)+1.0/face2DArea(i)*temp
      	endif
	 end do
	end do

	end subroutine

! ----------------------------------------------------------

!   calculate the edge gradient of u and v by the area weighted average
!   of neighbour triangles
	subroutine edgeGradientUV
	USE COMMON_MODULE,ONLY: nEdges,gradU,gradV,edgeFaces,face2DArea,edgeGradU,edgeGradV
	implicit none

	integer i
	integer face1,face2
	real*8  faceArea1,faceArea2

	do i = 1, nEdges
		face1=edgeFaces(i,1)
		face2=edgeFaces(i,2)
		
		if(face1.le.0) then
			edgeGradU(i,1)=gradU(face2,1)
			edgeGradU(i,2)=gradU(face2,2)
		else if(face2.le.0) then
			edgeGradU(i,1)=gradU(face1,1)
			edgeGradU(i,2)=gradU(face1,2)
		else
			faceArea1=face2DArea(face1)
			faceArea2=face2DArea(face2)
			edgeGradU(i,1)=(gradU(face1,1)*faceArea1+gradU(face2,1)*faceArea2)&
							/(faceArea1+faceArea2)
			edgeGradU(i,2)=(gradU(face1,2)*faceArea1+gradU(face2,2)*faceArea2)&
							/(faceArea1+faceArea2)
		end if
	end do

	end subroutine


!------------------------------------

!   calculate the edge water depth (Q1)
	subroutine calc_edgeQ1
	USE COMMON_MODULE,ONLY: nEdges,edgeFaces,face2DArea,edgeQ1,Q1,Qb1
	implicit none

	integer i,pos1,pos2
	integer face1,face2
	real*8  faceArea1,faceArea2
	integer edgePositionInFace

	do i = 1, nEdges
		face1=edgeFaces(i,1)
		face2=edgeFaces(i,2)
		
		if(face1.le.0) then
			pos2=edgePositionInFace(face2,i)
			edgeQ1(i)=Qb1(face2,pos2)
		else if(face2.le.0) then
			pos1=edgePositionInFace(face1,i)
			edgeQ1(i)=Qb1(face1,pos1)
		else
			faceArea1=face2DArea(face1)
			faceArea2=face2DArea(face2)
			pos1=edgePositionInFace(face1,i)
			pos2=edgePositionInFace(face2,i)
			edgeQ1(i)=(Qb1(face1,pos1)+Qb1(face2,pos2))/2 
		end if
	end do

	end subroutine


!	calculate the face limiters for Q and eta
	subroutine calc_faceLimiters
	USE COMMON_MODULE,ONLY: nFaces,edgeFaces,faceLimiters,gradQ1,gradQ2,gradQ3,Q1,&
							Q2,Q3,faceEdges,faceCenters,VLARGE,gQ1,gQ2,gQ3,faceEdgesNum,&
							edgeCenterCoor,etaLimiter,eta,gEta,gradEta,CALC_FLAG,drydeep
	implicit none

	real*8 rtemp(3), gradTemp(3), Qedge, temp, phi(3), phimax(3), phimaxEta,phiEta
	real*8 Q1FN,Q2FN,Q3FN,EtaFN
	integer i,j

	do i=1,nFaces
		phi=0.0
		phimax=VLARGE
		phimaxEta=VLARGE
	 if(Q1(i)>drydeep)then
		do j=1,faceEdgesNum(i)
			call neighbor_value(eta,gEta,i,j,EtaFN)
			call neighbor_value(Q1,gQ1,i,j,Q1FN)
			call neighbor_value(Q2,gQ2,i,j,Q2FN)
			call neighbor_value(Q3,gQ3,i,j,Q3FN)

			!r vector: from face center to edge center
			rtemp(1)=edgeCenterCoor(faceEdges(i,j),1)-faceCenters(i,1)
			rtemp(2)=edgeCenterCoor(faceEdges(i,j),2)-faceCenters(i,2)
			rtemp(3)=0.0

			!gradEta
			gradTemp(1)=gradEta(i,1)
			gradTemp(2)=gradEta(i,2)
			gradTemp(3)=0.0
			call vec_dot(gradTemp,rtemp,temp)
			!unlimited value
			Qedge=Eta(i) + temp

			call limiter(Eta(i),Qedge,EtaFN,phiEta)

			if(phimaxEta>phiEta)then
				phimaxEta=phiEta
			end if

			!gradQ1
			gradTemp(1)=gradQ1(i,1)
			gradTemp(2)=gradQ1(i,2)
			gradTemp(3)=0.0
			call vec_dot(gradTemp,rtemp,temp)
			!unlimited value
			Qedge=Q1(i) + temp

			call limiter(Q1(i),Qedge,Q1FN,phi(1))

			if(phimax(1)>phi(1))then
				phimax(1)=phi(1)
			end if

			!gradQ2
			gradTemp(1)=gradQ2(i,1)
			gradTemp(2)=gradQ2(i,2)
			gradTemp(3)=0.0
			call vec_dot(gradTemp,rtemp,temp)
			!unlimited value
			Qedge=Q2(i) + temp

			call limiter(Q2(i),Qedge,Q2FN,phi(2))

			if(phimax(2)>phi(2))then
				phimax(2)=phi(2)
			end if

			!gradQ3
			gradTemp(1)=gradQ3(i,1)
			gradTemp(2)=gradQ3(i,2)
			gradTemp(3)=0.0
			call vec_dot(gradTemp,rtemp,temp)
			!unlimited value
			Qedge=Q3(i) + temp

			call limiter(Q3(i),Qedge,Q3FN,phi(3))

			if(phimax(3)>phi(3))then
				phimax(3)=phi(3)
			end if
		end do
		etaLimiter(i)=max(0.0,phimaxEta)
		faceLimiters(i,1)=max(0.0,phimax(1))
		faceLimiters(i,2)=max(0.0,phimax(2))
		faceLimiters(i,3)=max(0.0,phimax(3))
	 endif
	end do
	
	end subroutine

!   update the values of Q for ghost cells
	subroutine update_ghostcell_value
	USE COMMON_MODULE,ONLY:nBoundaryEdges,ghostCellsNeighbor,gZB,gEta,gQ1,faceCenters
	implicit none
	
	integer i

	!for ghost cells
	do i=1,nBoundaryEdges
		gZB(i)=faceCenters(ghostCellsNeighbor(i),3)
		gEta(i)=gQ1(i)+gZB(i)
	enddo

	end subroutine

!	calculate the Courant number
	subroutine calc_Courant
	USE COMMON_MODULE,ONLY:u,v,alpha,faceCourantNumber,faceMinR,&
						   maxCourantNumber,minCourantNumber,nFaces,&
						   VSMALL,VLARGE,faceEdgesNum,dt
	implicit none
	
	real*8 maxSpeed
	integer i,j

	maxCourantNumber = VSMALL
	minCourantNumber = VLARGE
	do i=1, nFaces 
	   maxSpeed = VSMALL
	   do j=1, faceEdgesNum(i)
		  maxSpeed = max(maxSpeed, dsqrt(u(i,j)**2+v(i,j)**2)+alpha(i,j))
	   enddo
	   faceCourantNumber(i)=dt*maxSpeed/faceMinR(i)
	   maxCourantNumber = max(maxCourantNumber,faceCourantNumber(i))
	   minCourantNumber = min(minCourantNumber,faceCourantNumber(i))
	enddo

!	write(*,*) 'Max Courant=', maxCourantNumber
!	write(*,*) 'Min Courant=', minCourantNumber

	end subroutine


!   calculate the edge water depth (Q1)
	subroutine wetdry_edge
	USE COMMON_MODULE,ONLY: nEdges,edgeFaces,face2DArea,edgeQ1,Q1,Qb1,edgedrywet,&
	                        drydeep,eta
	implicit none

	integer i,pos1,pos2
	integer face1,face2
	real*8  faceArea1,faceArea2
	integer edgePositionInFace

	do i = 1, nEdges
		face1=edgeFaces(i,1)
		face2=edgeFaces(i,2)
		
		if(face1.le.0) then
           if(Q1(face2).gt.drydeep)then
             edgedrywet(i)=2
		   else
		     edgedrywet(i)=0
		   endif
		else if(face2.le.0) then
           if(Q1(face1).gt.drydeep)then
             edgedrywet(i)=2
		   else
		     edgedrywet(i)=0
		   endif
		else
			pos1=edgePositionInFace(face1,i)
			pos2=edgePositionInFace(face2,i)
           if(Q1(face1).gt.drydeep.and.Q1(face2).gt.drydeep)then
                edgedrywet(i)=2

		   elseif(Q1(face1).gt.drydeep.and.Q1(face2).le.drydeep)then
		       if(eta(face1).gt.eta(face2))then
                  edgedrywet(i)=1
			   else
                  edgedrywet(i)=0
			   endif
		   elseif(Q1(face1).le.drydeep.and.Q1(face2).gt.drydeep)then
		       if(eta(face1).lt.eta(face2))then
                  edgedrywet(i)=1
			   else
                  edgedrywet(i)=0
			   endif
		   elseif(Q1(face1).le.drydeep.and.Q1(face2).le.drydeep)then
              edgedrywet(i)=0
	
		   endif

		end if
	end do

	end subroutine
!*****************************************************************************


