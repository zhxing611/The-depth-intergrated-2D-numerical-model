!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   sediment transport: one step
	subroutine sediment()
	USE COMMON_MODULE
	implicit none
  
    integer i,j
    integer*4 dryindex,wetindex
	real*8  Q1FN,EtaFN


	 call csed_Limiter 		
		      
	 call calc_sedimentFlux
	 call calc_divSedimentFlux		
	

	do i=1,nFaces			
		curFace=i
				
		if(Q1(i)>drydeep)then
	       call calc_FV4
           call calc_FdotN4      		
		end if
	end do
		 
	do i=1,nFaces			
		
		curFace=i		
		if(Q1(i)>drydeep)then
           call calc_sedRofQi
		end if
	end do
    
	call calc_integrate


	!update the mesh information
	call updateMeshData

    call cellCenterToNodes(csed,nodecsed)
	    
	return
	end

!   read sediment properties

!xu-----------------------------------------------------      
!   calculate sediment flux
	subroutine calc_sedimentFlux()
	USE COMMON_MODULE
	implicit none

	integer i,face1,face2,j,pos,ghostCell
	real*8 dzdx,dzdy,temp,dtemp1,dtemp2
	real*8 v1(2),c1(3),c2(3),edgeCenter(3)
    real*8 v4(3),curEdgeNormal(3)
	real*8 dist2points
	real*8 utemp,vtemp,temp1
	real*8 csedFN,Qb2FN,Qb3FN
	integer edgePositionInFace
	
	  
    do i =1, nFaces
   
	 do j = 1, faceEdgesNum(i)

	   	  curEdge = faceEdges(i, j)

		 call neighbor_value(csed,gcsed,i,j,csedFN)

		 UQ(i,j)=Qb2(i,j)*faceEdgeNormals(i,j,1) +Qb3(i,j)*faceEdgeNormals(i,j,2)
       
	    if(UQ(i,j).gt.0.0)then
		   CBFlux(i,j,1)=Qb2(i,j)*(csed(i)+0.5*LimiterR(curEdge)*(csedFN-csed(i)))
           CBFlux(i,j,2)=Qb3(i,j)*(csed(i)+0.5*LimiterR(curEdge)*(csedFN-csed(i)))
		else
           CBFlux(i,j,1)=Qb2(i,j)*(csedFN-0.5*LimiterR(curEdge)*(csedFN-csed(i)))
           CBFlux(i,j,2)=Qb3(i,j)*(csedFN-0.5*LimiterR(curEdge)*(csedFN-csed(i)))
		endif
		   		 
		 if(faceNeighbors(i,j).lt.0) then !boundary cells
			ghostCell=boundaryEdgeGhostCells(faceEdges(i,j))
  
              CBFlux(i,j,1)=gQb2(ghostCell)*gcsed(ghostCell)
              CBFlux(i,j,2)=gQb3(ghostCell)*gcsed(ghostCell)
    
         endif

	  enddo
	enddo
    
	return
	end
!************************************************************************
!****************************************************************************

!c  This subroutine calculates the intercell flux based on Upwind funtion
!	calculate divergence of sediment flux at face centers
	subroutine calc_divSedimentFlux()
	USE COMMON_MODULE
	implicit none

	integer i,j,pos,k
    real*8 c1(3),v1(3),v2(3),v3(3),v4(3),curEdgeNormal(3)
	real*8 temp
	integer edgePositionInFace

    do i = 1, nFaces
         FdotCFlux(i) = 0.0

         ! for each edge
         do j = 1, faceEdgesNum(i)
            curEdge = faceEdges(i, j)

			pos=edgePositionInFace(i,curEdge)
			if(pos==-1)then
				write(*,*) 'Something wrong!'
				stop
			else
				curEdgeNormal(1) = faceEdgeNormals(i,pos,1)
				curEdgeNormal(2) = faceEdgeNormals(i,pos,2)
				curEdgeNormal(3) = 0.0	
							
			end if

            !current edge sediment flux vector(expand to 3D vector for vec_dot)
            v4(1) =CBFlux(i,j,1)   !x    CBFlux(curEdge, 1)
            v4(2) =CBFlux(i,j,2)   !y    CBFlux(curEdge, 2)
            v4(3) = 0.0                          !z

            !contribution of divergence from current edge integration
            call vec_dot(v4, curEdgeNormal, temp)
		
            FdotCFlux(i) = FdotCFlux(i) + temp*edgeLength(curEdge) 
 
         end do      
	  end do

	return
	end
!****************************************************************************
	subroutine calc_FV4()
	USE COMMON_MODULE
	implicit none

	integer j,edgeNum
	real*8 ht
	real*8 nx,ny
	
	do j=1,faceEdgesNum(curFace)
	  
	  curEdge = faceEdges(curFace, j)
      
	  if(Qb1(curFace,j)>drydeep)then
	  
	      visc=0.3*Qb1(curFace,j)*dsqrt((u(curFace,j)**2+v(curFace,j)**2))
		
    	  ht=Qb1av(curFace,j)
		
		  nx=faceEdgeNormals(curFace,j,1)
		  ny=faceEdgeNormals(curFace,j,2)
	
		  edgeNum=faceEdges(curFace,j)

		  FV4(curFace,j)=visc*(ht*edgeGcsed(edgeNum,1)*nx+ht*edgeGcsed(edgeNum,2)*ny)   !

	  end if
	end do

	end subroutine
!*********************************************************************
	subroutine calc_FdotN4
	USE COMMON_MODULE,ONLY: FdotN1,FdotN2,FdotN3,FI1,FI2,FI3,nFaces,curFace,&
		FIofQb1,FIofQb2,FIofQb3,FV1,FV2,FV3,sedimentctl,t,dt,ts,td,FdotN4,FV4,&
		scalc,sedInterval,binfo,edgeLength,faceEdges,ELEDGES,faceEdgesNum,gradcsed
	implicit none

	real*8 length(ELEDGES)
	integer i,j
	
	do i=1,faceEdgesNum(curFace)
		length(i)=edgeLength(faceEdges(curFace,i))
	end do

	FdotN4(curFace)=0.0

	do i=1,ELEDGES
		FdotN4(curFace)=FdotN4(curFace)+FV4(curFace,i)*length(i)

	end do

	return       
	end
 
!********************************************************************
	subroutine calc_sedRofQi()
	USE COMMON_MODULE,ONLY: maxfaces_,Rem4,Rem2,Rem3,oldRem4,oldRem2,oldRem3,Qs,&
		FdotN4,FdotN2,FdotN3,HiVi1,HiVi2,HiVi3,curFace,t,dt,ts,td,drydeep,mindeep&
		,Q1,Q2,Q3,sedimentctl,scalc,sedInterval,nFaces,bednet,face2DArea,FdotCFlux
	implicit none
   

	 Rem4(curFace)=(-FdotCFlux(curFace)+FdotN4(curFace))*(dt/face2DArea(curFace))  !+ &   !
	                !Qs(curFace)*dt   !  
 
	 if(Q1(curFace)<=drydeep)then  
		Rem4(curFace)=0.0		
				
	 end if

	return
	end

!*******************************************************************
	subroutine calc_integrate()
	USE COMMON_MODULE,ONLY: Q1,Q2,Q3,face2DArea,curFace,a,b,d,&
		Rem1,Rem2,Rem3,oldRem1,oldRem2,oldRem3,&
		t,dt,ts,td,nFaces,Sox,Soy,&
		tauwx,tauwy,Swx,Swy,nb,coarse,frctl,g,drydeep,mindeep,&
		sedimentctl,scalc,sedInterval,Rem4,VSMALL,poro,&
		faceNeighbors,binfo,faceEdgesNum,csed,FdotCFlux,FdotN4

	implicit none

	integer i
  
    do i = 1, nFaces


	 if(Q1(i)>drydeep)then

	   csed(i)=csed(i)*Q1(i)+ Rem4(i)    
	   csed(i)=csed(i)/Q1(i)
       csed(i)=max(0.0,csed(i))
     endif

    enddo

	return
	end

!********************************************************************

!*************************************************
!计算界面垂向速度UB1-left垂向，UB2-right垂向

    subroutine csed_Limiter

	USE COMMON_MODULE,ONLY:nFaces,g,u,v,faceEdges,nEdges,edgeFaces,&
		curFace,Q2in,FI1,FI2,FI3,FIofQb1,FIofQb2,FIofQb3,&
		modA,R,modLAMBDA,L,faceEdgeNormals,t,dt,ts,td,boundaryEdgeGhostCells,&
		Qb1,Qb2,Qb3,binfo,drydeep,mindeep,inctr,outctr,inletH,gUQ,&
		gFIofQb1,gFIofQb2,gFIofQb3,gQb1,gQb2,gQb3,faceEdgesNum,UQ,&
		eta,gEta,faceCenters,gQ1,gQ2,gQ3,faceNeighbors,LimiterR,&
		csed,gcsed,csedLimiter,edgeLimiter,csedR,edgeLength,edgePoints,facePoints,&
		pcoor,nodecsed
	
	implicit none
	integer i,j,k,pointk
	real*8 Qb1FN,Qb2FN,Qb3FN,csedFN
    real*8  temp,temp1,temp2,temp3
    real*8  c1(3), c2(3), v1(3)
	integer pos1,pos2
	integer face1,face2
	integer edgePositionInFace
!   Functions out side of main program file
    real*8 dist2points 
			

  !interpolate the face center values to edges
     do i = 1, nEdges

		face1=edgeFaces(i,1)
		face2=edgeFaces(i,2)
		
		if(face1.le.0) then
			pos2=edgePositionInFace(face2,i)
            edgeLimiter(i) = csedLimiter(face2)
			edgeLimiter(i) = 1.0
		else if(face2.le.0) then
			pos1=edgePositionInFace(face1,i)
			edgeLimiter(i) = csedLimiter(face1)
			edgeLimiter(i) = 1.0			
		else
			
			pos1=edgePositionInFace(face1,i)
			pos2=edgePositionInFace(face2,i)
 !判断流量流出单元
            temp=Qb2(face1,pos1)*faceEdgeNormals(face1,pos1,1)+ &
			         Qb3(face1,pos1)*faceEdgeNormals(face1,pos1,2)

!边界处变化梯度R--            
            if(temp.gt.0.0)then

			  if(facePoints(face1,1).ne.edgePoints(i,1).and.facePoints(face1,1).ne.edgePoints(i,2))then
			           pointk=facePoints(face1,1)
			  elseif(facePoints(face1,2).ne.edgePoints(i,1).and.facePoints(face1,2).ne.edgePoints(i,2))then
			           pointk=facePoints(face1,2)
              elseif(facePoints(face1,3).ne.edgePoints(i,1).and.facePoints(face1,3).ne.edgePoints(i,2))then
			           pointk=facePoints(face1,3)
              endif

 !face1 face2 he pointk de 距离			 
                   c1(1)=faceCenters(face1,1)	
        		   c1(2)=faceCenters(face1,2)
        		   c1(3)=0.0

                   c2(1)=faceCenters(face2,1)	
        		   c2(2)=faceCenters(face2,2)
        		   c2(3)=0.0

                   v1(1)=pcoor(pointk,1)	
        		   v1(2)=pcoor(pointk,2)
        		   v1(3)=0.0

                edgeLimiter(i)=(csed(face1)-nodecsed(pointk))/(csed(face2)-csed(face1))*dist2points(c1,c2)/dist2points(c1,v1)
                            
			else

			  if(facePoints(face2,1).ne.edgePoints(i,1).and.facePoints(face2,1).ne.edgePoints(i,2))then
			           pointk=facePoints(face2,1)
			  elseif(facePoints(face2,2).ne.edgePoints(i,1).and.facePoints(face2,2).ne.edgePoints(i,2))then
			           pointk=facePoints(face2,2)
              elseif(facePoints(face2,3).ne.edgePoints(i,1).and.facePoints(face2,3).ne.edgePoints(i,2))then
			           pointk=facePoints(face2,3)
              endif

 !face1 face2 he pointk de 距离			 
                   c1(1)=faceCenters(face1,1)	
        		   c1(2)=faceCenters(face1,2)
        		   c1(3)=0.0

                   c2(1)=faceCenters(face2,1)	
        		   c2(2)=faceCenters(face2,2)
        		   c2(3)=0.0

                   v1(1)=pcoor(pointk,1)	
        		   v1(2)=pcoor(pointk,2)
        		   v1(3)=0.0

                edgeLimiter(i)=(csed(face2)-nodecsed(pointk))/(csed(face1)-csed(face2))*dist2points(c1,c2)/dist2points(c2,v1)
			endif  						         
		  endif

          LimiterR(i)=max(0.0,min(1.0,2.0*edgeLimiter(i)),min(2.0,edgeLimiter(i)))  

      enddo

      end subroutine

!---------------------------------------------------------
 