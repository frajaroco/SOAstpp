C     Francisco J. Rodriguez-Cortes, September 2016

C     This function provides an edge-corrected estimate
C     of the anisotropic spatio-temporal inhomogeneous K-function.
      
      subroutine astkfunct(x,y,txy,n,lambda,ag,s,ns,t,nt,vol,
     +      wbi,wbimod,wt,correc,astkf)

      implicit real*8(a-h,o-z)

      integer n,ns,nt,iu,iv,i,j,k,correc
      double precision astkf,two,vij,wij,hij,tij,angij,ag,xi,yi,ti
      double precision x,y,txy,xc,yc,wbi,wbimod,wt,lambda
      dimension wbi(n,ns,nt),wbimod(n,ns,nt),wt(n,n),correc(4),lambda(n)
      dimension x(n),y(n),txy(n),s(ns),t(nt),astkf(ns,nt),xc(n),yc(n)

      two=2d0      
      pi=3.14159265d0
      do iu=1,ns
      do iv=1,nt
       do i=1,n
       xi=x(i)
       yi=y(i)
       ti=txy(i)

      do k=1,n
       xc(k)=x(k)-xi 
       yc(k)=y(k)-yi
      end do
         do j=1,n
          if (j.ne.i) then
            hij=dsqrt((xc(j))*(xc(j)) + (yc(j))*(yc(j)))
            tij=dabs(ti-txy(j))
C Quadrant one
            if ((xc(j).gt.0d0).and.(yc(j).gt.0d0)) then
               angij=datan(yc(j)/xc(j))
C Quadrant two		   
             else if ((xc(j).lt.0d0).and.(yc(j).gt.0d0)) then
               angij=datan(yc(j)/xc(j))+pi
C Quadrant three
             else if ((xc(j).lt.0d0).and.(yc(j).lt.0d0)) then
               angij=datan(yc(j)/xc(j))+pi
C Quadrant four
             else if ((xc(j).gt.0d0).and.(yc(j).lt.0d0)) then   
               angij=datan(yc(j)/xc(j))+(two*pi)
      end if
      if ((tij.le.t(iv)).and.(hij.le.s(iu)).and.(angij.le.ag)) then            
C none
      if (correc(1).eq.1) then
         vij=1d0
         wij=vij/(vol*lambda(i)*lambda(j))
         astkf(iu,iv)=astkf(iu,iv)+wij
      end if
C border
      if (correc(2).eq.1) then
         vij=wbi(i,iu,iv)
         wij=vij/(vol*lambda(i)*lambda(j))
         astkf(iu,iv)=astkf(iu,iv)+wij
      end if
C modified border
      if (correc(3).eq.1) then
         vij=wbimod(i,iu,iv)
         wij=vij/(vol*lambda(i)*lambda(j))
         astkf(iu,iv)=astkf(iu,iv)+wij
      end if
C translate
      if (correc(4).eq.1) then
         vij=wt(i,j)
         wij=vij/(vol*lambda(i)*lambda(j))
         astkf(iu,iv)=astkf(iu,iv)+wij
      end if
           end if          
             end if
        
        end do
        end do
        end do
        end do
      return
      end
