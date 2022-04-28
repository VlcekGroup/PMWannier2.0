subroutine awf_prep
      use OMP_LIB
      use commvar
      implicit none
      integer::l,i,j,k
      real*8::dis_sqr=0d0,sum_itmd=0d0

      !grid point coordinate
      l=1
      do i=1,nz
        do j=1,ny
          do k=1,nx
            GR(l,1) = (k-nx/2-1)*dx
            GR(l,2) = (j-ny/2-1)*dy
            GR(l,3) = (i-nz/2-1)*dz
            l=l+1
          end do
        end do
      end do

     !density of atom
     !$OMP parallel firstprivate(dis_sqr,sum_itmd)
     !$OMP do
      do j=1,nx*ny*nz
        do i=1,noa
          do k=1,3
            dis_sqr = dis_sqr+(GR(j,k)-COA(i,k))**2.0
          end do
          if(dsqrt(dis_sqr).lt.rcut) then
            DOA(j,i) = (Nel(i)/(gamma_a*sqrt_2pi))*exp(-dis_sqr/(2.0*gamma_a**2.0))
          else
            DOA(j,i) = 0d0
          end if
            dis_sqr = 0d0
        end do
        do k=1,noa
          sum_itmd = sum_itmd + DOA(j,k)
        end do
        do k=1,noa
          if(DOA(j,k).gt.0d0) then
           AWF(j,k) = DOA(j,k)/sum_itmd
          else
           AWF(j,k) = 0d0
          end if
        end do
        sum_itmd = 0d0
      end do
      !$OMP end do
      !$OMP end parallel

end subroutine awf_prep
