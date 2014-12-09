
c-----------------------------------------------------------------------
      subroutine bxinv (m, n, a, b, c)
      implicit none
      integer m, n
      real*8 a(n,n), b(m,n), c(m,n)
c
c     does the operation c :=  - b * inv(a)
c     where a has already been factored by Gauss.
c
      integer i, j, k
      real*8 sum
c
c     c = b (LU)**(-1) = b U\inv x L \inv
c     get c := b U \inv
c
      do i=1, m
c
c     U solve -- solve for row number k :   c(k,*) U = b(k,*)
c
         c(i,1) = - b(i,1) * a(1,1)
         do j=2, n
            sum = - b(i,j)
            do k=1, j-1
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum*a(j,j)
         enddo
      enddo
c
      do i=1, m
c
c     L- solve -- solve for row number i :   c(i,*) U = b(i,*)
c
         do j=n-1, 1, -1
            sum = c(i,j)
            do k=j+1, n
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum
         enddo
      enddo
c    
      end
c--------------

c-----------------------------------------------------------------------
      subroutine zbxinv (m, n, a, b, c)
      implicit none
      integer m, n
c      real*8 a(n,n), b(m,n), c(m,n)
      double complex a(n,n), b(m,n), c(m,n)
c
c     does the operation c :=  - b * inv(a)
c     where a has already been factored by Gauss.
c
      integer i, j, k
c      real*8 sum
      double complex sum
c
c     c = b (LU)**(-1) = b U\inv x L \inv
c     get c := b U \inv
c
      do i=1, m
c
c     U solve -- solve for row number k :   c(k,*) U = b(k,*)
c
         c(i,1) = - b(i,1) * a(1,1)
         do j=2, n
            sum = - b(i,j)
            do k=1, j-1
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum*a(j,j)
         enddo
      enddo
c
      do i=1, m
c
c     L- solve -- solve for row number i :   c(i,*) U = b(i,*)
c
         do j=n-1, 1, -1
            sum = c(i,j)
            do k=j+1, n
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum
         enddo
      enddo
c
      end
c--------------
        subroutine qsplit(a,ind,n,ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end

      subroutine gauss (n,a,ierr)
c-----------------------------------------------------------------------
      implicit none
      integer n, ierr
      real*8 a(n,n)
c
c     does the Gaussian factorization a := LU
c
      integer i, j, k
      real*8 piv
c-----------------------------------------------------------------------
      ierr = 0
      do k=1, n
         if (a(k,k) .eq. 0.0) then
            ierr = 1
            return
         endif
c
         a(k,k) = 1.0/a(k,k)
         do i=k+1, n
            piv = a(i,k) * a(k,k)
            do j=k+1, n
               a(i,j) = a(i,j) - piv*a(k,j)
            enddo
            a(i,k) = piv
         enddo
      enddo
      return
      end

      subroutine zgauss (n,a,ierr)
c-----------------------------------------------------------------------
      implicit none
      integer n, ierr
c      real*8 a(n,n)
      double complex a(n,n)
c
c     does the Gaussian factorization a := LU
c
      integer i, j, k
c      real*8 piv
      double complex piv
c-----------------------------------------------------------------------
      ierr = 0
      do k=1, n
         if (a(k,k) .eq. 0.0) then
            ierr = 1
            return
         endif
c
         a(k,k) = 1.0/a(k,k)
         do i=k+1, n
            piv = a(i,k) * a(k,k)
            do j=k+1, n
               a(i,j) = a(i,j) - piv*a(k,j)
            enddo
            a(i,k) = piv
         enddo
      enddo
      return
      end
