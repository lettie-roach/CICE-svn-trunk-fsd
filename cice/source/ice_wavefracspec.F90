! Lettie Roach, NIWA/VUW, June 2016
! Some of this modified from Chris Horvat's Matlab original
! New module for CICE to generate waves and break the FSD
!
      module ice_wavefracspec

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nfsd, ncat, max_ntrcr, max_blocks
 
      implicit none
      private
      public :: wave_frac, wave_frac_fsd 



      real, parameter  :: &
                swh_minval = 0.01_dbl_kind, & ! minimum value of wave height (m)
                straincrit =0.00003,        & ! critical strain
                D=10000.,                   & ! domain size
                dx = c1,                    & ! domain spacing
                threshold=10                  ! peak-finding threshold 
                                              ! (points are defined to be extrema
                                              ! if they are a local maximum or minimum
                                              ! over a distance of 10m on both sides, 
                                              ! based on the observations of Toyota et al. 
                                              ! (2011) who find this to be the order of
                                              ! the smallest floe size affected by wave fracture 
      integer (kind=int_kind) :: &
          nx = 10000                          ! number of points in domain
 
      logical (kind=log_kind), public :: &
         wave_spec      

!=======================================================================

      contains

!=======================================================================

     function get_damfstd_wave(amfstd_init, fracture_hist, frac) &
                              result(d_amfstd)

      real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         amfstd_init, fracture_hist

     real (kind=dbl_kind), dimension (nfsd,nfsd), intent(in) :: &
         frac

     real (kind=dbl_kind), dimension (nfsd) :: &
         d_amfstd, loss, gain, omega
 
      integer (kind=int_kind) :: k



      do k = 1,nfsd
           ! fracture_hist is already normalized
           omega(k) = amfstd_init(k)*SUM(fracture_hist(1:k-1)) 
      end do

      if (SUM(omega).gt.c1+puny) stop &
          'omega cannot be greater than 1, waves'
                        
      loss = omega

      do k =1,nfsd
           gain(k) = SUM(omega*frac(:,k)) 
      end do

      if (gain(nfsd).gt.puny) stop 'largest cat cannot gain, waves'
      if (loss(1).gt.puny) stop 'smallest cat cannot lose, waves'

       d_amfstd(:) = gain(:) - loss(:)
 
      if (SUM(d_amfstd(:)).gt.puny) stop 'area not cons, waves'



      end  function get_damfstd_wave

!=======================================================================

     function get_subdt_wave(amfstd_init, d_amfstd) &
                              result(subdt)

      real (kind=dbl_kind), dimension (nfsd), intent(in) :: &
         amfstd_init, d_amfstd

      real (kind=dbl_kind), dimension (nfsd) :: &
         check_dt 
 
      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: subdt

      check_dt(:) = bignum 
      do k = 1, nfsd
          if (d_amfstd(k).gt.puny) check_dt(k) = (1-amfstd_init(k))/d_amfstd(k)
          if (d_amfstd(k).lt.-puny) check_dt(k) = amfstd_init(k)/ABS(d_amfstd(k))
      end do 
            
      subdt = MINVAL(check_dt)
 

      end function get_subdt_wave

!=======================================================================
!  Author: Lettie Roach, NIWA, 2018
! 
! Given fracture histogram computed from local wave spectrum, evolve FSD

     subroutine wave_frac_fsd

     use ice_calendar, only: dt
     use ice_communicate, only: my_task
     use ice_domain, only: nblocks, blocks_ice
     use ice_blocks, only: block, get_block, nx_block, ny_block
     use ice_state, only: nt_fsd, aice, aicen, vice, &
                          vicen, trcrn, nt_qice
     use ice_flux, only: wave_spectrum, wave_hs_in_ice, dfreq
     use ice_fsd, only: d_afsd_wave, d_amfstd_wave
 
     ! local variables
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &  
        iblk, n, k, i, j, ks, t, &
        ilo, ihi, jlo, jhi, & ! beginning and end of physical domain
        nsubt ! number of subcycles 

      real (kind=dbl_kind), dimension (nfsd) :: &
        fracture_hist

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
        frac    

      real (kind=dbl_kind) :: &
        hbar, elapsed_t, subdt

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat) :: &
           afstd, &           ! joint floe size and ice thickness area distribution 
           amfstd_init, &     ! tracer array
           afstd_init         ! aicen*trcrn

      real (kind=dbl_kind), dimension (nfsd) :: &
         amfstd_tmp, d_amfstd_tmp, &
         stability_test, &   ! check timestep stability
         omega, &
         check_dt, &
         gain, loss, &
         tempfracs

     !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,hbar,fracture_hist,amfstd_init,amfstd_tmp,d_amfstd_tmp,n,k,ks,loss,gain,omega,check_dt,nsubt,subdt,elapsed_t)
     do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi


            d_afsd_wave(i,j,:,iblk) = c0

            do n=1, ncat
                        if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk).lt.c0)) stop 'neg b4-wb'
            end do


            if ((aice(i,j,iblk).gt.p01).and.(.NOT.(ALL(wave_spectrum(i,j,:,iblk).lt.puny)))) then

                ! save for diagnostics
                wave_hs_in_ice(i,j,iblk) = c4*SQRT(SUM(wave_spectrum(i,j,:,iblk)*dfreq))

                hbar = vice(i,j,iblk) / aice(i,j,iblk)

                call  wave_frac(hbar, wave_spectrum(i,j,:,iblk), fracture_hist)

                if (.not. ALL(fracture_hist.lt.puny)) then

                    do n = 1, ncat
                    if ((aicen(i,j,n,iblk).gt.puny).and.&
                       (SUM(trcrn(i,j,nt_fsd+1:nt_fsd+nfsd-1,n,iblk)).gt.puny).and.&
                       (.NOT.(trcrn(i,j,nt_fsd+1,n,iblk).ge.c1))) then
                       
                        amfstd_init(i,j,:,n)=trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk)

                        if (ABS(SUM(amfstd_init(i,j,:,n))-c1).gt.puny) stop &
                                'init mFSTD not norm, wave'
                        
                        ! protect against small numerical errors
                        WHERE (amfstd_init(i,j,:,n).lt.puny) amfstd_init(i,j,:,n) = c0
                        amfstd_init(i,j,:,n) = amfstd_init(i,j,:,n) / SUM(amfstd_init(i,j,:,n))

                        amfstd_tmp =  amfstd_init(i,j,:,n)

                       ! frac does not vary within subcycle
                        frac(:,:) = c0
                        do k = 2, nfsd
                                frac(k,:k-1) = fracture_hist(:k-1)
                                frac(k,k:) = c0 
                        end do
                        
                        do ks=1,nfsd
                            if (SUM(frac(ks,:)).gt.c0) frac(ks,:) = frac(ks,:)/SUM(frac(ks,:))
                        end do

                        ! adaptive sub-timestep
                        elapsed_t = c0
                        DO WHILE (elapsed_t.lt.dt)

                             ! calculate d_amfstd using current afstd
                             d_amfstd_tmp = get_damfstd_wave(amfstd_tmp, fracture_hist, frac)
                             WHERE (ABS(d_amfstd_tmp).lt.puny) d_amfstd_tmp = c0 
 
                             ! timestep required for this
                             subdt = get_subdt_wave(amfstd_tmp, d_amfstd_tmp)
                             subdt = MIN(subdt, dt)

                             ! update amfstd and elpased time
                             amfstd_tmp = amfstd_tmp + subdt * d_amfstd_tmp(:)

                             ! check conservation and negatives
                            if (ANY(amfstd_tmp.lt.-puny)) stop &
                                     'wb, <0 in loop'

                             if (ANY(amfstd_tmp.gt.c1+puny)) stop &
                                     'wb, >1 in loop'

                             ! in case of small numerical errors
                             WHERE (amfstd_tmp.lt.puny) amfstd_tmp = c0
                             amfstd_tmp = amfstd_tmp/SUM(amfstd_tmp)

                            ! update time
                             elapsed_t = elapsed_t + subdt 

                        END DO

                        ! was already normalized in loop
                        trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk) = amfstd_tmp
     
                        if (ABS(SUM(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk))-c1).gt.puny) stop 'not 1 wb'
                        if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk).lt.c0)) stop 'neg wb'
                        if (ANY(trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk).gt.c1)) stop '>1 wb'

                        ! for diagnostics
                        d_amfstd_wave(i,j,:,n,iblk) = trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,n,iblk) - amfstd_init(i,j,:,n)  
                        d_afsd_wave(i,j,:,iblk) =  d_afsd_wave(i,j,:,iblk) + aicen(i,j,n,iblk)*d_amfstd_wave(i,j,:,n,iblk)

                    end if ! aicen>puny
                    end do ! n

                end if ! not all frac zero
            end if ! aice>p01

         end do ! i
         end do  !j

     end do ! iblk
     !$OMP END PARALLEL DO


     end subroutine wave_frac_fsd


!=======================================================================
! Author: Lettie Roach, NIWA, 2018
!
! Based on MatLab code from Horvat & Tziperman (2015). Calculates functions 
! to describe the change in the FSD when waves fracture ice, given a wave
! spectrum (1D frequency, 25 frequency bins) in ice. 
! We calculate extrema and if these are successive maximum, minimum, maximum or
! vice versa, and have strain greater than a critical strain, break ice and 
! create new floes with lengths equal to these distances

     subroutine wave_frac(hbar, spec_efreq, frac_local)

     use ice_communicate, only: my_task, master_task
     use ice_fsd, only: floe_rad_l, floe_binwidth, floe_rad_c, floe_rad_h
     use ice_domain_size, only: nfsd
     use ice_flux, only: freq, dfreq
 
     real (kind=dbl_kind),  intent(in) :: &
         hbar   ! mean ice thickness
! TEMP INOUT
     real (kind=dbl_kind), dimension(25), intent(inout) :: &
         spec_efreq
 
     real (kind=dbl_kind), dimension (nfsd), intent(out) :: &
         frac_local

     ! local variables

     integer (kind=int_kind) :: i, j, k

     integer, parameter :: &
         loopcts = 1    ! number of SSH realizations

     real (kind=dbl_kind), dimension(25) :: &
         spec_elambda, &
         reverse_spec_elambda, &
         reverse_lambda, lambda,&! wavelengths (m)
         reverse_dlambda, dlambda, &
         spec_coeff, &
         phi, rand_array, summand

     real (kind=dbl_kind), dimension(:), allocatable :: &
         fraclengths

     real (kind=dbl_kind), dimension(nx) :: &
          X, eta

     real (kind=dbl_kind), dimension(nfsd) :: &
                frachistogram 

     logical (kind=log_kind) :: &
           e_stop          ! if true, stop and return zero omega and fsdformed


     ! switch that aborts fracture calc if true
     e_stop=.false.
 
     ! spatial domain
     do j=1,nx
        X(j)= j*dx
     end do


     ! dispersion relation
     reverse_lambda(:) = gravit/(c2*pi*freq(:)**c2)
     reverse_dlambda(:) = gravit/(c2*pi*dfreq(:)**c2)

     ! convert to lambda spectrum
     reverse_spec_elambda(:) = spec_efreq(:) * &
                          (p5 * (gravit/(2*pi*reverse_lambda(:)**c3) )**p5 )    
 
     ! reverse lambda
     lambda(:) = reverse_lambda(25:1:-1)
     dlambda(:) = reverse_dlambda(25:1:-1)
     spec_elambda(:) = reverse_spec_elambda(25:1:-1) 
 
     ! spectral coefficients
     spec_coeff = sqrt(c2*spec_elambda*dlambda)
     if (ANY(ISNAN(spec_coeff))) stop 'NaN spec_coeff' 

     ! initialize array to hold fracture lengths
     allocate(fraclengths(1))
     fraclengths(1)=c0
     
     ! loop over n. realizations of SSH
     do i=1,loopcts

         ! random phase for each Fourier component
         ! varies in each j loop
         call RANDOM_NUMBER(rand_array)
         phi = c2*pi*rand_array
 
         do j=1,nx
             !SSH field in space (sum over wavelengths, no attenuation)
             summand = spec_coeff*COS(2*pi*X(j)/lambda+phi)
             eta(j)=SUM(summand)
         end do
         
         if ((.NOT.(ALL(eta.eq.c0))).and.(hbar.gt.puny)) then 
            !if (my_task.eq.32) print *, 'now get frac lengths'
            call get_fraclengths(X, eta, fraclengths, hbar , e_stop)
            !if (my_task.eq.32) print *, 'got fraclengths'
         end if
     end do
 
     if (SIZE(fraclengths).eq.1) e_stop = .true.
     if (ALL(fraclengths(:).eq.c0)) e_stop = .true.

     frachistogram(:)=c0


     if (.not. e_stop) then

                ! convert from diameter to radii
                fraclengths(:)=fraclengths(:)/c2

                ! bin into FS cats
                ! highest cat cannot be fractured into
                do j=1,size(fraclengths)
                        do k=1,nfsd-1
                                if ((fraclengths(j).ge.floe_rad_l(k)).and.(fraclengths(j).lt.floe_rad_l(k+1))) then
                                        frachistogram(k)=frachistogram(k)+1
                                end if
                        end do

                        if (fraclengths(j).lt.(floe_rad_l(1)-puny)) stop 'fractures too small'
                
                end do
     end if
 
     deallocate(fraclengths)

     do k=1,nfsd
         frac_local(k)=floe_rad_c(k)*frachistogram(k)
     end do

     ! normalize to one (not D/2, as no longer multiply by fraction of domain reached by waves)              
     if (SUM(frac_local).ne.c0) frac_local(:) = frac_local(:) / SUM(frac_local(:))


     end subroutine wave_frac

!===========================================================================
!  Given the (attenuated) sea surface height, find the strain across triplets
!  of max, min, max or min, max, min (local extrema within 10m).
!  If this strain is greater than the  critical strain, ice can fracture
!  and new floes are formed with sizes equal to the distances between
!  extrema
!
        subroutine get_fraclengths(X, eta, fraclengths, hbar, e_stop)

        use ice_communicate, only: my_task

        real (kind=dbl_kind) :: &
                hbar

        real (kind=dbl_kind), intent(in), dimension (nx) :: &
                X, eta
 
        real (kind=dbl_kind), intent(inout), dimension (:), allocatable :: &
                fraclengths
 
        logical (kind=log_kind), intent(inout) :: &
                e_stop          ! if true, stop and return zero omega and fsdformed


        ! local
        integer (kind=int_kind) :: &
                spcing, &       ! distance in dx over which to search for extrema on each side of point
                j, k, &         ! indices to iterate over domain
                first, last, &  ! indices over which to search for extrema
                j_neg, &        ! nearest extrema backwards
                j_pos, &        ! nearest extrema forwards
                n_above         ! number of points where strain is above
                                ! critical strain
        

         real (kind=dbl_kind), dimension(nx) :: &
                strain, &        ! the strain between triplets of extrema
                frac_size_one, & !
                frac_size_two

        logical (kind=log_kind), dimension(nx) :: &
                is_max, is_min, &       ! arrays to hold whether each point is a local max or min
                is_extremum, &          ! or extremum
                is_triplet              ! or triplet of extrema

         real (kind=dbl_kind) :: &
                delta, &        ! difference in x between current and prev extrema
                delta_pos       ! difference in x between next and current extrema

       integer (kind=int_kind), dimension(1) :: &
        maxj, minj  ! indices of local max and min


       ! -------equivalent of peakfinder2
       !given eta and spcing, compute extremelocs in ascending order
       spcing=nint(threshold/dx)

       is_max = .false.
       is_min = .false.
       is_extremum = .false.
       is_triplet = .false.
       strain = c0
       frac_size_one = c0
       frac_size_two = c0
       j_neg = 0
       j_pos = 0      
 
       ! search for local max and min within spacing of
       ! 10m on either side of each point

       do j = 1, nx

                first = MAX(1,j-spcing)
                last = MIN(nx,j+spcing)


                maxj = MAXLOC(eta(first:last))
                minj = MINLOC(eta(first:last))

                if (COUNT(eta(first:last).eq.MAXVAL(eta(first:last))).gt.1) &
                        stop 'more than one max'
                if (COUNT(eta(first:last).eq.MINVAL(eta(first:last))).gt.1) &
                        stop 'more than one min'

                if (maxj(1)+first-1.eq.j) is_max(j) = .true.
                if (minj(1)+first-1.eq.j) is_min(j) = .true.

                if (is_min(j).and.is_max(j)) then
                         print *, 'X ',X
                         print *, 'eta ',eta
                         print *, 'frst last' ,first, last
                         print *, 'maxj, minj ',maxj,minj
                         stop &
                        'error in extrema'
                end if
                if (is_min(j).or.is_max(j)) is_extremum(j) = .true.
        end do

        do j = 2, nx-1
                if (is_extremum(j)) then
                        if (j.eq.2) then
                                if (is_extremum(1)) j_neg = 1
                        else
                            do k = j-1, 1, -1
                                if (is_extremum(k)) then
                                        j_neg = k
                                        EXIT
                                end if
                             end do
                        end if
                       
                        do k = j+1, nx
                                if (is_extremum(k)) then
                                        j_pos = k
                                        EXIT
                                end if
                        end do
                       
                        if ((j_neg.gt.0).and.(j_pos.gt.0)) then 
                            if (is_max(j_neg).and.is_min(j).and.is_max(j_pos)) &
                                is_triplet(j) = .true.
                            if (is_min(j_neg).and.is_max(j).and.is_min(j_pos)) &
                                is_triplet(j) = .true.
                        end if

                        ! calculate strain
                        if (is_triplet(j)) then
                                delta_pos = X(j_pos) - X(j)
                                delta = X(j) - X(j_neg)
                                strain(j) = (hbar/c2) *(eta(j_neg)*delta_pos - &
                                        eta(j)*(delta_pos+delta) + eta(j)*delta ) &
                                        /(delta*delta_pos*(delta+delta_pos))


                                if (strain(j).gt.straincrit) then
                                        frac_size_one(j) = X(j_pos) - X(j)
                                        frac_size_two(j) = X(j) - X(j_neg)
                                end if
                        end if

                end if

        end do

        n_above = COUNT(strain.gt.straincrit)
        if (n_above.gt.0) then
                DEALLOCATE(fraclengths)
                ALLOCATE(fraclengths(2*n_above))
                fraclengths(1:n_above) = PACK(frac_size_one,(frac_size_one.gt.c0))
                fraclengths(n_above+1:2*n_above) = &
                        PACK(frac_size_two,(frac_size_two.gt.c0))

                e_stop = .false.
        else
                e_stop = .true.

        end if

        end subroutine get_fraclengths




 
!=======================================================================
     
      end module ice_wavefracspec

!=======================================================================


