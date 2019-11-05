!=========================================================================
!
!  This module contains the subroutines required to define
!  a floe size distribution tracer for sea ice
!
!  Theory based on:
!
!    Horvat, C., & Tziperman, E. (2015). A prognostic model of the sea-ice 
!    floe size and thickness distribution. The Cryosphere, 9(6), 2119–2134.
!    doi:10.5194/tc-9-2119-2015
!
!  and implementation described in:
!
!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018). An emergent
!    sea ice floe size distribution in a global coupled ocean--sea ice model. 
!    Journal of Geophysical Research: Oceans, 123(6), 4322–4337. 
!    doi:10.1029/2017JC013692
!
!  with some modifications.
!
!  authors: Lettie Roach, VUW/NIWA
!           C. M. Bitz, UW
!  
!  2016: CMB started
!  2016-8: LR worked on most of it
!
     module ice_fsd

      use ice_domain_size, only: ncat, nfsd, max_blocks, nfreq
      use ice_blocks, only: nx_block, ny_block
      use ice_state, only: nt_fsd 
      use ice_kinds_mod
      use ice_constants

      implicit none

      private
      public :: init_fsd, init_fsd_bounds,     &
          write_restart_fsd, read_restart_fsd, &
          icepack_renormfsd, wave_dep_growth, partition_area, &
          icepack_mergefsd 

      logical (kind=log_kind), public :: & 
         restart_fsd      ! if .true., read fsd tracer restart file

      real(kind=dbl_kind), dimension(nfsd), save, public ::  &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth     ! fsd size bin width in m (radius)

      real (kind=dbl_kind), parameter,public :: &
         floeshape = 0.66_dbl_kind  ! constant from Steele (unitless)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,max_blocks), public, save :: &
        d_afsd_latg, d_afsd_latm, d_afsd_addnew, d_afsd_merge, d_afsd_wave

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat,max_blocks), public, save :: &
        d_amfstd_latg, d_amfstd_latm, d_amfstd_addnew, d_amfstd_merge, d_amfstd_wave

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), public, save :: &
        d_an_latg, d_an_latm, d_an_addnew

     real(kind=dbl_kind), dimension(nfsd) ::  &
         floe_rad_h,    &  ! fsd size higher bound in m (radius)
         floe_area_l,   &  ! fsd area at lower bound (m^2)
         floe_area_h,   &  ! fsd area at higher bound (m^2)
         floe_area_c       ! fsd area at bin centre (m^2)

      integer(kind=int_kind), dimension(nfsd, nfsd) ::  &
         iweld




!=======================================================================

      contains

!=======================================================================

!  Initialize ice fsd bounds (call whether or not restarting)
!  Define the bounds, midpoints and widths of floe size
!  categories in area and radius
!
!  Note that radius widths cannot be larger than twice previous
!  category width or floe merging will not have an effect
!
!  Note also that the bound of the lowest floe size category is used
!  to define the lead region width and the domain spacing for wave breaking
!
!  authors: Lettie Roach, NIWA/VUW and C. M. Bitz, UW
!
       subroutine init_fsd_bounds

        use ice_domain_size, only: nfsd,ncat
        use ice_constants, only: puny, c2
        use ice_communicate, only: my_task, master_task
        use ice_calendar, only: dt

       integer (kind=int_kind) :: k, a, b, c

       real (kind = dbl_kind) :: test

       real (kind = dbl_kind), dimension (nfsd+1) :: &
         area_lims, area_lims_scaled ! local variables
                                              
       real (kind = dbl_kind), dimension(:), allocatable :: &
           lims


        if (nfsd.eq.24) then

            allocate(lims(24+1))

            lims =   (/  6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                         5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                         3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                         9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                         3.35434988e+03,   4.55051413e+03,   6.17323164e+03,   8.37461170e+03, &
                         1.13610059e+04,   1.54123510e+04,   2.09084095e+04,   2.83643675e+04, &
                         3.84791270e+04 /)
        
        else  if (nfsd.eq.16) then

            allocate(lims(16+1))

            lims =   (/  6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                         5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                         3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                         9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03, &
                         3.35434988e+03 /)
        
        else if (nfsd.eq.12) then

            allocate(lims(12+1))

            lims =   (/  6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                         5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, &
                         3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, &
                         9.45812834e+02/)
 
        else if (nfsd.eq.6) then

            allocate(lims(6+1))

            lims =   (/  6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, &
                         5.24122136e+01,   8.78691405e+01,   1.39518470e+02 /)
 
        else

            stop 'floe size categories not defined for this nfsd'

        end if

        floe_rad_l = lims(:nfsd)
        floe_rad_h = lims(2:)
        floe_rad_c = (floe_rad_h+floe_rad_l)/c2

        floe_area_l = 4.*floeshape*floe_rad_l**2.
        floe_area_c = 4.*floeshape*floe_rad_c**2.
        floe_area_h = 4.*floeshape*floe_rad_h**2.

        floe_binwidth=(floe_rad_h-floe_rad_l)

        
        if (my_task == master_task) then
            write(*,*)&
            'floe size bin info: low, high, center, width, area_c'
            write(*,*) floe_rad_l
            write(*,*) floe_rad_h
            write(*,*) floe_rad_c
            write(*,*) floe_binwidth
            write(*,*) floe_area_l
            write(*,*) floe_area_h   
            write(*,*) floe_area_c
        end if

        ! floe size categories that can combine during welding
        iweld(:,:) = - 999
        do a = 1, nfsd
        do b = 1, nfsd
            test = floe_area_c(a) + floe_area_c(b)

            do c = 1, nfsd-1
                if ((test.ge.floe_area_l(c)).and.(test.lt.floe_area_h(c))) &
                    iweld(a,b) = c 
            end do
            if (test.ge.floe_area_l(nfsd)) iweld(a,b) = nfsd
        end do
        end do


      end subroutine init_fsd_bounds

!=======================================================================
!
!  Initialize the FSD 
!
!  When growing from no-ice conditions, initialize to zero.
!  This allows the FSD to emerge, as described in Roach, Horvat et al. (2018)
!
!  Otherwise initalize with a power law, following Perovich
!  & Jones (2014). The basin-wide applicability of such a 
!  prescribed power law has not yet been tested.
!
!  Perovich, D. K., & Jones, K. F. (2014). The seasonal evolution of 
!  sea ice floe size distribution. Journal of Geophysical Research: Oceans,
!  119(12), 8767–8777. doi:10.1002/2014JC010136
!
!  authors: Lettie Roach, NIWA/VUW


      subroutine init_fsd(ice_ic, nx_block, ny_block, iblk, ncat, nfsd, trcrn)
        
        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

   character(len=char_len_long), intent(in) :: &
                ice_ic           ! method of ice cover initialization

        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             iblk , &
             ncat, nfsd

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
              trcrn     ! tracer array

	real (kind=dbl_kind) :: alpha
    
    real (kind=dbl_kind),  dimension (nx_block,ny_block,ncat) ::  totfrac

        integer (kind=int_kind) :: k

	real  (kind=dbl_kind), dimension (nfsd) :: &
             num_fsd,  &  ! number distribution of floes
             frac_fsd     ! fraction distribution of floes

    if (trim(ice_ic) == 'none') then
        do k=1,nfsd
            trcrn(:,:,nt_fsd+k-1,:)=c0
        end do
 
    else            ! Perovich (2014)
                                
        ! fraction of ice in each floe size and thickness category
        ! same for ALL cells (even where no ice) initially
        alpha = 2.1_dbl_kind
        totfrac = c0                                   ! total fraction of floes 
        do k = 1, nfsd
            num_fsd(k) = (2*floe_rad_c(k))**(-alpha-c1) ! number distribution of floes
            trcrn(:,:,nt_fsd+k-1,:) = num_fsd(k)*floe_area_c(k)*floe_binwidth(k) ! fraction distribution of floes
            totfrac = totfrac + trcrn(:,:,nt_fsd+k-1,:)
        enddo
        do k = 1,nfsd                                                                                             
            trcrn(:,:,nt_fsd+k-1,:) = trcrn(:,:,nt_fsd+k-1,:)/totfrac
        end do
                                                                                                                                            endif ! ice_ic
      end subroutine init_fsd

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_fsd()

      use ice_domain_size, only: ncat,nfsd
      use ice_restart, only: write_restart_field
      use ice_state, only: trcrn, nt_fsd
      use ice_fileunits, only: nu_dump_fsd

      ! local variables

      character*2 ck
      logical (kind=log_kind) :: diag
      integer k

      diag = .true.

      !-----------------------------------------------------------------
      do k=1,nfsd
        write(ck,'(i2.2)') k
        call write_restart_field(nu_dump_fsd,0, trcrn(:,:,nt_fsd+k-1,:,:), &
                            'ruf8','fsd'//'_'//ck,ncat,diag)
      enddo

      end subroutine write_restart_fsd

!=======================================================================

! Reads all values needed for an ice fsd restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_fsd()

      use ice_domain_size, only: ncat,nfsd
      use ice_restart, only: read_restart_field
      use ice_state, only: trcrn, nt_fsd
      use ice_fileunits, only: nu_restart_fsd

      ! local variables

      character*2 ck
      logical (kind=log_kind) :: diag
      integer k

      diag = .true.

      do k=1,nfsd
        write(ck,'(i2.2)') k
        call read_restart_field(nu_restart_fsd,0,trcrn(:,:,nt_fsd+k-1,:,:), &
                 'ruf8','fsd'//'_'//ck,ncat,diag, &
                 field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      end subroutine read_restart_fsd

!=======================================================================
!
! Wrapper for renorm_mfstd
!
      subroutine icepack_renormfsd ( nx_block, ny_block, &
                                     aicen,    trcrn     )

        use ice_state, only: nt_fsd
        use ice_domain_size, only: max_ntrcr

        integer(kind=int_kind), intent(in) :: &
             nx_block , ny_block

        real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &  ! needs to be max_ntrcr here
             intent(inout) :: &
             trcrn     ! tracer array

        real(kind=dbl_kind), dimension(nx_block,ny_block,  ncat), &
             intent(in) :: aicen

        ! local variables

        integer (kind=int_kind) :: &
             i,j


        do j = 1,ny_block
        do i = 1,nx_block
        
            call renorm_mfstd(aicen(i,j,:), trcrn(i,j,nt_fsd:nt_fsd+nfsd-1,:))
            
        enddo !i
        enddo !j

      end subroutine icepack_renormfsd


!=======================================================================
!
! Normalize the floe size distribution so it sums to one in cells with ice.
! The FSD is zero is cells with no ice
!
! Includes some sanity checks for negative numbers
!
      subroutine renorm_mfstd(aicen,trcrn)


        real (kind=dbl_kind), dimension (nfsd,ncat), &  ! needs to be max_ntrcr here
             intent(inout) :: &
             trcrn     ! tracer array

        real(kind=dbl_kind), dimension(ncat), &
             intent(in) :: aicen

        ! local variables

        integer (kind=int_kind) :: &
             n  ! ice thickness category index


        do n = 1,ncat

            ! sanity checks
            if (aicen(n).lt.c0) stop 'negative aice'
                
            if (ANY(trcrn(:,n).lt.c0-1000*puny)) then
                print *, 'mFSTD ',trcrn(:,n)
                print *, 'aicen ',aicen(n)
                stop 'negative mFSTD'
            end if

            ! tiny mFSTD set to zero
            WHERE(trcrn(:,n).lt.puny) trcrn(:,n) = c0        
 
            ! mFSTD is zero when there is no ice
            if (aicen(n).le.puny) then 
                trcrn(:,n) = c0
            else
                if (SUM(trcrn(:,n)).lt.puny) then
                    print *, aicen(n)
                    print *, trcrn(:,n)
                    print *, SUM(trcrn(:,n))
                    stop 'mFSTD zero for non-zero aicen'
                end if
                if (ABS(SUM(trcrn(:,n))-c1).gt.c0) then
                    trcrn(:,n) = trcrn(:,n) / SUM(trcrn(:,n))
                end if
            end if

        enddo !n

      end subroutine renorm_mfstd

!=======================================================================
! 
!  Given a wave spectrum, calculate size of new floes based on 
!  tensile failure, following Shen et al. (2001)
!
!  The tensile mode parameter is based on in-situ measurements
!  by Roach, Smith & Dean (2018).
!
!  authors: Lettie Roach, NIWA/VUW
!


      subroutine wave_dep_growth (local_wave_spec, & 
                                       new_size )

       use ice_flux, only: dfreq, freq

      real (kind=dbl_kind), dimension(nfreq), intent(in) :: &
           local_wave_spec ! e(f), dimension set in ice_forcing

      integer (kind=int_kind), intent(out) :: &
           new_size ! index of floe size category in which new floes will growh

      ! local variables
      real (kind=dbl_kind), parameter :: &
          tensile_param = 0.167_dbl_kind

      real (kind=dbl_kind)  :: &
          mo,   &   ! zeroth moment of the spectrum (m)
          h_s,  &   ! significant wave height (m)
          w_a,  &   ! wave amplitude (m)
          f_p,  &   ! peak frequency (s^-1)
          w_l,  &   ! wavelength from peak freqency (m) 
          d_max, &  ! d_max from tensile failure mode
          r_max     ! radius

      integer (kind=int_kind) :: k
      

      ! zeroth moment
      mo = SUM(local_wave_spec*dfreq)

      ! sig wave height and amplitude
      h_s = c4*SQRT(mo)
      w_a = h_s/c2

      ! peak frequency
      f_p = freq(MAXLOC(local_wave_spec, DIM=1))

      ! wavelength from peak freq
      w_l = gravit / (c2*pi*f_p**c2)

      ! tensile failure
      if (w_a.gt.puny) then
          d_max = SQRT(c2*tensile_param*w_l**c2/(pi**c3*w_a*gravit*rhoi))
          r_max = d_max/c2
      else
          r_max = bignum
      end if

      new_size = nfsd
      do k = 1, nfsd - 1
          if (r_max.le.floe_rad_h(k)) then
              new_size = k
              EXIT
          end if
      end do


      end subroutine wave_dep_growth

!=======================================================================
! 
!  Given the joint ice thickness and floe size distribution, calculate
!  the lead region and the total lateral surface area following Horvat
!  and Tziperman (2015)
!
! author: Lettie Roach, NIWA/VUW

      subroutine partition_area (        aice,     &
                               aicen,    vicen,    &
                               trcrn,    lead_area,&
                               latsurf_area )

      
      real (kind=dbl_kind), intent(in) :: &
         aice        ! ice concentration
        
      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen, &      ! fractional area of ice 
         vicen         ! volume per unit area of ice

      real (kind=dbl_kind), dimension(nfsd,ncat), intent(in) :: &
         trcrn       ! tracer array

      real (kind=dbl_kind), intent(out) :: &
          lead_area, &  ! the fractional area of the lead region
          latsurf_area  ! the area covered by lateral surface of floes
      
      ! local variables
      
      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k                  ! floe size index

      real (kind=dbl_kind) :: &
        width_leadreg, &   ! width of lead region
        thickness          ! actual thickness of ice in thickness cat

    ! Initialize
    lead_area=c0
    latsurf_area=c0

    ! Set the width of the lead region to be the smallest
    ! floe size category, as per Horvat & Tziperman (2015)
    width_leadreg=floe_rad_c(1)
      
    if (aice.gt.puny) then

        ! lead area = sum of areas of annuli around all floes
        do n=1,ncat       
            do k=1,nfsd
                lead_area = lead_area + aicen(n) * trcrn(k,n) * &
                          ( c2*width_leadreg/floe_rad_c(k) + &
                            width_leadreg**c2/floe_rad_c(k)**2 )
             enddo !k
        enddo !n
       
        ! cannot be greater than the open water fraction
        lead_area = MIN(lead_area,(c1-aice))
      
        ! sanity checks
        if (lead_area.gt.c1) stop 'lead_area not frac!'
        if (lead_area.ne.lead_area) stop 'lead_a NaN'
        if (lead_area.lt.c0) then
                if (lead_area.lt.(c0-puny)) then
                        stop 'lead_area lt0 in partition_area'
                else
                        lead_area=MAX(lead_area,c0)
                end if
        end if

        ! Total fractional lateral surface area in each grid (per unit ocean area)
        do n = 1, ncat
            thickness = c0

            if (aicen(n).gt.c0) thickness = vicen(n)/aicen(n)

                do k=1,nfsd
                    latsurf_area = latsurf_area + trcrn(k,n) * aicen(n) * & ! FSD
                                   c2 * thickness/floe_rad_c(k)
                end do ! k
        end do ! n

        ! sanity checks
        if (latsurf_area.lt.c0) stop 'negative latsurf_ area'
        if (lead_area.ne.lead_area) stop 'latsurf_ area NaN'

    end if ! aice

              end subroutine partition_area

!=======================================================================
!
! Wrapper for floe_merge_thermo
!
      subroutine icepack_mergefsd ( &
                               iblk, nx_block, ny_block, &
                               icells, indxi, indxj,     & 
                               dt, aicen, frzmlt,        &
                               trcrn,                    &
                               d_afsd_merge, d_amfstd_merge)

      integer (kind=int_kind), intent(in) :: &
         iblk, &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj         ! compressed i/j indices

     real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen   ! concentration of ice
 
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         frzmlt    ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd,ncat), &
         intent(inout) :: &
         trcrn, &
         d_amfstd_merge

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfsd), &
         intent(inout) :: &
         d_afsd_merge

      ! local variables

      integer (kind=int_kind) :: &
        ij, i, j

      
    do ij = 1, icells
        
        i = indxi(ij)
        j = indxj(ij)

        call floe_merge_thermo(dt, aicen(i,j,:),       & ! in
                               frzmlt(i,j),            & ! in
                               trcrn(i,j,:,:),         & ! inout
                               d_afsd_merge(i,j,:),    & ! out
                               d_amfstd_merge(i,j,:,:) ) ! out


    end do!ij 

    
            end subroutine icepack_mergefsd


!=======================================================================
!
!  Floes are perimitted to weld together in freezing conditions, according
!  to their geometric probability of overlap if placed randomly on the 
!  domain. The coagulation equation is solved using the method of Filbet
!  & Laurencot (2004). The rate per unit area c_weld is the total number 
!  of floes that weld with another, per square meter, per unit time, in the 
!  case of a fully covered ice surface (aice=1), equal to twice the reduction
!  in total floe number. See Roach, Smith & Dean (2018).
!
!  Filbet, F., & Laurençot, P. (2004). Numerical simulation of the Smoluchowski 
!  coagulation equation. SIAM Journal on Scientific Computing, 25(6), 2004–2028. 
!  doi:10.1137/S1064827503429132
!
!  authors: Lettie Roach, NIWA/VUW
!

        subroutine floe_merge_thermo (dt, aicen, frzmlt, &
                                      areal_mfstd,       &
                                      d_afsd_merge,    &
                                      d_amfstd_merge       )

                       

     real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (ncat), &
         intent(in) :: &
         aicen   ! concentration of ice
 
      real (kind=dbl_kind), intent(in) :: &
         frzmlt    ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nfsd,ncat), &
         intent(inout) :: &
         areal_mfstd, &
         d_amfstd_merge

      real (kind=dbl_kind), dimension (nfsd), &
         intent(inout) :: &
         d_afsd_merge


      ! local variables

      real (kind=dbl_kind), parameter :: &
         aminweld = p1 ! minimum ice concentration likely to weld

      real (kind=dbl_kind), parameter :: &
         c_weld = 1.0e-8_dbl_kind     
                          ! constant of proportionality for welding
	 	          ! total number of floes that weld with another, per square meter,
			  ! per unit time, in the case of a fully covered ice surface
	 		  ! units m^-2 s^-1, see documentation for details


      integer (kind=int_kind) :: &
        t, &
        n, k, m, &
        kx, ky, kz, a

      real (kind=dbl_kind), dimension(nfsd) :: &
         stability , & ! check for stability
         nfsd_tmp  , & ! number fsd
         amfstd_init , & ! initial values
         amfstd_tmp  , & ! work array
         gain, loss    ! welding tendencies

      real(kind=dbl_kind) :: &
         prefac    , & ! multiples kernel
         kern      , & ! kernel
         subdt     , & ! subcycling time step for stability (s)
         elapsed_t     ! elapsed subcycling time


    stability = c0
    prefac = p5
  
    do n = 1, ncat
                                    
        d_afsd_merge(:) = c0
        d_amfstd_merge(:,n) = c0
      
        ! If there is some ice in the lower (nfsd-1) categories
        ! and there is freezing potential
        if ((frzmlt.gt.puny).and. & ! if freezing potential
            (aicen(n).gt.aminweld).and.  & ! skip low concentrations (unlikely to merge)
            (SUM(areal_mfstd(:nfsd-1,n)).gt.puny)) then ! some ice in lower (nfsd-1) categories

            ! save initial values
            amfstd_init(:) = areal_mfstd(:,n)
            amfstd_tmp = amfstd_init

            ! in case of minor numerical errors
            WHERE(amfstd_tmp.lt.puny) amfstd_tmp = c0
            amfstd_tmp = amfstd_tmp/SUM(amfstd_tmp)

            ! adaptive sub-timestep
            elapsed_t = c0 
            DO WHILE (elapsed_t < dt) 

               ! calculate sub timestep
               nfsd_tmp = amfstd_tmp/floe_area_c
               stability = nfsd_tmp/(c_weld*amfstd_tmp*aicen(n))
               WHERE (stability.lt.puny) stability = bignum
               subdt = MINVAL(stability)
               subdt = MIN(subdt,dt)

               loss(:) = c0
               gain(:) = c0

               do kx=1,nfsd ! consider loss from this category
               do ky=1,nfsd ! consider all interaction partners

                   k = iweld(kx,ky) ! product of kx and ky

                   if(k.gt.kx) then
                   
                       kern = c_weld * floe_area_c(kx) * aicen(n)
                      
                       loss(kx) = loss(kx) + kern*amfstd_tmp(kx)*amfstd_tmp(ky)

                       if (kx.eq.ky) prefac = c1 ! otherwise 0.5

                       gain(k) = gain(k) + prefac*kern*amfstd_tmp(kx)*amfstd_tmp(ky)

                   end if

               end do
               end do

               ! does largest category lose?
               if (loss(nfsd).gt.puny) stop 'weld, largest cat losing'
               if (gain(1).gt.puny) stop 'weld, smallest cat gaining'

               ! update afsd   
               amfstd_tmp(:) = amfstd_tmp(:) + subdt*(gain(:) - loss(:))

               ! in case of minor numerical errors
               WHERE(amfstd_tmp.lt.puny) amfstd_tmp = c0
               amfstd_tmp = amfstd_tmp/SUM(amfstd_tmp)

               ! update time
               elapsed_t = elapsed_t + subdt

               ! stop if all in largest floe size cat
               if (amfstd_tmp(nfsd).gt.(c1-puny)) exit 

            END DO ! time


           
            ! in case of small numerical errors
            areal_mfstd(:,n) = amfstd_tmp/SUM(amfstd_tmp)

            ! sanity checks
            if (ANY(areal_mfstd(:,n).lt.-puny)) stop 'neg, mrg'

            WHERE(areal_mfstd(:,n).lt.c0) areal_mfstd(:,n) = c0
            
            ! diagnostic
            d_amfstd_merge(:,n) = areal_mfstd(:,n) - amfstd_init

        end if ! try to weld

    end do !n
    
    ! diagnostics
    do k = 1, nfsd
        d_afsd_merge(k) = c0
        do n = 1, ncat
            d_afsd_merge(k) = d_afsd_merge(k)  + &
            aicen(n)* d_amfstd_merge(k,n)
        end do ! n
    end do ! k



       end subroutine floe_merge_thermo




!=======================================================================

      end module ice_fsd

!=======================================================================

