! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff

      implicit none

      include "test_suite_extras_def.inc"

! here are the x controls used below

!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

      contains

      include "test_suite_extras.inc"


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_alpha_mlt => alpha_mlt_routine
      end subroutine extras_controls


      subroutine alpha_mlt_routine(id, ierr)
         use chem_def, only: ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, h1
         real(dp) :: alpha_H, alpha_other, H_limit
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         alpha_H = s% x_ctrl(21)
         alpha_other = s% x_ctrl(22)
         H_limit = s% x_ctrl(23)
         h1 = s% net_iso(ih1)
         !write(*,1) 'alpha_H', alpha_H
         !write(*,1) 'alpha_other', alpha_other
         !write(*,1) 'H_limit', H_limit
         !write(*,2) 'h1', h1
         !write(*,2) 's% nz', s% nz
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
            !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k),
         end do
         !stop
      end subroutine alpha_mlt_routine


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (.not. restart) s% ixtra(1) = 0
         call apply_phase_controls(s)

         if (.not. s% x_logical_ctrl(37)) return

      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         use num_lib, only: find0
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt, m
         integer :: k, nz
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         write(*,'(A)')
         select case (s% x_integer_ctrl(5))
         case (6)
            ! put target info in TestHub output
            testhub_extras_names(1) = 'fe_core_mass'; testhub_extras_vals(1) = s% fe_core_mass

            if(s% fe_core_mass < 1d0) then
               write(*,1) "Bad fe_core_mass", s%fe_core_mass
            else
               if(s% fe_core_infall > s% fe_core_infall_limit) then
                  write(*,'(a)') 'all values are within tolerance'
               else
                  write(*,'(a)') "Bad fe core infall"
               end if
            end if
         end select
         call test_suite_after_evolve(s, ierr)
         if (.not. s% x_logical_ctrl(37)) return
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'zbar_div_abar'
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
      end subroutine data_for_extra_profile_columns

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call update_phase(s)
         if (.not. s% x_logical_ctrl(37)) return
         if (extras_finish_step == terminate) &
             s% termination_code = t_extras_finish_step
      end function extras_finish_step

      ! Phases replicate the six test-suite parts in one continuous run:
      !  0 pre-MS (part 1)         -> 1 at age >= 1d3 yr  (part 2: reset E err)
      !  1 to ZAMS (part 2)        -> 2 at X_c <= 0.68    (part 3: reset E err)
      !  2 to He depl. (part 3)    -> 3 at Y_c < 1d-6     (part 4: min_D_mix)
      !  3 to C depl. (part 4)     -> 4 at C_c < 1d-3     (part 5: E limits)
      !  4 to lgTmax 9.6 (part 5)  -> 5 at lgTmax >= 9.6  (part 6: CC controls)
      subroutine update_phase(s)
         type (star_info), pointer :: s
         integer :: p0
         p0 = s% ixtra(1)
         select case (s% ixtra(1))
         case (0)
            if (s% star_age >= 1d3) s% ixtra(1) = 1
         case (1)
            if (s% center_h1 <= 0.68d0) s% ixtra(1) = 2
         case (2)
            if (s% center_h1 < 1d-6 .and. s% center_he4 < 1d-6) s% ixtra(1) = 3
         case (3)
            if (s% center_c12 < 1d-3) s% ixtra(1) = 4
         case (4)
            if (s% log_max_temperature >= 9.6d0) s% ixtra(1) = 5
         end select
         if (s% ixtra(1) /= p0) then
            write(*,'(a,i2,a,i8)') 'phase ->', s% ixtra(1), ' at model', s% model_number
            if (s% ixtra(1) <= 2) s% cumulative_energy_error = 0d0
            call apply_phase_controls(s)
         end if
      end subroutine update_phase


      subroutine apply_phase_controls(s)
         type (star_info), pointer :: s
         if (s% ixtra(1) >= 3) then  ! parts 4-6
            s% set_min_D_mix = .true.
            s% min_D_mix = 1d-2
         end if
         if (s% ixtra(1) >= 4) then  ! parts 5-6
            s% limit_for_rel_error_in_energy_conservation = 1d-3
            s% hard_limit_for_rel_error_in_energy_conservation = 1d-2
         end if
         if (s% ixtra(1) >= 5) then  ! part 6
            s% op_split_burn_min_T = 1d9
            s% ignore_too_large_correction = .true.
            s% delta_lg_XH_cntr_min = 0.5d0
            s% delta_lg_XO_cntr_min = 0.5d0
            s% delta_lg_XC_cntr_min = 0.5d0
            s% delta_lg_XHe_cntr_min = 0.5d0
            s% delta_XSi_cntr_limit = 0.0025d0
            s% calculate_Brunt_B = .true.
            s% max_q_for_conv_timescale = 0.2d0
            s% max_X_for_conv_timescale = 1d-6
            s% photo_interval = 10
         end if
      end subroutine apply_phase_controls

      end module run_star_extras


