! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Module  for FF calculation with ML correction
module xtb_gfnff_ffml
  use xtb_mctc_accuracy, only : wp
  use xtb_gfnff_topology, only : TGFFTopology, Tffml 
  use xtb_type_molecule                        
  use xtb_type_environment
  use xtb_type_restart, only : TRestart
  use xtb_io_reader, only : readMolecule
  use xtb_io_writer, only : writeMolecule
  use mctc_io_filetype, only : filetype, getFileType => get_filetype
  use forpy_mod
  use xtb_gfnff_neighbourlist, only: TGFFNeighbourList

  implicit none
  private

  public :: calc_ML_correction

  !@thomas_ffml
  ! the Tffml type is now in topology due to circular dependencies problem

contains

! the actual routine for the ML correction
  !@thomas check if all input is needed
subroutine calc_ML_correction(env,ffml,fname,topo,nlist, mol)
  use forpy_mod
  type(TEnvironment),intent(inout)            :: env
  type(Tffml), intent(inout)        :: ffml
  character(len=*), intent(in)    :: fname
  type(TGFFTopology), intent(in)  :: topo
  type(TMolecule), intent(in)     :: mol
  type(TGFFNeighbourList), intent(in) :: nlist

  character(len=*), parameter :: source = "gfnff_ffml"
  type(TGFFTopology)              :: refTopo
  character(len=:), allocatable   :: cmd  !@thomas 
  integer :: ich, sdf_ftype, maxoptcycle, scciter, i,j
  type(TEnvironment) :: envtmp
  type(TMolecule)    :: moltmp
  type(TRestart) :: chktmp
  real(wp) :: egap, etemp, etot, stmp(3,3)
  real(wp), allocatable :: gtmp(:,:)
  integer :: ierror
  type(list)       :: paths      ! for adding path to python module
  real(wp) :: xyztmp(3,mol%n)
  real(wp), allocatable :: model_xyz(:,:)
  integer :: ic


  ! initialize forpy
  ierror = forpy_initialize()
  ! add path of python module file to paths
  !@thomas TODO relativer pfad!?! Geht nicht wenn jmd die xtb binary verwendet -> error handling
  ! ist halt nicht gegeben das der user das file überhaupt hat. Also unnötig das 
  ! irgendwie weiter zu suchen, ggf kann man gucken obs ne möglichkeit gibt das 
  ! zu "installieren" oder sowas
  ierror = get_sys_path(paths)
  ierror = paths%append(".") ! the module ml_mod.py should be in cwd

!@thomas delete start
!@thomas novel reference generator
! generate local surrounding for each atom and calculate interactions
!  ic=1
!  if(allocated(model_xyz)) deallocate(model_xyz)
!  call generate_atom_surrounding(ic, mol%n, mol%at, mol%xyz, topo%nb, model_xyz)
!@thomas delete end

  !@thomas testing interaction with python functions !@thomas_mark01
  call send_input_to_python_ML_receive_eg(env, mol%n, mol, topo, nlist, ffml)
  ! destroy forpy path object
  call paths%destroy
  ! finalize forpy 
  call forpy_finalize

end subroutine calc_ML_correction


subroutine send_input_to_python_ML_receive_eg(env,n, mol, topo, nlist, ffml)
  use forpy_mod
  use iso_fortran_env, only: real64
  type(TEnvironment),intent(inout)            :: env
  integer, intent(in)    :: n  ! number atoms
  type(TMolecule), intent(in)     :: mol
  type(TGFFTopology), intent(in)  :: topo
  type(Tffml), intent(in)        :: ffml
  type(TGFFNeighbourList), intent(in) :: nlist
  !real(wp), intent(in)   :: xyz(3,n) ! xyz coordinates
  real(wp)   :: xyz_local(3,n)  ! xyz coordinates copy
  integer, dimension(2) :: shape_grad
  real(wp) :: ml_energy
  real(wp) :: ml_gradient(3,n), ml_gradientT(n,3)
  real(wp) :: calc_grad(3,n) !@thomas delete debug
  integer :: i,j !@thomas delete debug
  real(wp), pointer ,dimension(:,:) :: tmp_point
  character(len=*), parameter :: source = "gfnff_ffml"
  ! forpy types and needed variables
  type(module_py)  :: ml_module    ! module with all the python functions etc for ML part
  type(object)     :: receive_obj  ! python object received from function call or similar
  type(object)     :: obj1, obj2
  type(dict)       :: kwargs, receive_dict
  type(tuple)      :: args  ! python tuple to feed into ml_function
  type(ndarray)    :: xyz_arr, tmp1, eatom_arr, q_arr, nb_arr, bpair_arr, alist_arr
  type(ndarray)    :: angle_arr, eangl_arr, hyb_arr, phi_arr
  type(ndarray)    :: blist_arr, tlist_arr, vtors_arr, vbond_arr, vangl_arr
  type(ndarray)    :: hblist1_arr, hblist2_arr, hblist3_arr, hbe1_arr, hbe2_arr, hbe3_arr, nlistq_arr
  integer          :: ierror  ! return value for forpy methods (error value)

  
  !@thomas_mark01 goto call


  ! import python ML module
  ierror = import_py(ml_module, "ml_mod") ! omit the .py
  write(*,'(a40,i5)') 'Imported ml_module with ierror=',ierror
  if(ierror.eq.-1)then ! tell user to get module if not found
    call env%error("Can not find ml_module.py. Please copy the file into your current working &
      &directory. It is available at https://github.com/grimme-lab/xtb/tree/main/src/gfnff",source)
  endif

  ! create tuple (args) containing arguments for ML python function
  ierror = tuple_create(args, 0)  !@thomas empty tuple needed for call_py
  ! create dict containing kwargs
  ierror = dict_create(kwargs)
  ! add xyz coordinates
  xyz_local=mol%xyz 
  ierror = ndarray_create(xyz_arr, xyz_local)
  ierror = kwargs%setitem("xyz", xyz_arr)
  ! add gfnff charges
  ierror = ndarray_create(q_arr, ffml%q) !@thomas feat_is 
  ierror = kwargs%setitem("ffmlq", q_arr)
  ! add neighbors !@thomas topo is from orig input, refTopo is for reference struc
  ierror = ndarray_create(nb_arr, topo%nb)
  ierror = kwargs%setitem("toponb", nb_arr)
  ! add hyb
  ierror = ndarray_create(hyb_arr, topo%hyb) !@thomas feat_is
  ierror = kwargs%setitem("topohyb", hyb_arr)
  ! add bpair
  ierror = ndarray_create(bpair_arr, topo%bpair)
  ierror = kwargs%setitem("topobpair", bpair_arr)
  ! add alist
  ierror = ndarray_create(alist_arr, topo%alist) !@thomas feat_is
  ierror = kwargs%setitem("topoalist", alist_arr)
  ! add angle ! angle around atom topo%alist(1, nbond) -> 2 and 3 are the neighbors
  ierror = ndarray_create(angle_arr, ffml%angle) !@thomas feat_is
  ierror = kwargs%setitem("ffmlangle", angle_arr)
  ! add energy for angle triplet
  ierror = ndarray_create(eangl_arr, ffml%eangl) !@thomas feat_is
  ierror = kwargs%setitem("ffmleangl", eangl_arr)
  ! add blist
  ierror = ndarray_create(blist_arr, topo%blist) !@thomas feat_is
  ierror = kwargs%setitem("topoblist", blist_arr)
  ! add tlist
  ierror = ndarray_create(tlist_arr, topo%tlist) !@thomas feat_is
  ierror = kwargs%setitem("topotlist", tlist_arr)
  ! add phi_tors ! dihedral angle for tors in tlist
  ierror = ndarray_create(phi_arr, ffml%phi_tors) !@thomas feat_is
  ierror = kwargs%setitem("ffmlphi_tors", phi_arr)
  ! add vtors
  ierror = ndarray_create(vtors_arr, topo%vtors)
  ierror = kwargs%setitem("topovtors", vtors_arr)
  ! add vbond
  ierror = ndarray_create(vbond_arr, topo%vbond)
  ierror = kwargs%setitem("topovbond", vbond_arr)
  ! add vangl
  ierror = ndarray_create(vangl_arr, topo%vangl)
  ierror = kwargs%setitem("topovangl", vangl_arr)
  ! add hblist1
  ierror = ndarray_create(hblist1_arr, nlist%hblist1)
  ierror = kwargs%setitem("nlisthblist1", hblist1_arr)
  ! add hblist2
  ierror = ndarray_create(hblist2_arr, nlist%hblist2)
  ierror = kwargs%setitem("nlisthblist2", hblist2_arr)
  ! add hblist3 << being nxb
  ierror = ndarray_create(hblist3_arr, nlist%hblist3)
  ierror = kwargs%setitem("nlisthblist3", hblist3_arr)
  ! add hbe1
  ierror = ndarray_create(hbe1_arr, nlist%hbe1)
  ierror = kwargs%setitem("nlisthbe1", hbe1_arr)
  ! add hbe2
  ierror = ndarray_create(hbe2_arr, nlist%hbe2)
  ierror = kwargs%setitem("nlisthbe2", hbe2_arr)
  ! add hbe3
  ierror = ndarray_create(hbe3_arr, nlist%hbe3)
  ierror = kwargs%setitem("nlisthbe3", hbe3_arr)
  ! add q < EEQ charges
  ierror = ndarray_create(nlistq_arr, nlist%q) !@thomas feat_is
  ierror = kwargs%setitem("nlistq", nlistq_arr)
  ! add atom wise energy
  ierror = ndarray_create(eatom_arr, ffml%eatoms) !@thomas feat_is
  ierror = kwargs%setitem("eatoms", eatom_arr)



  ! call ML function from the python module
  ierror = call_py(receive_obj, ml_module, "receive_ml_input_send_output", args, kwargs)
  write(*,'(a40,i5)') 'Called receive_send_fct with ierror=',ierror
!if(.false.)then !@thomas debug
  ! unpack received object
  ierror = dict_create(receive_dict)
  ierror = cast(receive_dict, receive_obj)
  ! retrieve energy
  ierror = receive_dict%getitem(obj1, "energy")
  ierror = cast(ml_energy, obj1)
  ! retrieve gradient
  shape_grad(1) = 3  ! define shape 3 dimensions (x y z)
  shape_grad(2) = n  ! n atoms
  ierror = ndarray_create_empty(tmp1, shape_grad)
  ierror = receive_dict%getitem(obj2, "gradient")  ! write gradient into object format
  ierror = cast(tmp1, obj2)                        ! cast object to ndarray
  ierror = tmp1%get_data(tmp_point,'A')            ! get pointer to ndarray data
  ml_gradient = tmp_point                          ! get ml_gradient through pointer

  ! calculate gradient directly in fortran
  do i=1, 3
    do j=1, n
      calc_grad(i,j) = xyz_local(i,j)*ffml%eatoms(j)
    enddo
  enddo

  ! warning for future me !@thomas
  write(*,*) ''
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) '  when using actual ML python module  '
  write(*,*) '    check what format gradient has    '
  write(*,*) '           (3,n) vs (n,3)             '
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) ''
!@thomas test for data transfer !@thomas important remove test for production version
  if(abs(ml_energy-SUM(abs(xyz_local))).gt.1.0d-9) write(*,*) 'Warning: Energy difference. Code: wd420' 
  if(SUM(abs(calc_grad - ml_gradient)).gt.1.0d-9) write(*,*) 'Warning: Gradient difference. Code: wd420' 
  
  call args%destroy
  call receive_obj%destroy
  call obj1%destroy
  call obj2%destroy
  call receive_dict%destroy
  call tmp1%destroy
  call xyz_arr%destroy
  call q_arr%destroy
  call ml_module%destroy

end subroutine send_input_to_python_ML_receive_eg


subroutine get_new_dist(mol, xyz, distances)
  type(TMolecule), intent(in)     :: mol
  real(wp), intent(in) :: xyz(3, mol%n)
  real(wp), intent(inout)  :: distances(mol%n*(mol%n+1)/2)
  real(wp) :: sqrab  (mol%n*(mol%n+1)/2)
  integer  :: i,j,k,kk
 
!  allocate( rab(mol%n*(mol%n+1)/2), source = 0.0d0 )
  distances = 0.0_wp
 ! allocate( sqrab(mol%n*(mol%n+1)/2), source = 0.0d0 )

      do i=1,mol%n
         kk=i*(i-1)/2
         do j=1,i-1
            k=kk+j
            sqrab(k)=(mol%xyz(1,i)-mol%xyz(1,j))**2+(mol%xyz(2,i)-mol%xyz(2,j))**2+(mol%xyz(3,i)-mol%xyz(3,j))**2
            distances(k)  =sqrt(sqrab(k))
         enddo
      enddo

end subroutine get_new_dist


end module xtb_gfnff_ffml
