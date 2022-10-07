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


  ! initialize forpy
  ierror = forpy_initialize()
  ! add path of python module file to paths
  !@thomas TODO relativer pfad!?! Geht nicht wenn jmd die xtb binary verwendet -> error handling
  ! ist halt nicht gegeben das der user das file überhaupt hat. Also unnötig das 
  ! irgendwie weiter zu suchen, ggf kann man gucken obs ne möglichkeit gibt das 
  ! zu "installieren" oder sowas
  ierror = get_sys_path(paths)
  ierror = paths%append(".") ! the module ml_mod.py should be in cwd

  !(1)! convert 3D to SMILE 
  write(*,*) 'Coordinate file: ',fname
  !@thomas added -xh to include H in output smile, xO for atom order
  cmd="obabel "//fname//" -O out.smi -b -xh -xO > atomOrder.out"   !write(*,*) 'cmd |>',cmd,'<|'
  call execute_command_line(cmd)
  ! get atom order from atomOrder.out file
  allocate(ffml%ref2o(mol%n),ffml%o2ref(mol%n), source=0)
  call atom_order_from_file(env, mol%n, ffml%ref2o, ffml%o2ref)

  !(2)! remove cis/trans info from SMILE (represented as / \ )
  !@thomas TODO not good to use external script here?!
  cmd="~/bin/rmCisTrans.sh"
  call execute_command_line(cmd)

  !(3)! convert SMILE to 2D   geo_2D.sdf
  ! file names are hardcoded from here on
  write(*,*) 'Convert SMILE to 2D'
  cmd="obabel out.smi -O geo_2D.sdf --gen2d -h"  ! using obabel SMILE 
!  cmd="obabel openff.smi -O geo_2D.sdf --gen2d -h"  ! using openFF SMILE
  call execute_command_line(cmd)

  !(4)! load geo_2D.sdf geometry
  ! need env, ftype for sdf is 6 (see comment below), 
  !sdf_ftype = getFileType('geo_2D.sdf')
  !write(*,*) 'sdf_ftype =',sdf_ftype  ! gives sdf_ftype = 6
  sdf_ftype = 6
  call init(envtmp)
  call open_file(ich, 'geo_2D.sdf', 'r')
  call readMolecule(envtmp, moltmp, ich, sdf_ftype)
  call close_file(ich)

  !(5)! run structure converter (2D to 3D) with deterministic (fixed) MD
  ! retrieved Data: maxscciter=300.0_wp, maxoptcycle=0,
  ! ml_struc_converter also optimizes structure:
  !(6)! optimize reference structure
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! IMPORTANT NOTE !!! IF OPTIMIZATION FAILS: gfnff_topo file might be wrong !
  !  maybe reference topo is too different from original topo
  !@thomas important !!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !@thomas change xyz sorting back to orig sorting so i can use orig topo
  xyztmp=0.0_wp
  do i=1, mol%n
    xyztmp(:,i)=moltmp%xyz(:,ffml%ref2o(i))
  enddo
  moltmp=mol
  moltmp%xyz=xyztmp
  !@thomas
  egap=0.0_wp
  etemp=0.0_wp
  scciter=250
  maxoptcycle=0
  etot=0.0_wp
  stmp=0.0_wp
  allocate(gtmp(3,moltmp%n), source=0.0_wp)
  call ml_struc_convert(envtmp, .false., moltmp, chktmp, egap, &
          etemp, scciter, maxoptcycle, etot, gtmp, stmp, ffml, topo, refTopo)

  !(7)! convert to original file format  ! use mctc file reader and writer
  call open_file(ich, 'ref_struc.xyz', 'w')
  call writeMolecule(moltmp, ich, format=fileType%xyz, energy=etot, &
   & gnorm=NORM2(gtmp))                                                   
  call close_file(ich)
write(*,'(a,f20.15,a)') ' Found reference structure with energy of ',etot,' Eh'
write(*,*) 'The optimized geometry is written to ref_struc.xyz'


  !@thomas eventually remove this test for production version TODO
  ! testing if mapping is correct; comparing atom type and number of neighbors
  do i=1, mol%n
    !@thomas important: Mapping might be correct, but retrieving values not trivial!
    ! e.g. in topo%nb neighbors are sorted by index -> first neighbor in original might
    !      be the second neighbor in reference topo (due to different index)
    if(mol%at(i).ne.moltmp%at(ffml%o2ref(i))) then
      call env%warning('Mapping from input to reference structure failed.',source)
      exit
    endif
  enddo
  !@thomas del sec end

  !@thomas testing interaction with python functions !@thomas_mark01
  call send_input_to_python_ML_receive_eg(env, mol%n, mol, topo, nlist, ffml)
  ! destroy forpy path object
  call paths%destroy
  ! finalize forpy 
  call forpy_finalize

end subroutine calc_ML_correction


subroutine ml_struc_convert( &
         & env,restart,mol,chk,egap,et,maxiter,maxcycle,&
         & etot,g,sigma, ffml, topo, refTopo)
  use xtb_mctc_accuracy, only : wp
  use xtb_gfnff_param
  use xtb_gfnff_setup
  use xtb_disp_dftd3param
  use xtb_type_environment
  use xtb_type_molecule
  use xtb_type_restart
  use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
  use xtb_type_data
  use xtb_restart
  use xtb_setmod
  use xtb_setparam
  use xtb_dynamic
  use xtb_geoopt
  use xtb_readin, only : xfind
  use xtb_gfnff_setup !@thomas delete? or use for updating topo%nb TODO
  implicit none
! Dummy -----------------------------------------------------------------------
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule),intent(inout)               :: mol
  type(TRestart),intent(inout)                :: chk
  integer,intent(in)                          :: maxiter
  integer,intent(in)                          :: maxcycle
  real(wp),intent(inout)                      :: etot
  real(wp),intent(in)                         :: et
  real(wp),intent(inout)                      :: egap
  real(wp),intent(inout)                      :: g(3,mol%n)
  real(wp),intent(inout)                      :: sigma(3,3)
  type(Tffml), intent(in)        :: ffml
  type(TGFFTopology), intent(in)              :: topo
  type(TGFFTopology), intent(out)             :: refTopo
  logical,intent(in)                          :: restart
  character(len=:),allocatable                :: fnv
! Stack -----------------------------------------------------------------------
  type(TEnvironment)                          :: env2
  type(TGFFCalculator)                        :: calc, calc2, empty_calc
  integer                                     :: ich
  integer                                     :: idum
  integer                                     :: mode_input
  real(wp)                                    :: time_in
  real(wp)                                    :: temp_in
  real(wp)                                    :: step_in
  real(wp)                                    :: dump_in
  real(wp)                                    :: hmass_in
  logical                                     :: exist
  logical                                     :: fail
  logical                                     :: search_done
  integer, allocatable :: opt_in
  character(len=*),parameter                  :: p_fname_param_gfnff = '.param_gfnff.xtb'
! loop geoopt -----------------------------------------------------------------
  integer                                     :: i, j, k, l, nruns
  integer                                     :: num_shift_runs !@thomas changed from 3 !!!!<<<<
  real(wp), allocatable                       :: shift(:)
  real(wp), allocatable                       :: mol_xyz_arr(:,:,:)
  real(wp), allocatable                       :: etot_arr(:)
  real(wp), allocatable                       :: g_arr(:)
  real(wp)                                    :: sign_threshold
  type(TMolecule)                             :: mol_shifted
  character(len=*), parameter :: source = "gfnff_ffml"
  real(wp), parameter :: init_etot = 999.0_wp
  integer :: o2r(mol%n),r2o(mol%n) ! atom index mapping from orig to ref and vice versa
  character(len=1024) :: emsg
  !integer :: nseed
  integer :: seed(2), normcn(86)
  real(wp) :: checkRNG
  integer :: vague_nnb(mol%n) !@thomas vague number of neighbors slight overestimation ...
  ! ...   only interested in finding optimized broken bonds
  real(wp)  :: distances(mol%n*(mol%n+1)/2), mchar(mol%n), cov_rad(mol%n*(mol%n+1)/2) !@thomas new 4refstruc
  integer :: nnb(mol%n)
!  real(wp) :: shvec(3)
!  real(wp) :: sign_threshold

  ! number of "random" shifts plus optimization 
  !@thomas TODO critical for most systems 3 is enough
  !    for odd bonding situation it needs about 100-500
  num_shift_runs=32 !@thomas TODO < If I can map sorting in 2d input to the original structure
    ! then I could use the orignial gfnff_topo and retrieve topo%nb/neigh%nb and use it for
    ! setting up the shift below
  allocate(shift(mol%n), source=0.0_wp)
  allocate(mol_xyz_arr(3, mol%n, num_shift_runs),etot_arr(num_shift_runs),g_arr(num_shift_runs), source=0.0_wp)
  !o2r=ffml%o2ref
  !r2o=ffml%ref2o
  etot_arr = init_etot !@thomas debug <<<< maybe leave at 9.01 or go down again
!------------------------------------------------------------------------------
! set up force field
  call struc_convert_header(env%unit)
  if (allocated(set%opt_engine)) then
    opt_in = set%opt_engine
  end if
  set%opt_engine = p_engine_rf
  mode_input = set%mode_extrun
  set%mode_extrun = p_ext_gfnff
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  call newGFFCalculator(env, mol, calc, fnv, restart, gffVersion%harmonic2020)
  ! use those parts of original topology that are really independent of coordinates
  !calc%topo=topo ! use whole topology from original input 
!write(*,*) 'topo from newGFF on 2D-mol'
write(*,*) 'nb20 from orig topo'
write(*,'(23i4)') topo%nb(20,:)
!write(*,*) 'qa'
!write(*,'(12f10.3)') calc%topo%qa
  calc%topo%bpair = topo%bpair
  calc%topo%alphanb = topo%alphanb
  calc%topo%nbond = topo%nbond
  calc%topo%blist = topo%blist
!  calc%topo%nb = topo%nb
!  calc%topo%hyb = topo%hyb
!  calc%topo%nangl = topo%nangl
!  calc%topo%alist = topo%alist
!  calc%topo%ntors = topo%ntors
!  calc%topo%tlist = topo%tlist
!  calc%topo%nbatm = topo%nbatm
!  calc%topo%b3list = topo%b3list
!  calc%topo%bond_hb_nr = topo%bond_hb_nr
!  calc%topo%bond_hb_AH = topo%bond_hb_AH
!  calc%topo%bond_hb_B = topo%bond_hb_B
!  calc%topo%bond_hb_Bn = topo%bond_hb_Bn
!  calc%topo% = topo%


!===============================
! Set Block
  time_in  = set%time_md
  set%time_md  = 5.0_wp              ! short 5 ps MD to move it in 3D
  temp_in  = set%temp_md
  set%temp_md  = 298.0_wp            ! md temperature 298 K
  step_in  = set%tstep_md
  set%tstep_md = 2.5_wp              ! md time step 2.5 fs
  dump_in  = set%dump_md2
  set%dump_md2 = 100.0_wp            ! md dump 100 fs
  hmass_in = set%md_hmass
  set%md_hmass = 4.0_wp              ! md hydrogen mass
  set%nvt_md = .true.                ! md thermostat
  ! set temperature to 0K  for deterministic md
  set%temp_md  = 0.0_wp    ! md temperature 0 K
!===============================
  if (allocated(set%opt_logfile)) then
    fnv = set%opt_logfile
  else
    deallocate(fnv)
  endif
  set%opt_logfile = 'convert.log'

  ! set fix random seed
  seed(1)=2140693379
  seed(2)=499
  call random_seed(put=seed)
  ! check if random seed worked (also in the future if on other PCs still same RNG)
  call RANDOM_NUMBER(checkRNG)
  if(checkRNG.ne.0.469239734525501_wp) call env%warning("State of random number &
     & generator is not as expected.", source)
  if(checkRNG.ne.0.469239734525501_wp) write(*,*) "State of random number &
     & generator is not as expected."

  !-----------------------------------------------------------------------------
  ! force field geometry optimization on shifted 2D structure
  ! loop runs 3 geoopt with different shifts in the new 3rd coordinate
  ! and then keeps the mol%xyz with lowest etot for md
  search_done=.false.
  nruns=1
  do i=1, num_shift_runs !@thomas changed from i=1
    mol_shifted = mol
    ! create array with deterministic alternating shifts between -1 and 1
    shift=0.0_wp
    call RANDOM_NUMBER(shift)
    shift=shift+0.1 !@thomas 
    do j=1, size(shift) !!
      call RANDOM_NUMBER(sign_threshold)
      if (sign_threshold.lt.0.5_wp) shift(j) = -1.0_wp*shift(j)
    enddo
    mol_shifted%xyz(3,:) = mol_shifted%xyz(3,:) + shift  ! apply shifts
!    ! set temperature to 0K
!    set%temp_md  = 0.0_wp    ! md temperature 0 K
    ! crude optimization of shifted coordinates 
  ! p_olev_crude: ethr=5.d-4  gthr=1.d-2  maxcycle=n  acc=3.00d0    
    call geometry_optimization &
        &     (env,mol_shifted,chk,calc,   &
        &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_lax,.false.,.true.,fail)
!        &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_crude,.false.,.true.,fail)
    mol_xyz_arr(:,:,i) = mol_shifted%xyz  ! store optimized xyz
    etot_arr(i) = etot                    ! store energy etot
    g_arr(i) = NORM2(g)                   ! store gradient norm
    ! check for plausible reference candidate every 4 iterations
    if(MODULO(i,4).eq.0) then
      j=minloc(etot_arr, DIM=1)
      if(etot_arr(j).lt.max(0.1_wp*real(mol%n,wp),0.5_wp)) then
        ! assume that structure is good
        search_done=.true.
      endif
    endif 
    nruns=nruns+1
    if(search_done) exit
  enddo

  if(search_done) then
    write(*,'(a,i4,a)') 'Found candidate for reference structure in ',nruns-1,' searches.'
    write(*,'(a,f20.15,a)') '    Crude energy of candidate:', minval(etot_arr), ' Eh'
  else
    write(*,'(a,i4,a)') 'Could not find plausible candidate for reference structure in ', &
         & nruns-1,' searches.'
    write(*,*) 'Lowest lying crude optimized structure is:'
    write(*,'(3f20.12)') mol_xyz_arr(:,:,minloc(etot_arr, DIM=1))
    !@thomas TODO env%warning does not show up, env%error leads to SIGSEGV
  endif
  mol%xyz(:,:) = mol_xyz_arr(:,:,minloc(etot_arr, DIM=1))  ! keep xyz with lowest etot
  deallocate(mol_xyz_arr)

    if (allocated(fnv)) then
      set%opt_logfile = fnv
    else
      deallocate(set%opt_logfile)
    endif
    write(*,*)
!------------------------------------------------------------------------------
! force field md simulation
  idum = 0
  call md                &
      &   (env,mol,chk,calc, &
      &    egap,set%etemp,maxiter,etot,g,sigma,0,set%temp_md,idum)
!------------------------------------------------------------------------------
! set all back to input
  set%time_md  = time_in
  set%temp_md  = temp_in
  set%tstep_md = step_in
  set%dump_md2 = dump_in
  set%md_hmass = hmass_in
  if (allocated(opt_in)) then
    set%opt_engine = opt_in
  else
    deallocate(set%opt_engine)
  end if
  ! optimize the final geometry
  etot=0.0_wp
  g=0.0_wp
  sigma=0.0_wp
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  call newGFFCalculator(env, mol, calc2, fnv, .false.)
  !@thomas TODO think about p_olev_crude vs _sloppy -> test difference
  !             also might need call gfnff_setup to get complete correct topology of candidate
  !             for final GFN-FF optimization below
write(*,*) 'in converter nb:'     
write(*,'(10i4)') calc2%topo%nb(20,:)
!  calc2%topo%bpair = topo%bpair
!  calc2%topo%alphanb = topo%alphanb
!  calc2%topo%nbond = topo%nbond
!  calc2%topo%blist = topo%blist
  ! GFN-FF optimization of best candidate
  call geometry_optimization &
      &     (env,mol,chk,calc2,   &
      &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_sloppy,.false.,.true.,fail)
  ! save reference topology
  refTopo=calc2%topo

!------------------------------------------------------------------------------
  write(*,*)
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,'(10x,"|           2D => 3D conversion done!             |")')
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,*)
  set%mode_extrun = mode_input
  mol%info%two_dimensional=.false.
  call gfnff_param_dealloc(calc%topo)

end subroutine ml_struc_convert


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
  ierror = ndarray_create(q_arr, ffml%q)
  ierror = kwargs%setitem("ffmlq", q_arr)
  ! add neighbors !@thomas topo is from orig input, refTopo is for reference struc
  ierror = ndarray_create(nb_arr, topo%nb)
  ierror = kwargs%setitem("toponb", nb_arr)
  ! add bpair
  ierror = ndarray_create(bpair_arr, topo%bpair)
  ierror = kwargs%setitem("topobpair", bpair_arr)
  ! add alist
  ierror = ndarray_create(alist_arr, topo%alist)
  ierror = kwargs%setitem("topoalist", alist_arr)
  ! add blist
  ierror = ndarray_create(blist_arr, topo%blist)
  ierror = kwargs%setitem("topoblist", blist_arr)
  ! add tlist
  ierror = ndarray_create(tlist_arr, topo%tlist)
  ierror = kwargs%setitem("topotlist", tlist_arr)
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
  ierror = ndarray_create(nlistq_arr, nlist%q)
  ierror = kwargs%setitem("nlistq", nlistq_arr)
  ! add atom wise energy
  ierror = ndarray_create(eatom_arr, ffml%eatoms)
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

end subroutine


subroutine atom_order_from_file(env, n, ref2o, o2ref)
  type(TEnvironment),intent(inout) :: env
  ! number atoms
  integer, intent(in) :: n
  ! atom index mapping: ref:reference  2->to  o:original
  integer, intent(out) :: ref2o(n), o2ref(n)
  ! does file for atom order exist?
  logical :: fexists
  logical :: sortingBroken(n)
  character(len=*), parameter :: source = "gfnff_ffml"
  character(len=13) :: fatmorder ! file with atom order
  character(len=10) :: buffer, buffer2
  character(:), allocatable :: filestr, trimfilestr
  integer :: bufsize, bufsize2, i, totsize, fstart, fend, iatm,idx

  fatmorder = "atomOrder.out"
  ref2o=0
  o2ref=0
  totsize=0
  inquire( file=fatmorder, exist=fexists )

  if(fexists)then
    open(unit=69, file=fatmorder)
    READ(69, "(A)", ADVANCE='NO', SIZE=bufsize, EOR=10, END=20) buffer 
    do while (.true.)
10    totsize = totsize + 10 ! mark for reading file; read chunks of 10 characters
      READ(69, "(A)", ADVANCE='NO', SIZE=bufsize, EOR=10, END=20) buffer 
    enddo
20  close(69)  ! mark for end of file
    allocate(character(len=totsize) :: filestr)
    filestr=''

    ! Now actually read the lines into string
    open(unit=69, file=fatmorder)
    READ(69, "(A)", ADVANCE='NO', SIZE=bufsize2, EOR=30, END=40) buffer2
    do while (.true.)
!30    filestr = trim(filestr)// buffer2
30    filestr = filestr// buffer2
      READ(69, "(A)", ADVANCE='NO', SIZE=bufsize2, EOR=30, END=40) buffer2 
    enddo
40  close(69)
    fstart = index(filestr, '<')+2
    fend = index(filestr, '>')-1
    allocate(character(totsize) :: trimfilestr)
    trimfilestr = filestr(fstart:fend) 
    filestr = ''
    idx=0
    sortingBroken=.true.
    do i=1, 3*n
      fend = index(trimfilestr,' ')
      if(fend.ne.1) then
        filestr = trimfilestr(1:fend)
        trimfilestr=trimfilestr(fend:)
        read(filestr,'(i)') iatm
        if(iatm.ne.0) then
          idx=idx+1
          ref2o(idx) = iatm
          o2ref(iatm) = idx
          sortingBroken(idx)=.false.
        else
          exit
        endif
      else  
        trimfilestr = trimfilestr(2:)
      endif
    enddo
  else ! fexists is false
    call env%error("Can not find atomOrder.out file.", source)
  endif

  if(any(sortingBroken)) call env%error("Could not retrieve atom order from file.", source)
  write(*,*) 'Retrieved atom order:'
  write(*,*) ref2o

end subroutine atom_order_from_file


subroutine smile_handling(env, mol, topo, ffml, ref2o_map, o2ref_map)
  use forpy_mod
  use iso_fortran_env, only: real64
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule), intent(in)     :: mol
  type(TGFFTopology), intent(in)  :: topo
  type(Tffml), intent(in)        :: ffml
  integer, intent(inout) :: ref2o_map(mol%n), o2ref_map(mol%n)  ! index mapping; ref:reference 2->to o:original
  !
  character(len=*), parameter :: source = "gfnff_ffml"
  type(module_py)  :: smile_mod    ! openff module for SMILE generation
  type(object)     :: receive_obj  ! python object received from function call or similar
  type(dict)       :: kwargs, receive_dict
  type(ndarray)    :: nb_arr, molat_arr, tmp1
  type(tuple)      :: args  ! python tuple to feed into ml_function
  type(object)     :: obj1, obj2
  integer*8, pointer ,dimension(:) :: tmp_poin
  integer :: ierror, i
!  character(len=:) :: smile
  write(*,*) 'topo%nb'
  write(*,'(20i4)') topo%nb
  write(*,*)
  write(*,'(10i4)') mol%at

  ierror = import_py(smile_mod, "smile_gen") ! omit the .py
  write(*,*) 'Imported smile_gen.py with ierror=',ierror
  if(ierror.eq.-1)then ! tell user to get module if not found
    call env%error("Can not find smile_mod.py. Please copy the file into your current working &
      &directory. It is available at https://github.com/grimme-lab/xtb/tree/main/src/gfnff",source)
  endif

  ! create tuple (args) containing arguments for ML python function
  ierror = tuple_create(args, 0)  !@thomas empty tuple needed for call_py
  ! create dict containing kwargs
  ierror = dict_create(kwargs)
  ! add atom types
  ierror = ndarray_create(molat_arr, mol%at)
  ierror = kwargs%setitem("molat", molat_arr)
  ! add neighbors
  ierror = ndarray_create(nb_arr, topo%nb)
  ierror = kwargs%setitem("toponb", nb_arr)
  
  ! Call python function that creates SMILE from loaded module
  ierror = call_py(receive_obj, smile_mod, "create_SMILE", args, kwargs)
  write(*,'(a40,i5)') 'Called create_SMILE with ierror=',ierror
  ! unpack received object
  ierror = dict_create(receive_dict)
  ierror = cast(receive_dict, receive_obj)
  write(*,*) 'Received dictionary from create_SMILE with ierror=',ierror
  ! retrieve SMILE
!  ierror = receive_dict%getitem(obj1, "smile")
!  ierror = cast(smile, obj1)
  ! retrieve index mapping (from 3D to smile atom sorting changes)
  ref2o_map=0
  ierror = ndarray_create_empty(tmp1, mol%n, 'int64')  ! 
!  write(*,*) 'ierror a=',ierror
  ierror = receive_dict%getitem(obj2, "idxmap")  ! write  into object format
!  write(*,*) 'ierror b=',ierror
  ierror = cast(tmp1, obj2)                        ! cast object to ndarray
!  write(*,*) 'ierror c=',ierror
  ierror = tmp1%get_data(tmp_poin)                ! get pointer to ndarray data
!  write(*,*) 'ierror d=',ierror
  ref2o_map = tmp_poin                               !
  write(*,'(a)') '--> Fortran mapping ref to orig <--'
  write(*,'(10i5)') ref2o_map
  write(*,'(a)') '-- -- --  --  -- -- --'

  do i=1, mol%n
    !
    o2ref_map(ref2o_map(i))=i
  enddo
  write(*,'(a)') '--> Fortran mapping orig to ref <--'
  write(*,'(10i5)') o2ref_map
  write(*,'(a)') '-- -- --  --  -- -- --'

  ! destroy module object
  call smile_mod%destroy
  call molat_arr%destroy
  call nb_arr%destroy
  call tmp1%destroy
  call obj2%destroy

end subroutine smile_handling

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

!@thomas 4refstruc
subroutine get_cov_rad(fq,n,at,cn,qa,metal,rab)
  implicit none
  real(wp), intent(in) :: fq
  integer,  intent(in) :: n                 ! number of atoms
  integer,  intent(in) :: at(n)             ! ordinal numbers
  real(wp), intent(in) ::  cn(n)             ! normcn
  real(wp), intent(in) :: qa(n)
  integer, intent(in) :: metal(86)
  real(wp), intent(inout) :: rab(n*(n+1)/2)     ! output bond lengths estimates

  real(wp) :: fat(86), fac

      real(wp) :: f1,f2

      integer m,i,j,k,ii,jj,ati,atj,ir,jr
      INTEGER iTabRow6,lin

      real*8 ra,rb,k1,k2,den,ff,p(6,2)
      real(wp), parameter :: en(86) = (/&
      2.30085633, 2.78445145, 1.52956084, 1.51714704, 2.20568300,&
      2.49640820, 2.81007174, 4.51078438, 4.67476223, 3.29383610,&
      2.84505365, 2.20047950, 2.31739628, 2.03636974, 1.97558064,&
      2.13446570, 2.91638164, 1.54098156, 2.91656301, 2.26312147,&
      2.25621439, 1.32628677, 2.27050569, 1.86790977, 2.44759456,&
      2.49480042, 2.91545568, 3.25897750, 2.68723778, 1.86132251,&
      2.01200832, 1.97030722, 1.95495427, 2.68920990, 2.84503857,&
      2.61591858, 2.64188286, 2.28442252, 1.33011187, 1.19809388,&
      1.89181390, 2.40186898, 1.89282464, 3.09963488, 2.50677823,&
      2.61196704, 2.09943450, 2.66930105, 1.78349472, 2.09634533,&
      2.00028974, 1.99869908, 2.59072029, 2.54497829, 2.52387890,&
      2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000,&
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,&
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,&
      2.00000000, 2.30089349, 1.75039077, 1.51785130, 2.62972945,&
      2.75372921, 2.62540906, 2.55860939, 3.32492356, 2.65140898,&
      1.52014458, 2.54984804, 1.72021963, 2.69303422, 1.81031095,&
      2.34224386/)
      real(wp), parameter :: r0(86) = (/&
      0.55682207, 0.80966997, 2.49092101, 1.91705642, 1.35974851,&
      0.98310699, 0.98423007, 0.76716063, 1.06139799, 1.17736822,&
      2.85570926, 2.56149012, 2.31673425, 2.03181740, 1.82568535,&
      1.73685958, 1.97498207, 2.00136196, 3.58772537, 2.68096221,&
      2.23355957, 2.33135502, 2.15870365, 2.10522128, 2.16376162,&
      2.10804037, 1.96460045, 2.00476257, 2.22628712, 2.43846700,&
      2.39408483, 2.24245792, 2.05751204, 2.15427677, 2.27191920,&
      2.19722638, 3.80910350, 3.26020971, 2.99716916, 2.71707818,&
      2.34950167, 2.11644818, 2.47180659, 2.32198800, 2.32809515,&
      2.15244869, 2.55958313, 2.59141300, 2.62030465, 2.39935278,&
      2.56912355, 2.54374096, 2.56914830, 2.53680807, 4.24537037,&
      3.66542289, 3.19903011, 2.80000000, 2.80000000, 2.80000000,&
      2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000,&
      2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000,&
      2.80000000, 2.34880037, 2.37597108, 2.49067697, 2.14100577,&
      2.33473532, 2.19498900, 2.12678348, 2.34895048, 2.33422774,&
      2.86560827, 2.62488837, 2.88376127, 2.75174124, 2.83054552,&
      2.63264944/)
      real(wp), parameter :: cnfak(86) = (/&
      0.17957827, 0.25584045,-0.02485871, 0.00374217, 0.05646607,&
      0.10514203, 0.09753494, 0.30470380, 0.23261783, 0.36752208,&
      0.00131819,-0.00368122,-0.01364510, 0.04265789, 0.07583916,&
      0.08973207,-0.00589677, 0.13689929,-0.01861307, 0.11061699,&
      0.10201137, 0.05426229, 0.06014681, 0.05667719, 0.02992924,&
      0.03764312, 0.06140790, 0.08563465, 0.03707679, 0.03053526,&
     -0.00843454, 0.01887497, 0.06876354, 0.01370795,-0.01129196,&
      0.07226529, 0.01005367, 0.01541506, 0.05301365, 0.07066571,&
      0.07637611, 0.07873977, 0.02997732, 0.04745400, 0.04582912,&
      0.10557321, 0.02167468, 0.05463616, 0.05370913, 0.05985441,&
      0.02793994, 0.02922983, 0.02220438, 0.03340460,-0.04110969,&
     -0.01987240, 0.07260201, 0.07700000, 0.07700000, 0.07700000,&
      0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000,&
      0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000,&
      0.07700000, 0.08379100, 0.07314553, 0.05318438, 0.06799334,&
      0.04671159, 0.06758819, 0.09488437, 0.07556405, 0.13384502,&
      0.03203572, 0.04235009, 0.03153769,-0.00152488, 0.02714675,&
      0.04800662/)

      p(1,1)=    29.84522887
      p(2,1)=    -1.70549806
      p(3,1)=     6.54013762
      p(4,1)=     6.39169003
      p(5,1)=     6.00000000
      p(6,1)=     5.60000000
      p(1,2)=    -8.87843763
      p(2,2)=     2.10878369
      p(3,2)=     0.08009374
      p(4,2)=    -0.85808076
      p(5,2)=    -1.15000000
      p(6,2)=    -1.30000000

! special hacks
  fat=1.0_wp
  fat( 1)=1.02_wp
  fat( 4)=1.03_wp
  fat( 5)=1.02_wp
  fat( 8)=1.02_wp
  fat( 9)=1.05_wp
  fat(10)=1.10_wp
  fat(11)=1.01_wp
  fat(12)=1.02_wp
  fat(15)=0.97_wp
  fat(18)=1.10_wp
  fat(19)=1.02_wp
  fat(20)=1.02_wp
  fat(38)=1.02_wp
  fat(34)=0.99_wp
  fat(50)=1.01_wp
  fat(51)=0.99_wp
  fat(52)=0.95_wp
  fat(53)=0.98_wp
  fat(56)=1.02_wp
  fat(76)=1.02_wp
  fat(82)=1.06_wp
  fat(83)=0.95_wp

  fac=1.05_wp
      do i=1,n
        f1=fq
        ati=at(i)
        if(metal(ati) > 0) f1 = f1 * 2.0d0
         do j=1,i-1
         k=lin(j,i)
         f2=fq
         atj=at(j)
         if(metal(atj) > 0) f2 = f2 * 2.0d0
         ir=itabrow6(ati)
         jr=itabrow6(atj)
         ra=r0(ati)+cnfak(ati)*cn(i)
         rb=r0(atj)+cnfak(atj)*cn(j)
         den=abs(en(ati)-en(atj))
         k1=0.005d0*(p(ir,1)+p(jr,1))
         k2=0.005d0*(p(ir,2)+p(jr,2))
         ff=1.0d0-k1*den-k2*den**2
         rab(k)=(ra+rb)*ff
         rab(k)=rab(k)-qa(i)*f1-qa(j)*f2
!@thomas ABOUT fac:  additional 5% since I am looking for optimized broken bonds
!  hopefully broken bonds and wrong additional bonds due to +5% do not cancel
         rab(k)=rab(k)*fat(ati)*fat(atj)*fac
         enddo
      enddo


end subroutine get_cov_rad 

subroutine get_nnb(mol, dist, rthr, nnb)
  type(TMolecule), intent(in)     :: mol
  real(wp), intent(in)  :: dist(mol%n*(mol%n+1)/2), rthr(86) !@thomas new 4refstruc
  integer, intent(inout) :: nnb(mol%n)

  integer :: i,j,k,kk

  nnb=0
  do i=1, mol%n
    kk=i*(i-1)/2
    do j=1, i-1
      k=kk+j
      if(dist(k).le.rthr(mol%at(i))) then
        nnb(i) = nnb(i) + 1
        nnb(j) = nnb(j) + 1
      endif
    enddo
  enddo


end subroutine get_nnb

end module xtb_gfnff_ffml
