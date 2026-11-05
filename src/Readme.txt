This code is for the following paper:

Mengze Wang & Tamer Zaki, 2021 State estimation in turbulent channel flow from limited observations. J. Fluid Mech.

Code prepared by Mengze Wang (mzwang@jhu.edu)

!---------------------------------------------------------------
!                     Content
!---------------------------------------------------------------
1. Code description
2. Compile the code
3. Run the code
4. Parameters used to generate results in the paper


!---------------------------------------------------------------
!                  1. Code description
!---------------------------------------------------------------
00_lib/
      Library of subroutines used in the DNS code
01_Gen_Cart_Grid/
      Generate Cartesian grid for channel flow
      Read: input.file
      Output:  grid.data
02_setup_separate_uc/
      Set up initial conditions for channel flow
      Read: input.setup   grid.data
      Output: *0000000 
03_Forward_subsample/
      Forward code, generate / read sub-sampled velocity data on the fly
      Used for section 3.1-3.4 of the paper
      Cost function: equation (2.1) and (3.7) of the paper
      Read: input_forward   param.inc
      Output: a.out
      NOTE: when "ifric" in input_forward is set to 1, the velocity at first grid from both walls are observed
03_Forward_walldata/
      Forward code, generate /read wall stress and pressure data on the fly
      Used for section 3.5 of the paper
      Read: input_forward   param.inc
      Output: a.out
      NOTE: the wall data are assumed to be fully resolved in this code
04_Adjoint_subsample/
      Adjoint code of 03_Forward_subsample
      Read: input_adjoint   param.inc
      Output: a.out 
04_Adjoint_walldata/
      Adjoint code of 03_Forward_subsample
      Read: input_adjoint   param.inc
      Output: a.out 
05_opt/
      L-BFGS code, calls forward and adjoint code at every iteration to minimize the cost function
      Read: param.inc
      Output: lbfgs.out
06_LSE_walldata/
      Linear stochastic estimation code
      Reference: Encinar & Jimenez, 2019 Logarithmic-layer turbulence: A view from the wall
      Read: grid.data     a number of velocity fields, 01uc00.*, 02vc00.*, 03wc00.*  to compute linear estimator
      Output: first guess for section 3.5 (01uc00.dble.0000000, 02vc00.dble.0000000, 03wc00.dble.0000000)


!---------------------------------------------------------------
!                  2. Compile the code
!---------------------------------------------------------------
00_lib/
      make all
      Outcome: a folder "lib/" under your home directory
01_Gen_Cart_Grid/
      ./make.bat
      Outcome: Gen_grid.data.exe 
02_setup_separate_uc/
      make chan
      Outcome: a.out
03_Forward_subsample/
03_Forward_walldata/
04_Adjoint_subsample/
04_Adjoint_walldata/
      make imdiffz
      Outcome: a.out
05_opt/
      ./make_command
      Outcome: lbfgs.out
06_LSE_walldata/
      No need to compile Matlab files


!---------------------------------------------------------------
!                  3. Run the code
!---------------------------------------------------------------

(0) Set working directory and copy source code
      mkdir my_da_case
      copy source_code under my_da_case/
      cd my_da_case/
(1) Generate grid
      go to 01_Gen_Cart_Grid/
      set input.file
      ./Gen_grid.data.exe

(2) Generate initial condition
      go to 02_setup_separate_uc/
      set input.setup
      cp ../01_Gen_Cart_Grid/grid.data ./
      ./a.out

(3) Run forward code to generate observation
      go to my_da_case/
      mkdir Truth
      cd Truth/
      cp ../source_code/02_setup_separate_uc/grid.data ./
      cp ../source_code/02_setup_separate_uc/*0000000 ./
      cp archive.0000000 archive
      cp ../source_code/03_Forward_??/a.out ./forward.out
      cp ../source_code/03_Forward_??/input_forward ./
      (NOTE: "??" could be "subsample" or "walldata", depending on the types of observations you consider)
      (NOTE: you might want to replace 01uc00.dble.0000000, 02vc00.dble.0000000, 03wc00.dble.0000000 
       by your "true" initial conditions. These files could be generated from simulation of fully-developed turbulent channel)
      vi data
      write 0 into data file
      set "iassim" in input_forward to be 1
      mkdir obs
      mpirun -np xx ./forward.out > history_obs

(4) Set forward code to compute the cost function
      go back to my_da_case/
      mkdir da_run
      cd da_run/
      cp ../Truth/* ./
      set "iassim" in input_forward to be 2
      replace 01uc00.dble.0000000, 02vc00.dble.0000000, 03wc00.dble.0000000 by the first guess
      mkdir for
      set 0 into data file
      cp archive.0000000 archive
      You can try "mpirun -np xx ./forward.out > history_forward" to check if the code is able to run
      The value of cost function is stored in file "cost_function"

(5) Set adjoint code to compute gradient of the cost function
      cp ../source_code/04_Adjoint_??/a.out ./adjoint.out
      cp ../source_code/04_Adjoint_??/input_adjoint ./
      Make sure that input_adjoint is consistent with input_forward (e.g. nsteps, observation setup, iPrj, etc.) 
      vi data
      write 0 into data file
      cp archive.0000000 archive
      You can try "mpirun -np xx ./adjoint.out > history_adjoint" to check if the code is able to run
      The gradient is stored in the following files:
      01ucad.dble.0000000
      02vcad.dble.0000000
      03wcad.dble.0000000

(6) Run L-BFGS to minimize the cost function
      mkdir LBFGS_backup 
      This folder will store the estimates of I.C. and their gradients at every iteration
      vi data_zero
      write 0 into data_zero
      cp ../source_code/05_opt/lbfgs.out ./
      ./lbfgs.out > history_LBFGS
      Now the code should start running
      At every iteration, the lbfgs.out code calls forward.out and adjoint.out 

(7) For wall data: use 06_LSE_walldata/  to generate first guess
      go to my_da_case/
      mkdir data
      cd data/
      use 03_Forward_walldata/ to run a long forward simulation to generate a number of restart files and observations (set iassim = 1)
      go to my_da_case/
      mkdir LSE/
      cp 06_LSE_wall_data/* LSE/
      run LSE_estimator.m (could take a long time)
      run LSE_predictor.m
      The *0000000 files are the first guess estimated using LSE      

!---------------------------------------------------------------
!      4. Parameters used to generate results in the paper
!---------------------------------------------------------------

Flow configuration (Re_tau = 180)
      01_Gen_Cart_Grid/input.file
            Lx           = 12.566370614
            Ly           = 1.0
            ngx          = 385
            ngy          = 257
            imeshY       = 1
            ratio        = 15.0  
      02_setup_separate_uc/input.setup
            nix          = 385
            niy          = 257
            niz          = 321
            xL           = 12.56637061
            yL           = 2.0
            zL           = 6.283185307
            imeshY       = 2
            ratio        = 15.0
            ubulk        = 1.0
            wbulk        = 1.0
            Re           = 2800 
            Other parameters are the same as the default setting
      03_Forward_??/input_forward
            nsteps       = 896   !  assimilation window T+ = 50
            Re           = 2800.
            isymproj     = 1     ! use symmetric projector
            dt           = 0.005
            iPrj         = 0     ! fixed pressure gradient
            GPxi         = -2.011191E-05   ! time average from fixed mdot code
      04_Adjoint_??/input_adjoint
            veloflag     = 2
            RestoreStart = 0
            RestoreEnd   = 896   ! should always equal to nsteps
            ReStoregap   = 1
            AdjointNt    = 1
            Adjointflag  = 2
            Other parameters should be consistent with input_forward
      
Section 3.1
      03_Forward_subsample/input_forward
      04_Adjoint_subsample/input_adjoint
            isubdomain   = 1
            ifric        = 0
            ObsStart     = 0
            ObsEnd       = 896
            ObsNt        = 8
            istart       = 1
            iend         = 384
            iskip        = 8
            jstart       = 1 
            jend         = 256
            jskip        = 8
            kstart       = 1
            kend         = 320
            kskip        = 8  
            sigma        = 0.0

Section 3.2
      03_Forward_subsample/input_forward
      04_Adjoint_subsample/input_adjoint
            Case N1, N2
            sigma        = 0.05, 0.1   
            Other parameters remain the same from Section 3.1

Section 3.3          
      03_Forward_subsample/input_forward
      04_Adjoint_subsample/input_adjoint
            Case RZ1-4
            iskip       = 1
            jskip       = 1
            kskip       = 4, 8, 16, 32
            Case RX1-4
            iskip       = 4, 8, 16, 32
            jskip       = 1
            kskip       = 1
            Case RT1-4
            ObsNt       = 16, 32, 64, 128
            iskip       = 32
            jskip       = 1
            kskip       = 1
            Other parameters remain the same from Section 3.1

Section 3.4
      03_Forward_subsample/input_forward
      04_Adjoint_subsample/input_adjoint
            ifric       = 1
            iskip       = 8
            jstart      = 60
            jend        = 196
            jskip       = 8
            kskip       = 8
            Other parameters remain the same from Section 3.1

Section 3.5
      03_Forward_walldata/input_forward
      04_Adjoint_walldata/input_adjoint
           ifric        = 1
           ipre         = 1
           ObsStart     = 1
           ObsEnd       = 896
           ObsNt        = 1



