Tropopause folding calculator:
(*original reference to be added*)

original work by B.Skerlak (bojan.skerlak@env.ethz.ch)
extended for netcdf data by A.Pozzer (andrea.pozzer@mpic.de)
Extended for ERA5 data by S. Bartusek (samuel.bartusek@columbia.edu):
 - Edits were made to "3d_labelling_and_fold_id.f90"
 - See Bartusek et al. 2023 Supplementary Material for details
 - "3d_labelling_and_fold_id_bartusek_etal_2023.exe" is the .exe used for Bartusek et al. 2023
 - A new "3d_labelling_and_fold_id.exe" should be generated according to directions below

------------------------------------------------------------------
To use the algorithm:

1) Extract the zip file 
2) Modify the Makefile to the current machine set-up
3) Compile to code with the command "make".
   Other options possible with "make help".
4) Create a namelist file (filename.nml) with the specific of the
   netcdf files to be analyzed (see EMAC.nml for example).
5) run the code with the command "3d_labelling_and_fold_id.exe filename.nml"

-------------------------------------------------------------------
NML format:
the namelist file should contain the following fields:
------------------------------------------------------------------
&CTRL
! input file (path included)
file_input = '/ptmp/andrep/ZANIS/ZANIS__________20000801_0000_fold.nc'
! output file (path included) 
file_output = 'test.nc' 
! name of the longitude dimension in the file
X_name = 'lon'  
! name of the latitude dimension in the file
Y_name = 'lat'  
! name of the level dimension in the file
Z_name = 'lev'  
! name of the time dimension in the file
T_name = 'time'  
! name of the surface pressure
APS_name = 'aps'  
! name of the hybrid coefficient A
HYAM_name = 'hyam'  
! name of the hybrid coefficient B
HYBM_name = 'hybm'  
! name of the specific humidity field
Q_name = 'qm1'  
! name of the potential vorticity field
PV_name = 'PV' 
! name of the potential temperature field
PT_name = 'tpot' 
/
------------------------------------------------------------------
WARNINGS:
1) The code works automatically only if all the needed data (i.e. Q,PV,PT)
   are included in the same input file.
2) Only hybrid pressure coordinates included in the pressure
   calculation.
------------------------------------------------------------------
