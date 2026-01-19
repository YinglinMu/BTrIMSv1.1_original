This repository contains the Lagrangian moisture source identification model BTrIMSv1.1 
(Back-Trajectory for the Identification of Moisture Sources), which is designed to identify 
the moisture sources of precipitation.

The original BTrIMS model was developed by Prof. Jason Evans  
(https://github.com/jpevans/BTrIMS), and was further extended by Dr. Chiara Holgate  
(https://github.com/chiaraholgate/QIBT_shared). The present version has been further modified 
and developed by the author (https://github.com/YinglinMu/BTrIMSv1.1_original).

Different branches of this repository correspond to different assumptions or configuration 
choices used in the model.

All branches share the same core components, including:
- `global_data.f90`, which defines global parameters and variables;
- `util.f90`, which contains utility subroutines and functions independent of specific assumptions;
- `bt_subs.f90` and `back_traj.f90` (the main module), which include subroutines implementing 
  different physical assumptions.

The files `bt_subs.f90` and `back_traj.f90` therefore differ among branches, while the other 
components remain identical.

The `main` branch serves as the baseline configuration. The differences among branches are 
summarized as follows:

- `main`, `isentropic`, and `equivalent_potential_temp` represent different assumptions for 
  vertical movement of air parcels;
- `main`, `cloud_hydrometeors`, and `water_vapor` represent different choices of parcel release height;
- `main` and `bicubic` correspond to different within-grid interpolation methods;
- `main`, `WaterSip`, `splitting`, and `splitting_convection` represent different 
  along-trajectory moisture identification approaches.

Further details on these assumptions and configurations can be found in the following paper:
https://egusphere.copernicus.org/preprints/2025/egusphere-2025-2833/



How to use:
  
1. change setting up in global_data.f90,
2. change run_module.sh according to your environment
3. summit a job running the program
This project is licensed under the [GNU General Public License v3.0](LICENSE)
