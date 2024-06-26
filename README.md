# autoreg_glommod
A set of functions and their application for the modeling of glomerular mechanics under different autoregulatory conditions

for reference

Richfield, O., R. Cortez, L.G. Navar "" 2024

R scripts are used to generate data which is subsequently plotted as figures by make_figs.m

R scripts:
        tubule_schem.R                  generates data for figure 2
        check_takenaka.R                generates data for figures 4 and 5, computes NSC values in table 4
        autoreg_compliant_glom.R        generates data for figures 6, 7 and 8
        param_opt_pass.R                estimates unknown model parameters
        parms.R                         known model parameters
        libs.R                          R libraries
        funcs_3_glomSS.R                set of functions that allow for parameterization param_opt_pass.R
        funcs_4.R                       set of functions that use the estimated unknown parameters
        train_data.R                    data generated from literature - Takenaka and Bell
        shea_anatomy_pressure.R         anatomical data used to construct the model glomerulus
        prep_anat.R                     generates anatomical variables used in running the glomerulus model
        
Data included as prerecquisite for use in the R scripts:

      surr_glom_*.RDS are lookup tables generated by surr_glom_1.R to reduce computational cost 
