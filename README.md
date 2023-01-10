# Rn_mass_balance

Click on the green `Code` button on GitHub.com to download the repository. The `FINIFLUX2.0_2022` folder contains the `Matlab` based program "Finite element method for quantifying groundwater fluxes to streams using Radon", along with installation instructions and a user manual in pdf format. All other files and folders in the repository are associated with `R` scripts that quantify Radon budgets in coastal, estuarine, and lake environments. The rest of this page focuses on these `R` scripts.

The code is freely available to download and use. Please include a reference to the source in the resulting work (paper, report, presentation). APA and BibTeX citation formats are provided in the right hand column of the GitHub page of the repository. The authors made their best efforts but do not guarantee successful application of the code on all platforms and all data formats. 

Use the `Rn_mass_balance.Rproj` file in the local project folder to launch the project in `RStudio` (you need to have `R` and `RStudio` already installed on your computer). **Verify in the top right corner of `RStudio` that you are working within the project.** When the project is launched the first time, the [`renv` package maneger](https://rstudio.github.io/renv/articles/collaborating.html "collaborating with renv") should automatically bootstrap itself, downloading and installing the appropriate version of `renv` into the project library. After this has completed, use `renv::restore()` to restore the project library locally on your machine. Once the project is successfully loaded, open the R script of interest though the File menu. For instance, for the coastal time-series mass balance, select `File -> Open File -> sgd_coastal_ts.R`

The folder structure is the following:

-   `input` - input files (i.e. data in csv format)
-   `output` - output files (i.e. data in csv format)
-   `R` - the `R` code for the project
-   `renv` - contains packages used in the project (should not be edited by a user)

The `R/setup.R` and `R/util_funs.R` scripts are not intended to be edited by the user.

The other scripts in the `R` folder analyze radon mass balance in coastal, estuarine and lake environments as implied by their names. Each step in the scripts is explained thoroughly via comments.

Input file format – `csv`. The file name is specified near the top of the `R` script, for example: 
````
# input file name
csv_file_in <- "sgd_coastal_ts_data.csv"
````

The column headings in the `csv` files list the input variables used in the analysis. The `csv` files in the `input` folder are examples of acceptable formats. Don't change the column headings, and remember that `R` is case sensitive, so `time` $\neq$ `Time`; they are different variables. The descriptions of input variables (some used in multiple `R` scripts) are below:

-   `time` – Time and date, acceptable formats: `yyyy-mm-dd hh:mm:ss`, `mm/dd/yyy hh:mm`, `mm/dd/yyy hh:mm:ss`

-   `Rn_air__Bqm3` – atmospheric radon in air activity in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable atmospheric Rn in air values throughout the period of the measurements

-   `Rn_gw__Bqm3` – radon activity in groundwater in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in groundwater values throughout the period of the measurements. This variable is used for unstratified systems. for stratified estuary option see `Rn_gw_btm__Bqm3` and `Rn_gw_surf__Bqm3`

-   `Rn_gw_btm__Bqm3`  – radon activity in groundwater representing the groundwater end-member discharging below the pycnocline in units of Bq/m3

-   `Rn_gw_surf__Bqm3`  – radon activity in groundwater representing the groundwater end-member discharging above the pycnocline in units of Bq/m3

-   `Rn_offshore__Bqm3` – radon activity in offshore ocean in units of Bq/m3, either determined directly or approximated from offshore dissolved 226Ra measurements, however this latter will result in over-correction as it neglects radon evasion, `Rn_offshore__Bqm3` is used to correct for radon input from offshore brought to the coastline by flood tide

-   `Rn_exch__Bqm3` – radon measured using `RAD-Aqua` in an air-water gas exchanger in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn values throughout the period of the measurements. None of the values need to be filled in if `Rn_wat__Bqm3` are provided.

-   `Rn_wat__Bqm3` – radon activity in water if available, otherwise `Rn_wat_Bqm3` will be calculated using `Rn_exch_Bqm3` and water `sal_wat` and `temp_wat__C`. None of `Rn_wat__Bqm3` values need to be filled in if `Rn_exch__Bqm3` are provided but one of these needs to be provided. The code checks if any `Rn_wat_Bqm3` is provided and if yes, it will only calculate radon mass balance for those rows where it is provided. Do not provide any numbers in `Rn_wat_Bqm3`if `Rn_exch__Bqm3` is to be used.

-   `Rn_wat_ups__Bqm3` – radon activity at the upstream estuarine boundary in units of Bq/m3

-   `Rn_wat_dws__Bqm3` – radon activity at the downstream estuarine boundary in units of Bq/m3

-   `Rn_wat_surf_ups__Bqm3`  – radon activity in surface estuarine water above the pycnocline at the upstream estuarine boundary in units of Bq/m3.

-   `Rn_wat_surf_dws__Bqm3`  – radon activity in surface estuarine water above the pycnocline at the downstream estuarine boundary in units of Bq/m3.

-   `Rn_wat_btm_ups__Bqm3`  – radon activity in bottom estuarine water below the pycnocline at the upstream estuarine boundary in units of Bq/m3.

-   `Rn_wat_btm_dws__Bqm3`  – radon activity in bottom estuarine water below the pycnocline at the downstream estuarine boundary in units of Bq/m3.

-   `Ra226_wat__Bqm3` – dissolved 226Ra in water in units of Bq/m3, this is used to calculate excess 222Rn and ingrowth of 222Rn from `Ra226_wat__Bqm3`dissolved in the water column; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable 226Ra values throughout the period of the measurements

-   `Ra226_wat_surf__Bqm3`  – dissolved 226Ra in surface estuarine water above the pycnocline in units of Bq/m3.

-   `Ra226_wat_btm__Bqm3`  – dissolved 226Ra in bottom estuarine water below the pycnocline in units of Bq/m3.

-   `q_ups__m3d` – river discharge at the upstream estuarine boundary in units of m3/d

-   `q_dws__m3d` – river discharge at the downstream estuarine boundary in units of m3/d

-   `layerID` – stratified lakes will have 3 layers, use `epi` for epilimnion, `meta` for metalimnion, and 'hypo' for hypolimnion to identify variables belonging to each layer

-   `depth__m` – water depth measurement in units of `m` if water column is fully mixed, water layer depth if coastal model brackish water surface plume is considered; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable values throughout the period of the measurements.

-   `d_box__m` – box depth in units of m used in cases of stratified water layers

-   `a_box__m2` – box area in units of m2 used in cases of stratified water layers

-   `temp_wat__C` – water temperature in degrees Celsius

-   `sal_wat` – water salinity

-   `sal_surf` – water salinity above the estuarine pycnocline

-   `sal_btm` – water salinity below the estuarine pycnocline

-   `sal_wat_ups` – water salinity at the upstream estuarine boundary

-   `sal_wat_dws` – water salinity at the downstream estuarine boundary

-   `sal_wat_surf_ups` – water salinity above the estuarine pycnocline at the upstream estuarine boundary

-   `sal_wat_surf_dws` – water salinity above the estuarine pycnocline at the downstream estuarine boundary

-   `sal_wat_btm_ups` – water salinity below the estuarine pycnocline at the upstream estuarine boundary

-   `sal_wat_btm_dws` – water salinity below the estuarine pycnocline at the downstream estuarine boundary

-   `wind__ms` – wind measurements in m/s

-   `wat_current__cms`currents measured in estuary in units of cm/s

-   `f_mix_exp__Bqm2hr` - radon mixing losses in the coastal model may be measured directly using current meters or residence time estimates. If `f_mix_exp__Bqm2hr` are not provided then losses by mixing in teh coastal ocean are set to equal negative `f_Rn_net__Bqm2hr`, this is a conservative approach providing minimal estimate of mixing loss and more representative radon budgets may result by direct experimental measurements of mixing losses provided as `f_mix_exp__Bqm2hr`. No values should be filled in if experimental estimates are not available.

-   `f_dif__Bqm2hr` it is up to the user what method they want to use to derive/estimate Rn diffusion from bottom sediments, this program does not calculate radon fluxes by diffusion from sediments but it uses a user provided value in Bq/m2/hr; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable diffusion values throughout the period of the measurements; set to 0 if unknown or negligible, alse set to 0 if coastal model brackish water surface plume, stratified estuary surface layer or epilimnion and metalimnion in lakes are considered.


