# Rn_mass_balance

Click on the green `Code` button on GitHub.com to download the repository. The `FINIFLUX2.0_2022` folder contains the `Matlab` based program "Finite element method for quantifying groundwater fluxes to streams using Radon", along with installation instructions and a user manual in pdf format. All other files and folders in the repository are associated with `R` scripts that quantify Radon budgets in coastal, estuarine, and lake environments. The rest of this page focuses on these `R` scripts.

The code is freely available to download and use. Please include a reference to the source in the resulting work (paper, report, presentation). APA and BibTeX citation formats are provided in the right hand column of the GitHub page of the repository. The authors made their best efforts but do not guarantee successful application of the code on all platforms and all data formats. 

Use the `Rn_mass_balance.Rproj` file in the local project folder to launch the project in `RStudio` (you need to have `R` and `RStudio` already installed on your computer). **Verify in the top right corner of `RStudio` that you are working within the project.** When the project is launched the first time, the [`renv` package maneger](https://rstudio.github.io/renv/articles/collaborating.html "collaborating with renv") should automatically bootstrap itself, downloading and installing the appropriate version of `renv` into the project library. After this has completed, use `renv::restore()` to restore the project library locally on your machine. 

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

The column headings in the `csv` files list the input variables used in the analysis; their descriptions are below:

-   Time and date – acceptable formats: `yyyy/mm/dd hh:mm:ss`, `mm/dd/yyy hh:mm`, `mm/dd/yyy hh:mm:ss`

-   `Rn_air__Bqm3` – atmospheric radon in air activity in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable atmospheric Rn in air values throughout the period of the measurements

-   `Rn_gw__Bqm3` – radon activity in groundwater in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in groundwater values throughout the period of the measurements

-   `Rn_offshore__Bqm3` – radon activity in offshore ocean in units of Bq/m3, either determined directly or approximated from offshore dissolved 226Ra measurements, however this latter will result in over-correction as it neglects radon evasion, `Rn_offshore__Bqm3` is used to correct for radon input from offshore brought to the coastline by flood tide

-   `Rn_dws_Bqm3` – radon activity at the downstream boundary in units of Bq/m3

-   `q_dws__m3d` – river discharge at the downstream boundary in units of m3/d

-   `Rn_ups__Bqm3` – radon activity at the upstream boundary in units of Bq/m3

-   `q_ups__m3d` – river discharge at the upstream boundary in units of m3/d

-   `v_box__m3` – box volume in units of m3

-   `a_box__m2` – box area in units of m2

-   `Ra226_wat__Bqm3` – dissolved 226Ra in water in units of Bq/m3, this is used to calculate excess 222Rn; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable 226Ra values throughout the period of the measurements

-   `Ra226_sed__Bqg` – 226Ra in sediments in units of Bq/g, this is used to calculate radon input by diffusion from sediments based on empirical relationship from experimental data by Burnett et al (2003):

-   `f_dif__Bqm2hr` = (495 x `Ra226_sed__Bqg` \* 60 + 18.2) / 24; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable 226Ra values throughout the period of the measurements

-   `Rn_exch__Bqm3` – radon measured using `RAD-Aqua` in an air-water gas exchanger in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn values throughout the period of the measurements. None of the values need to be filled in if `Rn_wat__Bqm3` are provided.

-   `Rn_wat__Bqm3` – radon activity in water if available, otherwise `Rn_wat_Bqm3` will be calculated using `Rn_exch_Bqm3` and water `sal_wat` and `temp_wat__C`. None of `Rn_wat__Bqm3` values need to be filled in if `Rn_exch__Bqm3` are provided but one of these needs to be provided. The code checks if any `Rn_wat_Bqm3` is provided and if yes, it will only calculate radon mass balance for those rows where it is provided. Do not provide any numbers in `Rn_wat_Bqm3`if `Rn_exch__Bqm3` is to be used.

-   `Rn_wat_bot__Bqm3` – radon activity in bottom estuarine water below the pycnocline in units of Bq/m3
-   `Rn_wat_surf__Bqm3`  – radon activity in surface estuarine water above the pycnocline in units of Bq/m3

-   `Rn_gw_surf__Bqm3`  – radon activity in groundwater representing the groundwater end-member discharging above the pycnocline in units of Bq/m3

-   `Rn_gw_bot__Bqm3`  – radon activity in groundwater representing the groundwater end-member discharging below the pycnocline in units of Bq/m3

-   `temp_wat__C` – water temperature in degrees Celsius

-   `sal_wat` – coastal water salinity

-   `wind__ms` – wind measurements in m/s

-   `depth__m` – water depth measurement in units of `m` if water column is fully mixed, mixed layer depth if water is stratified. all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable values throughout the period of the measurements.

-   `f_mix_exp__Bqm2hr` - radon mixing losses may be measured directly using current meters or residence time estimates. If `f_mix_exp__Bqm2hr` are not provided then losses by mixing are set to equal negative `f_Rn_net__Bqm2hr`, this is a conservative approach providing minimal estimate of mixing loss and more representative radon budgets may result by direct experimental measurements of mixing losses provided as `f_mix_exp__Bqm2hr`. None of the values need to be filled in if experimental estimates are not available.

- `wat_current__cms`currents measured in estuary in units of cm/s

The `csv` files in the `input` folder are examples of acceptable formats.

