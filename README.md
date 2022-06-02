# Rn_mass_balance

Click on the green `Code` button on GitHub.com to clone, fork or download the repository. Use the `Rn_mass_balance.Rproj` file in the local project folder to launch the project in `RStudio` (you need to have `R` and `RStudio` already installed on your computer). **Verify in the top right corner of `RStudio` that you are working within the project.** When launched the first time, the [`renv` package maneger](https://rstudio.github.io/renv/articles/collaborating.html "collaborating with renv") should automatically bootstrap itself, downloading and installing the appropriate version of `renv` into the project library. After this has completed, use `renv::restore()` to restore the project library locally on your machine. The folder structure is the following:

-   `input` - input files (e.g. data in csv format)
-   `output` - output files (e.g. data in csv format)
-   `R` - the `R` code for the project
-   `renv` - contains packages used in the project (should not be edited by a user)

The code in the `analysis.R` script analyzes coastal radon mass balance from coastal radon time series measurements - at least 2 such measurements of the following variables are required: radon in water, radon in air, radon in groundwater, water temperature and salinity, water depth, and wind speed. The code is freely available to download and use. The use of the code assumes that a reference of the source will be included in resulting work (paper, report, presentation). APA and BibTeX citation formats are provided in the right hand column of the GitHub page. The authors made their best efforts but do not guarantee successful application of the code on all platforms and all data formats. The `setup.R` and `util_funs.R` scripts are not intended to be edited by the user. The structure of the `csv` input file is described below.

Input file format – `csv`. The file name is specified near the top of the `analysis.R` script, for example: 
````
# input file name
csv_file_in <- "sgd_ts_data_RADAquaMixDif.csv"
````

The `csv` input file needs to have the following headings:

-   Time and date – acceptable formats: `yyyy/mm/dd hh:mm:ss`, `mm/dd/yyy hh:mm`, `mm/dd/yyy hh:mm:ss`

-   `Rn_air__Bqm3` – radon in air measurement in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in air values throughout the period of the measurements

-   `Rn_gw__Bqm3` – radon in groundwater measurement in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in groundwater values throughout the period of the measurements

-   `Rn_offshore__Bqm3` – radon in offshore ocean measurement in units of Bq/m3, either determined directly or approximated from offshore 226Ra measurements, however this will result in overestimate as it neglects radon evasion, this is used to subtract radon input from offshore by flood tide

-   `Ra226_wat__Bqm3` – 226Ra in coastal ocean in units of Bq/m3, this is used to calculate excess 222Rn; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable 226Ra values throughout the period of the measurements

-   `Ra226_sed__Bqg` – 226Ra in sediments in units of Bq/g, this is used to calculate radon input by diffusion from sediments based on empirical relationship from experimental data by Burnett et al (2003):

-   `f_dif__Bqm2hr` = (495 x `Ra226_sed__Bqg` \* 60 + 18.2) / 24; all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable 226Ra values throughout the period of the measurements

-   `Rn_exch__Bqm3` – radon measured using `RAD-Aqua` in air-water gas exchanger in units of Bq/m3, all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in groundwater values throughout the period of the measurements. None of the values need to be filled in if `Rn_wat__Bqm3` are provided.

-   `Rn_wat__Bqm3` – radon measured in water if available, otherwise `Rn_wat_Bqm3` will be calculated using `Rn_exch_Bqm3` and water `sal_wat` and `temp_wat__C`. None of the values need to be filled in if `Rn_exch__Bqm3` are provided.

-   `temp_wat__C` – coastal water temperature in degrees Celsius

-   `sal_wat` – coastal water salinity

-   `wind__ms` – wind measurements in m/s

-   `depth__m` – coastal water depth measurement in units of m if water column is fully mixed, mixed layer depth if water is stratified. all rows have to be filled even if the same number is applicable for all measurements but this format allows the use of variable Rn in groundwater values throughout the period of the measurements.

-   `f_mix_exp__Bqm2hr` - radon mixing losses may be measured directly using current meters or residence time estimates. If `f_mix_exp__Bqm2hr` are not been provided then losses by mixing are set to equal negative `f_Rn_net__Bqm2hr`, this is a conservative approach providing minimal estimate of mixing loss and more representative radon budgets may result by direct experimental measurements of mixing losses `f_mix_exp__Bqm2hr`. None of the values need to be filled in if experimental estimates are not available.

The `csv` files in the `input` folder are examples of acceptable formats.
