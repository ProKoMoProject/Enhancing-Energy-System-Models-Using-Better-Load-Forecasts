$ontext
Dispatch model project ProKoMo

Chair of Energy Economics BTU Cottbus-Senftenberg

target:
    price forecast for 24 hours of the next day

method:
    - rolling horizon including a time intervall before and after the forecasted day
    - for storages a fixed (water)price is implemented
    - data is uploaded by annual data sets
    - years, months, days and hours are mapped
    - UTC time is applied for all data
    - data is uploaded only for the respective year
    - the model loops 374 runs to generate price estimators for one year

$offtext


*#############################  DEFAULT OPTIONS  #############################
$eolcom #

$setglobal Startup ""     # if "*" the startup functions excluded, if "" startup functions included
$setglobal Flow   ""      # if "*" the trade excluded, if "" trade included
$setglobal CHP    ""      # if "*" the trade excluded, if "" trade included
$setglobal ConPow ""      # if "*" Control Power excluded, if "" Control Power included


$ifthen "%Startup%" == ""   $setglobal exc_Startup "*"
$else                       $setglobal exc_Startup ""
$endif



*#####################  DIRECTORIRY and FILE MANAGEMENT  #####################

$setglobal YearonFocus "2019"

*Location of input files
$setglobal datadir                data\
$setglobal DataIn_yearly              InputData%YearonFocus%
$setglobal DataIn_general             InputData_allyears

*Location of output files
$setglobal output_dir   output\
$setglobal result       Results_year%YearonFocus%

set
    daily_window  all days of the model horizon /day1*day374/

    t      all hours                       / t1*t8928  /
;

*#############################   DATA LOAD     ###############################

$include 01_declare_parameters.gms

*#############################   REPORTING INPUT ###############################

*execute_unload '%datadir%Input_final.gdx'
*$stop

*#############################   MODEL     #####################################

$include 02_MODEL.gms

*#############################   SOLVING     ###################################

$include 03_loop.gms

*#############################   results     #################################

$include 04_aftermath.gms












