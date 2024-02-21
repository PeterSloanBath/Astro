# Astro
Set of jupyter notebooks for astronomy data reduction and photometry primarily for variable stars

file structure is

- \name of target\date\bias
- \name of target\date\dark
- \name of target\date\flat\filter
- \name of target\date\light\filter

Most of these note books need some of the funciton in the My_Functions_v3.py file. See headers.

Run in this order :

- Bias_stats_master
  -Takes all the bias frames and combines into a master bias frame

- Dark_stats_master
  - Takes all the dark frames and combines into a master dark frame in flux (per sec). Has master bias removed.

- Flat_stats_master
  - Takes all the flat frames for a filter and combines into a master flat frame in flux. Has master bias removed and dark counts.
    
- Data_Reduction_v8
  - Reads all raw light frames taken with a filter, reduces as (Light - Mster_Dark*Exposuretime - Master_Bias) / Master_Flat_filter
  - This version has the option of output as singles (not doubles) and also bin the data to produce smaller output files.
  - Also saves an "image" that has a value of 1 if the orignal raw image was not saturated, or zero if it was.

- Find_Stars_aperture_Photometry_v4
  - This will open all the reduced files and find the stars, perform Aperture photometry and then save all that as a .csv
  - Needs the aperture function in the My_Functions_v3.py file

BY HAND
- Look at each reduced file
  - if it no good, dim, blury, bumped, remove it to a sub folder called "NoGood" and so that file will be ignores in the following analysis.
  - I use the https://github.com/siyu6974/QuickLook.Plugin.FitsViewer/blob/master/README.md plugin to easily do this.

If you have lots of new data then you may want to...
- Make_Star_Master_v4
  - Take all the images in a filter band and make them into a composite image, send to astrometry.net and then SIMBAD so align and solve.
  - Can be memory hungry as it opens all the images at once.
  - Saves a weight-mean image and a (non-weighted) median image. Excludes pizels that are saturated.
  - Produces an compsite saturated "image"

- Find_Normalisation stars_v1
  - Cycles through all the found stars in the master image, look at them to see if they have neighbours, or are close to the edge etc. and asks the user if they are worthy of being a normlisation star or not.
  - Run for each filter

- Star_align_v4
  - Reads in the various output and make a new registered csv with the id_ref for each star common to all the images.
  - Run each filter separetly


- Normalise_form_Star_align
 - Reads in the registered csv's and find a common set of local-standard and works out a normalisation constant per image. Applies and outputs various pdfs and csv.
