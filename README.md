# BACCHUS_Tools
Python functions for automating BACCHUS

Step 1: open bacchus_tools.py and change BACCHUS_DIR to point to your BACCHUS directory

Step 2: open up a terminal and cd to your BACCHUS_DIR (cd /path/to/BACCHUS)

Step 3: make sure your terminal is searching for executables in your current working directory (e.g., type PATH=$PATH:./ or add your BACCHUS_DIR to your .bash_profile)

Step 4: open init.com (e.g., vim init.com) and ensure everything is as you'd like (e.g., are you solving for Teff, logg, [M/H], vmic, and convol? are your linelists and lineselections correct? etc.)

Step 5: open stellar_parameters.tab (e.g., vim stellar_parameters.tab) and add in a line corresponding to the star whose spectrum you'd like to measure. Save edits and return to terminal.

  ****NOTE: stellar_parameters.tab MUST be tab-delimited or some of these BACCHUS tools won't work.  so each line should read:
  
  <star_name>\t<path_to_spectrum>\t<Teff_init>\t<logg_init>\t<[M/H]_init>\t<vmic_init (-99 if unknown) >\t<convolution_init>\t<rv_init (MUST BE 0 since spectrum must be rv-corrected) >

Step 6: (my preferred method) Start an iPython instance (type ipython)

Step 7: import your functions:

>>>> import sys
>>>> sys.path.append("/path/to/BACCHUS_Tools")
>>>> import bacchus_tools as b

Step 8: let's get the stellar parameters for a star that you named "star" in stellar_parameters.tab:

>>>> b.get_star_param(star)

Step 9: let's get its abundances in Fe, Mg, and Nd:
>>>> b.get_abund(star, elements=['Fe', 'Mg', 'Nd'])

COMING SOON: making summary result tables
