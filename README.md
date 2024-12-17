# BACCHUS_Tools
Python functions for automating BACCHUS

Step 1: open bacchus_tools.py (e.g., <code>vim bacchus_tools.py</code>) and change BACCHUS_DIR to point to your BACCHUS directory

Step 2: open up a terminal and cd to your BACCHUS_DIR (<code>cd /path/to/BACCHUS</code>)

Step 3: make sure your terminal is searching for executables in your current working directory (e.g., type <code>PATH=$PATH:./</code> or add your BACCHUS_DIR to your .bash_profile)

Step 4: open init.com (e.g., <code>vim init.com</code>) and ensure everything is as you'd like (e.g., are you solving for Teff, logg, [M/H], vmic, and convol? are your linelists and lineselections correct? etc.)

Step 5: open stellar_parameters.tab (e.g., <code>vim stellar_parameters.tab</code>) and add in a line corresponding to the star whose spectrum you'd like to measure. Save edits and return to terminal.

  ****NOTE: stellar_parameters.tab MUST be tab-delimited or some of these BACCHUS tools won't work.  so each line should read:
  
  <star_name>\t<path_to_spectrum>\t<Teff_init>\t<logg_init>\t<[M/H]_init>\t<vmic_init (-99 if unknown) >\t<convolution_init>\t<rv_init (MUST BE 0 since spectrum must be rv-corrected) >

  WHERE \t IS LITERALLY A TAB

Step 6: (my preferred method) Start an iPython instance (type <code>ipython</code>)

Step 7: import your functions:
<code>
import sys
sys.path.append("/path/to/BACCHUS_Tools")
import bacchus_tools as b
</code>

Step 8: let's get the stellar parameters for a star that you named "star" in stellar_parameters.tab:
<code>
b.get_star_param(star)
</code>


Step 9: let's iterate on the stellar parameters until they converge (giving up after 10 attempts)
<code>
for i in range(10):
&emsp; b.redo_if_necessary(star)
</code>

Step 10: let's get the star's abundances in Fe, Mg, and Nd:
<code>
b.get_abund(star, elements=['Fe', 'Mg', 'Nd'])
</code>

Step 11: now, let's get differential stellar parameters for star2 with respect to star to get differential stellar parameters (see, e.g., Yong et al. 2023).  First edit stellar_parameters.tab to include info for star2, and then run:
<code>
b.run_star_diff(star2, star)
for i in range(10):
&emsp; b.redo_diff_if_necessary(star2, star)
</code>

COMING SOON: making summary result tables
