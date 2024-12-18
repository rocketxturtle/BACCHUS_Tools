# BACCHUS_Tools
Python functions for automating BACCHUS

<b>Step 1:</b> open bacchus_tools.py and change BACCHUS_DIR to point to your BACCHUS directory

<b>Step 2:</b> open up a terminal and cd to your BACCHUS_DIR (<code>cd /path/to/BACCHUS</code>)

<b>Step 3:</b> make sure your terminal is searching for executables in your current working directory (e.g., type <code>PATH=$PATH:./</code> or add your BACCHUS_DIR to your .bash_profile)

<b>Step 4:</b> open init.com (e.g., <code>vim init.com</code>) and ensure everything is as you'd like (e.g., are you solving for Teff, logg, [M/H], vmic, and convol? are your linelists and lineselections correct? etc.)

<b>Step 5:</b> open stellar_parameters.tab (e.g., <code>vim stellar_parameters.tab</code>) and add in a line corresponding to the star whose spectrum you'd like to measure. Save edits and return to terminal.

&emsp; &emsp; NOTE: stellar_parameters.tab MUST be tab-delimited or some of these BACCHUS tools won't work.  so each line should read:
  
&emsp; &emsp; <star_name>\t<path_to_spectrum>\t<Teff_init>\t<logg_init>\t<[M/H]_init>\t<vmic_init (-99 if unknown) >\t<convolution_init>\t<rv_init (MUST BE 0 since spectrum must be rv-corrected) >

&emsp; &emsp; WHERE \t IS LITERALLY A TAB

<b>Step 6:</b> Start an iPython instance (type <code>ipython</code>)

<b>Step 7:</b> import your functions:
<code>
import sys
sys.path.append("/path/to/BACCHUS_Tools")
import bacchus_tools as b
</code>

<b>Step 8:</b> let's get the stellar parameters for a star that you named "star" in stellar_parameters.tab:
<code>
b.get_star_param(star)
</code>


<b>Step 9:</b> let's iterate on the stellar parameters until they converge (giving up after 10 attempts)
<code>
for i in range(10):
&emsp; b.redo_if_necessary(star)
</code>

<b>Step 10:</b> let's get the star's abundances in Fe, Mg, and Nd:
<code>
b.get_abund(star, elements=['Fe', 'Mg', 'Nd'])
</code>

<b>Step 11:</b> now, let's extract all the abundances and produce some summary tables:
<code>
b.get_bracket_abunds(star)
</code>

<b>Step 12:</b> now, let's get differential stellar parameters for star2 with respect to star (see, e.g., Yong et al. 2023).  First edit stellar_parameters.tab to include info for star2, and then run:
<code>
b.run_star_diff(star2, star)
for i in range(10):
&emsp; b.redo_diff_if_necessary(star2, star)
</code>

<b>Step 13:</b> finally, let's perform a line-by-line differential abundance analysis of star2 against star and get out abundance differences between the two stars in each line:
<code>
r.get_diff(star, BACCHUS_DIR, star2, BACCHUS_DIR)
</code>
