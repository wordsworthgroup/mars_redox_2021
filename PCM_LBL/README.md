Planetary Climate Model (Line-By-Line) v1.0
----------------------------------------------

Written by R. Wordsworth (2017-2021), with additional input from F. Ding.

PCM_LBL is a 1D radiative-convective code designed to simulate the climates of diverse planetary atmospheres. The code is written in modular modern Fortran and uses a 'brute-force' spectral approach where absorption coefficients are computed on a fixed spectral grid directly from line data. This allows climate calculations to be performed more simply and at higher accuracy than in a correlated-k approach. 

If you use this code in your research, please cite: Wordsworth, R., Knoll, A. H., Hurowitz, J., Baum, M., Ehlmann, B. L., Head, J. W., Steakley, K.. A coupled model of episodic warming, oxidation and geochemical transitions on early Mars. Nature Geoscience (2021).

To run the code, you first need to download this repository to your home directory. Then, navigate to the HITRAN site (https://hitran.org/) and request a free account. Once your account is activated, download line-by-line absorption data for the species that are present in your planet's atmosphere. By default, the code is set up to run with H2O and CO2, so you should download these species at a minimum. Rename the .par HITRAN files HITRAN_H2O.par and HITRAN_CO2.par for H2O and CO2, respectively, and place them in the data/ directory of your repository. If you want to include more species, you can refer to source/cross_section.f90 and modify the source code directly.

Once you have put the HITRAN gas absorption data files in the data/ directory, navigate to example_run/. By default, the source and object directory lines in the Makefile are:

dir = /Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/source

obj = /Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/objects

You will need to change these to match the location of PCM_LBL on your own system. In addition, you must update the following definition in source/params_consts.f90:

character(len=100), save :: datadir='/Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/data/'

to match the data directory on your system where you saved the HITRAN data. Once this is done, type 'make all' while in the example_run/ directory to compile the code. 

After compiling, you can type './PCM_LBL.e' to run the resulting executable. 

The namelist file 'input.nml' contains a list of code parameters that can be modified. By default, the code is set up to run an early Mars case with a CO2-dominated 1 bar atmosphere. The code outputs results as ascii files to saved_sigma_data/ and results/.

One of the strengths of this model is that it allows the user to iterate rapidly between fast, lower accuracy calculations and slow high accuracy calculations. By default, the model is set up to run fairly fast at moderate resolution. If you want to increase the accuracy of the calculation, you can modify the values of the integer constants nS and nTem in source/params_consts.f90 (see comments in that portion of the code for further details).
