There are three files for the three reaction models:
	Chemical_Degradation_final.ipynb
	Electrochemical_Degradation_final.ipynb
	Cuppled_Degradation_final.ipynb

All of these files are structured identically:
	importing libraries and extern functions (e. g. rate_const_H2O2)
	calculating values, which are later unchanging (number of reaction sites, fluoride concentration etc.)
	definition of constants and intern functions (rxn, gas_trans, deg_rate, H2_price)
	execution of the simulation with different sets of parameters in different visualizations

In order to execute the simulation, the function deg_rate is executed. It takes temperature T, current density curden, cathodic pressure P_cat, anodic pressure P_an and membrane thickness h_memb as input parameters. Depending on the file in use, deg_rate also takes different parameters relevant to this instance of the simulation (e. g. f_Fe, changing the iron source rate relevant to the electrochemical part but irrelevant to the chemical part). It is important to understand, that deg_rate only takes single values as input, never lists or arrays! Simulation over a list of a parameter has to be solved with a loop.
Depending on the input parameters, deg_rate calculates all dependant values and executes the reaction network rxn. Output of deg_rate is a list of different values, usually with attack_rate at position Zero with unit mol/(l*s).
Furthermore, Cuppled_degradation_final.ipynb contains the function H2_price. This takes temperature T, current density curden, cathodic pressure P_cat, anodic pressure P_an, membrane thickness h_memb, active membrane area A, lifetime_h (in hours !) and energy_price (in €/Wh !) as input parameters. Output is a list of the H2-price per kilogram, total mass of produced H2, CAPEX, OPEX and energy cost.


More files are:
	H2O2_formation.ipynb; This simulates and validates the cell model. The function rate_const_H2O2 is created and takes temperature T, cathodic pressure P_cat, anodic pressure P_an and membrane thickness h_memb as input parameters. Output is a list with another list a position Zero. This inner list contains the reaction rate constant for the H2O2-formation in dependance of 10,000 values of the current density between 0 and 3 A/cm². To obtain the reaction rate constant at a certain current density, the value of the 10,000 is determined, which lies closest to the input value. The index of this approximation value is used to obtain the reaction rate constant at that value. > k=rate_const_H2O2(T,P_cat,P_an,h_memb)[0][abs(j-curden).argmin()]
	gas_trans_v2.ipynb; This simulates the mass transport through the membrane following Schalenbach. It contains the function gas_trans used by deg_rate in the simulations to calculate certain values. It takes temperature T, current density curden, cathodic pressure P_cat, anodic pressure P_an and membrane thickness h_memb as input parameters.


At different occasions data from the literature is being used for comparison. This data is contained in the following text files:
	Chandesris_FER.txt; This contains x-values [0:6] and y-values [6:12] of FER vs. curden at 60 °C as well as x-values [12:15] and y-values [15:18] of FER vs. curden at 80 °C. Fig. 3 from DOI: 10.1016/j.ijhydene.2014.11.111
	Hellsing_desorption.txt; This contains x-values [0:19] and y-values [19:38] of OH-desorption rates vs. alpha as well as x-values [38:67] and y-values [67:96] of H2O-desorption rates vs. alpha. Fig. 1 from DOI: 10.1016/0021-9517(91)90258-6
	Laconti_FER.txt; This contains x-Values [0:6] and y-values [6:12] of FER vs. T as well es two x-values [12:14] and two y-values [14:16] for the regression line. Fig. 9 from DOI: 10.1149/1.2214554
	Schalenbach_crossover.txt; For H2 in O2 vs. curden this contains x-values [0:10] and y-values [10:20] at equal pressure as well as x-values [20:29] and y-values [29:38] at differential pressure. Fig. 3 aus DOI: 10.1016/j.ijhydene.2013.09.013