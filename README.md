# checking and comparing chemical networks

In this directory, several python routines are provided to check the consistency of chemical networks, and compare chemical networks with each other. I re-wrote these routines in python using pandas dataframes to significantly increase the computational efficiency. As an example, the chemical network used in Taquet, Furuya, Walsh, & van Dishoeck (2016) and updated with some new reactions together with a corrupted one are provided. 


## check_network.py

This routine checks the consistency between the species and reaction files that are usually used by astrochemical models to describe a chemical network. You need to specify the name of the species and reactions files, and the suffix name of the output files as inputs and  check_network.py (```python3 check_network.py```) will check that:
- there is no duplicated species in the species file
- all species in the reactions file exist in the species file
- the mass conservation of each reaction
- there is no duplicated reactions in the reactions file. Two reactions are considered as duplicated if 1) they have the same reactants and products, 2) they are of similar type, 3) their temperature ranges overlap.

The species and reactions files start with an arbitrary number of commented lines, starting with a ! character. 
The species file must start with a line describing the different columns and usually consists of: the number, the name, the charge, the number of elements, and the initial abundance of each species. 
Similarly, the reactions file starts with a line describing the different columns and must consist of: two or three reactants, two to four products, followed by the alpha, beta, gamma parameters, and the type of reaction. Optionally, the temperature range, the type of reaction and formula, and the reaction number.

For both files, the python routine will use this line to determine the location of each column in the ASCII file, in particular regarding the length of the species characters. 

The routine will generate a series of .csv files that list the duplicated species, duplicated reactions, reactions that are not conservative, and reactions containing reactants/products not listed in the species file. 


## compare_networks.py

Once the chemical networks have been checked with the check_networy.py routine, it might be useful to compare them to have a summary of the missing reactions among the networks or the rates or temperature ranges that have been updated. To this aim, run compare_networks.py with python3 (```python3 compare_networks.py```) to check:
- the reactions that are missing in one of the two networks
- the reactions with different rates (through the alpha, beta, and gamma parameters)
- the reaction with different temperature ranges.

Note that the chemical networks to be tested need to be checked previously in order to avoid duplicated (i.e. same reactants, products, and reaction types) reactions within each network. As for check_network.py, the reactions files start with a line describing the different columns and must consist of: two or three reactants, two to four products, followed by the alpha, beta, gamma parameters, and the type of reaction. Optionally, the temperature range, the type of reaction and formula, and the reaction number. Obviouslyhe two chemical networks to be test need to have the same columns.

The routine will generate a series of .csv files that list the reactions that are missing in each chemical network with respect to the other one, the reactions with different rates (alpha, beta, gamma parameters), and reactions with different temperature ranges. 