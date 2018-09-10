# networks: Checking and comparing chemical networks

In this directory, several python routines are provided to check the consistency of chemical networks, and compare chemical networks with each other. I re-wrote these routines in python using pandas dataframes to significantly increase the computational efficiency. 

## check_network.py

This routine checks the consistency between the species and reaction files that are usually used by astrochemical models to describe a chemical network. You need to specify the name of the species and reactions files, and the suffit name of the output files as inputs and run check_network.py with python3 (```python3 check_network.py```) will check that:
- there is no duplicated species in the species file
- all species in the reactions file exist in the species file
- the mass conservation of each reaction
- there is no duplicated reactions in the reactions file. Two reactions are considered as duplicated if 1) they have the same reactants and products, 2) they are of similar type, 3) their temperature ranges overlap.

The species and reactions files start with an arbitrary number of commented lines, starting with a ! character. 
The species file must start with a line describing the different columns and usually consists of: the number, the name, the charge, the number of elements, and the initial abundance of each species. 
Similarly, the reactions file starts with a line describing the different columns and must consist of: two or three reactants, two to four products, followed by the alpha, beta, gamma parameters, and the type of reaction. Optionally, the temperature range, the type of reaction and formula, and the reaction number.

For both files, the python routine will use this line to determine the location of each column in the ASCII file, in particular regarding the length of the species characters. 



## compare_networks.py
