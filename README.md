# ANOVA
## Functionality
Takes a dataframe containing columns with the names "Sample", "Group", "compound_names", and "into" (for peak area in this example). Runs ANOVA and outputs the results of Tukey Test and Kruscal-Wallis test in csv by sample and test. Each sample will generate a KW and Tukey file. This function does not return an value. Function also takes an optional parameter for an output_filepath. By default the output will be generated in the same folder as the input. 

## Usage
* Update input CSV and save it. 
* Set the working directory to the folder containing the input csv and the data files
* Adjust the locations and filenames in lines 80, 81, and 84
* Run lines 6-73 to load in the needed functions
* Run lines 80-87 for input set up
* Run 88 to read in the input file
* (Optional) Run 91 to check the input data that has been read in.
* Run line 95 to generate output CSVs

## Hegeman Lab - University of Minnesota Twin-Cities
This code was developed for use in the Hegeman Lab at the University of Minnesota Twin-Cities. If you use this script in your research, please don't forget to acknowledge or cite publication. Additionally, if there are any questions about how to use this code, feel free to contact [Adrian Hegeman](mailto:hegem007@umn.edu). 
