# Takes a dataframe containing columns with the names "Sample", "Group", "compound_names", and "into" (for peak area). 
# Runs ANOVA and outputs the results of Tukey Test and Kruscal-Wallis test in csv by sample and test. 
# Each sample will generate a KW and Tukey file. This function does not return an value. 
# Function also takes an optional parameter for an output_filepath. 
# By default the output will be generated in the same folder as the input. 
anova_tukey_kw <- function(rawData, output_filepath)
{
  # if an output location has been specified, sets wd to new output location. And gets input directory. 
  if(!missing(output_filepath))
  {
    startDirectory <- getwd()
    setwd(output_filepath)
  }
  
  # loop over each sample in the experiment
  for (sample in unique.default(rawData$Sample))
  {
    # set up dataframes to hold output data to be dumped into csv 
    resultsAnova <- data.frame()
    resultsTurkey <- data.frame()
    resultsKruscalWallis <- data.frame()
    
    # Split data according to sample
    sampleSubsetter <- rawData$Sample == sample
    sampleSubset <- rawData[sampleSubsetter,]
    
    # loop over each compound in the sample
    for (compound in unique.default(sampleSubset$compound_names))
    {
      # Makes a data subset for each compound as it runs through the loop
      compoundSubsetter <- sampleSubset$compound_names == compound
      compoundSubset <- sampleSubset[compoundSubsetter, ]
      
      # Runs ANOVA using peak area (into) as output variable and Group as input
      aov.out <- aov(into ~ Group, compoundSubset)
      Sum_Sq_Groups <- summary(aov.out)[[1]][[1,"Sum Sq"]]
      Sum_Sq_Residuals <- summary(aov.out)[[1]][[2,"Sum Sq"]]
      Var_exp_Groups <- Sum_Sq_Groups/sum(Sum_Sq_Groups, Sum_Sq_Residuals)*100
      Var_exp_Residuals <- 100 - Var_exp_Groups
      anova.row <- cbind(compound, summary(aov.out)[[1]][[1,"Df"]],summary(aov.out)[[1]][[1,"F value"]],summary(aov.out)[[1]][[1,"Pr(>F)"]],
      round(Var_exp_Groups,1), round(Var_exp_Residuals,1))

      # Post Hoc test - Tukey HSD - default
      turkeyTemp <- TukeyHSD(aov.out)
      turkey.row <- cbind(rownames(turkeyTemp$Group), compound, turkeyTemp$Group)
      
      # Kruskal-Wallis Test
      kw.out <- kruskal.test(into ~ Group, compoundSubset)
      kw.row <- cbind(compound, kw.out$statistic, kw.out$parameter, kw.out$p.value)
      
      # Store results in dataframe outside loop
      resultsAnova <- rbind(resultsAnova, anova.row)
      resultsTurkey <-  rbind(resultsTurkey, turkey.row)
      resultsKruscalWallis <- rbind(resultsKruscalWallis, kw.row)
    }
    
    # format output data
    colnames(resultsAnova) <- c("compound", "df", "F value", "p-value", "Var. explained_groups (%)","Var. explained_residulas (%)")
    colnames(resultsKruscalWallis) <- c("compound", "chi-squared", "df", "p-value")
    colnames(resultsTurkey) <- c("comparison", "compound", "diff", "lwr", "upr", "p-adj")
    
    # write CSV
    write.csv(resultsAnova, row.names = FALSE, file = paste0(sample, "-anova.csv"))
    write.csv(resultsKruscalWallis, row.names = FALSE, file = paste0(sample, "-kw.csv"))
    write.csv(resultsTurkey, row.names = FALSE, file = paste0(sample, "-tukey.csv"))
  }
  
  # reset to original directory
  if(!missing(output_filepath))
  {
    setwd(startDirectory)  
  }
}

##########################################################################################################################
# Example Usage for anova_tukey_kw function                                                                              #
##########################################################################################################################

# assign variable for input path and input file name
inPath <- "C:/Users/Hegeman Lab/Desktop/Data/ANOVA"
dataFile <- "input.csv"

# assign variable for output path
outPath <- "C:/Users/Hegeman Lab/Desktop/Output/ANOVA_results"

# get input data. note, the column names for the input data need to be "into", "Group", "compound_name", and "Sample"
setwd(inPath)
rawData <- read.csv(dataFile)

# View input data to see proper column names
View(rawData)

# Call function. This call uses outpath but the parameter is optional. 
# By default the csvs will be written to the current directory. 
anova_tukey_kw(rawData = rawData, output_filepath = outPath)

