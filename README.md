# ArboEpiSim
ArboEpiSim is an arbovirus agent-based epidemiological simulator taking into account within-host dynamics in R.

This piece of code has been used in the following publication (please refer to it for more informations):
- Albin Fontaine, Sebastian Lequime, Isabelle Moltini-Conclois, Davy Jiolle, Isabelle Leparc-Goffart, Robert Charles Reiner, Louis Lambrechts (2018) Epidemiological significance of dengue virus genetic variation in mosquito infection
dynamics. PLoS Pathogens (in press)

The three different files can be found in the repository. They differ in the way "Mosquito to human" (or vertebrate host) transmission probability is computed:
- ArboEpiSim-3Parameter.R : the probability of transmission is depend on time according to a 3 parameter model (see above paper for more details). Each parameter set is uniquely drawn for each individual mosquito in the population.
- ArboEpiSim-Fixed_Threshold.R : the probability of transmission is constant (100%) once the EIP (fixed, same for every individual in the population) is reached.
- ArboEpiSim-Variable_Threshold.R : the probability of transmission is constant (100%) once the EIP (uniquely drawn from a distribution, according for measured parameters of the population, for each mosquito individual) is reached. 

This code runs with the help of the following packages:
- foreach: Microsoft & Steve Weston (2017). foreach: Provides Foreach Looping Construct for R. R package version 1.4.4. https://CRAN.R-project.org/package=foreach
- doParallel: Microsoft & Steve Weston (2017). doParallel: Foreach Parallel Adaptor for the 'parallel' Package. R package version 1.0.11. https://CRAN.R-project.org/package=doParallel
- plyr: Hadley Wickham (2011) The Split=Apply=Combine Strategy for Data Analysis. Journal of Statistical Software 40(1):1=29
- dplyr: Hadley Wickham, Romain Francois, Lionel Henry, Kirill Müller (2017) dplyr: A Grammar of Data Manipulation. R package version 0.7.4. https://CRAN.R-project.org/package=dplyr
- reshape2: Hadley Wickham (2007) Reshaping Data with the reshape Package. Journal of Statistical Software 21(10):1=20
- data.table:Matt Dowle & Arun Srinivasan (2017) data.table: Extension of 'data.frame'. R package version 1.10.4=3. https://CRAN.R-project.org/package=data.table
- stringr:Hadley Wickham (2017) stringr: Simple, Consistant Wrappers for Common String Operations. R package version 1.2.0 https://CRAN.R-project.org/package=stringr