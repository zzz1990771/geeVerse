data("simuGene")
sim_data <- generateData(nsub = 1000, nobs = rep(5,1000), p = 50,
                         beta0 = c(rep(1,9),rep(0,41)), rho = 0.6, correlation = "exchangeable",
                         ka = 0.5, SNPs = simuGene[,1:25])
y = sim_data$y
x = sim_data$X

yx = cbind(rep(1:5,each = 1000),rep(1:5,1000),y,x)
colnames(yx)[1:3] = c("subject","observe time", "response")

summary_table = as.data.frame(yx[1:7,c(1,2,3,4,5,8,10,11,17)])
colnames(summary_table)[c(5,8)] = c("...",'...')
summary_table[,c(5,8)] = "..."

# Load required libraries
library(dplyr)
library(knitr)
library(kableExtra)

# Create the table
table <- summary_table %>%
  kable("html", caption = "Simulated Longitudinal Genetics Data", align= "rrrrcrrcr",
        ,row.names = FALSE, digits = 2) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  add_header_above(c(" " = 2, " " = 1, "Controls" = 3, "Genetic Markers" = 3)) %>%
  column_spec(1:2, bold = TRUE) %>%
  column_spec(3:9, width = "70px")
  #pack_rows("subject 1", 1, 5)  # Assuming 5 rows for subject 1


# Print the table
print(table)




# Create the table for AIDS data
# load the "aids" dataset from the JM package
aids <- JM::aids
aids_table <- aids[1:8, c(1,5,6,7,8,9,4)]

table <- aids_table %>%
  kable("html", caption = "Longitudinal AIDS Data", align= "lrrrrrr",
        ,row.names = FALSE, digits = 2) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2:7, width = "70px")
#pack_rows("subject 1", 1, 5)  # Assuming 5 rows for subject 1
# Print the table
print(table)
