##Problem set 7
#1-a) Get the length of my vector:
myVector <- c(2,10,34)
length(myVector)

#1-a) Get the length of my list:
myList <- list(30,10,4,5)
length(myList)

#1-b) Create range(a,b), with a = 1 and b = 10:
myRange <- seq(1,10)
myRange

#1-c) Initialize a vector (of length 10) of all 0s:
VectorZero <- rep(0,10)
VectorZero

#2-a)
v1 <- c(3,6,7,3,1)
v2 <- c(6,3,0,6,1)

#2-b)
v1+v2 # Outcome: 9 9 7 9 2
v1*v2 # Outcome: 18 18 0 18 1
v1-v2 # Outcome: -3 3 7 -3 0
v1/v2 # Outcome: 0.5 2.0 Inf 0.5 1.0, Inf is because I divided a number by 0.

#2-c) I do get the kind of result I expect.
10^2 # Outcome: 100
sqrt(100) # Outcome: 10
log(24) # Outcome: 3.178054

#3-a)
MatrixMult1 <- v1%*%v2
myExpectation1 <- 3*6+3*6+7*0+3*6+1*1
MatrixMult1 == expectation # TRUE

#3-b)
m <- matrix(c(3,7,6,1), nrow=2, ncol=2)
v = c(3,1)

#3-c)
MatrixMult2 <- v%*%m
myExpectation2 <- c(3*3+1*7,3*6+1*1)
MatrixMult2 == myExpectation2 # TRUE

MatrixMult3 <- m%*%v
myExpectation3 <- c(3*3+6*1,7*3+1*1)
MatrixMult3 == myExpectation3 # TRUE

#4-c)
# Download library:
library(data.table)

# Set local path:
readPath <- setwd( "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW7")

# Read data table:
data <- fread("HW7_4_data.txt") #creates a data.table
data

#4-d)
# Compute the log-fold change of every gene (i.e., log2(BY_expression/RM_expression)):
log_change <- data[, .(gene_name, log_change = log2(BY_expression/RM_expression))] 
head(log_change) # We get a lot of NaN (i.e., values divided by "0") because a lot of the RM expression values are "0".

#4-e)
# Filter out the "bad" rows, and redo part d:
log_change_filtered <- data[with(data, RM_expression > 0 & BY_expression > 0), .(gene_name, log_change = log2(BY_expression/RM_expression))] 
head(log_change_filtered)
# If I wanted to update data table by removing RM_expression and BM_expression rows with value = 0:
data_filtered <- data[with(data, RM_expression > 0 & BY_expression > 0), .(gene_name, gene_length, BY_expression, RM_expression)]

#4-f)
# Add pseudo-counts to every gene:
data_pseudo <- data[, .(gene_name, gene_length, BY_pseudo = BY_expression + 1,RM_pseudo = RM_expression + 1)] 
head(data_pseudo)

#data_pseudo_2 <- data[, BY_pseudo := BY_expression + 1] # := adds this information as a new column to the data.table "data"

#4-g)
# Re-compute the log-fold change using pseudo-counted data:
log_change_pseudo <- data_pseudo[, .(gene_name, log_change = log2(BY_pseudo/RM_pseudo))] 
head(log_change_filtered)

#4-h)
# Compute FPKM values for every gene:
FPKM <- data_filtered[, .(gene_name, FPKM_BY = BY_expression/(gene_length*sum(BY_expression)*(10^9)), FPKM_RM = RM_expression/(gene_length*sum(RM_expression)*(10^9)))]
head(FPKM)

#4-i)
# Compute the log-fold change of the FPKM values:
log_change_FPKM <- FPKM[, .(gene_name, log_change_FPKM = log2(FPKM_BY/FPKM_RM))] 
head(log_change_FPKM)
