######----- INSTALLING AND LOADING NEEDED PACKAGES---- 

#install.packages("tidyverse")
library(tidyverse)
library(dplyr)

#install.packages("coil")
library(coil)

#install.packages("seqinr")
library(seqinr)

#Install Biostrings
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library("Biostrings")

#install.packages("data.table")
library(data.table)



######----- ACQUIRING TARDIGRADA DATA FROM BOLD---- 

Tardigrada <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Tardigrada&format=tsv")
#Result: 1635 observations with 80 variables


######----- FILTERING TARDIGRADA COI-5P SEQUENCES FROM BOLD---- 

#Criteria:
#1)Folmer region of the sequence must be greater that 600 base pairs
#2)The taxonomy of the sequence must be known at least to the genus level
#3)There are no missing base pairs in the sequence


Tardigrada.COI <- Tardigrada %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[BDEFHIJKLMNOPQRSUVWXYZ]", TRUE)) %>% #Criteria 3 + removing any nucleotide sequences that contain any letters besides ACGT
  filter(!is.na(genus_name)) %>% #Criteria 2
  filter(nchar(nucleotides) > 600) #Criteria 1 (Maybe, as this is the whole COI gene not simply the Folmer region)
  

#Here I remove dashes in COI sequences on the edges and the inside (to be able to access the Folmer region later on)
#Tardigrada.COI$nucleotides <- gsub("-", "", Tardigrada.COI$nucleotides)
#I don't do this anymore

#Histogram of nucleotide sequence lengths
hist(nchar(Tardigrada.COI$nucleotides),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I apply the coil pipeline

Tardigrada.COI_coi5p_pipe <- lapply(1:length(Tardigrada.COI$processid), function(i){
  coi5p_pipe(Tardigrada.COI$nucleotides[i],
             name = Tardigrada.COI$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(Tardigrada.COI_coi5p_pipe)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Tardigrada_coi5p_df <- flatten_coi5p(Tardigrada.COI_coi5p_pipe)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Tardigrada_coi5p_df$framed),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES----

#Here I extract the framed sequences into their own dataframe

Tardigrada_coi5p_framed <- subset(Tardigrada_coi5p_df, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Tardigrada_coi5p_framed$framed <- toupper(Tardigrada_coi5p_framed$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters
Tardigrada_coi5p_framed$framed <- gsub("^.{0,50}", "", Tardigrada_coi5p_framed$framed)

#Here I remove sequences containing dashes 
Tardigrada_coi5p_framed <- Tardigrada_coi5p_framed %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Tardigrada_coi5p_framed$framed),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


#Here I generate the FASTA file

#This creates an empty character vector, with length twice the length of the dataframe; then puts the values from column1 in every second position starting at 1, and the values of column2 in every second position starting at 2.

Tardigrad_fasta <- character(nrow(Tardigrada_coi5p_framed) * 2)
Tardigrad_fasta[c(TRUE, FALSE)] <- paste0(">", Tardigrada_coi5p_framed$name)
Tardigrad_fasta[c(FALSE, TRUE)] <- Tardigrada_coi5p_framed$framed

#Then write using writeLines:
  
writeLines(Tardigrad_fasta, "Tardigrada_coi5p.fasta") 


######----- BAR PLOTS FOR # OF DIFFERENT TAXANOMIC GROUPS----

#Class

barplot(table(Tardigrada.COI$class_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Classes of Tardigrada")
#Order

barplot(table(Tardigrada.COI$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Tardigrada")

#Family

barplot(table(Tardigrada.COI$family_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Families of Tardigrada",
        las=2,
        cex.names = 0.6,
        srt = 60)

#Better looking ggplot barplots (these contain N/A taxanomic groups):

#Class

ggplot(data = Tardigrada.COI) + 
  geom_bar(mapping = aes(x = class_name), stat = "count", fill = "turquoise") +
  labs(title = "Classes of Tardigrada", x = "Taxanomic Group", y = "Counts")

#Order

ggplot(data = Tardigrada.COI) + 
  geom_bar(mapping = aes(x = order_name), stat = "count", fill = "turquoise") +
  labs(title = "Orders of Tardigrada", x = "Taxanomic Group", y = "Counts")

#Family

ggplot(data = Tardigrada.COI) + 
  geom_bar(mapping = aes(x = family_name), stat = "count", fill = "turquoise") +
  labs(title = "Families of Tardigrada", x = "Taxanomic Group", y = "Counts") +
  theme(axis.text.x = element_text(angle = 90))



######----- REMOVING 2 MINOR ORDERS and BALANCING THE QUANTITY OF SEQUENCES FOR THE 2 MAJOR ORDERS----

#There are 4 Orders: Apochela (minor), Arthrotardigrada (minor), Echiniscoidea (major), Parachaela (major) (& a couple NAs)

Tardigrada.COI_2Orders <- Tardigrada.COI %>%
  filter(str_detect(order_name, "Apochela", TRUE)) %>%
  filter(str_detect(order_name, "Arthrotardigrada", TRUE))

#Class

barplot(table(Tardigrada.COI_2Orders$class_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Classes of Tardigrada")
#Order

barplot(table(Tardigrada.COI_2Orders$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Tardigrada")

#Family

barplot(table(Tardigrada.COI_2Orders$family_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Families of Tardigrada",
        las=2,
        cex.names = 0.6,
        srt = 60)

table(Tardigrada.COI_2Orders$order_name)
#There are 259 Echiniscoidea and 575 Parachaela


###Here I duplicate the Echiniscoidea column

#This creates the new column named "duplicate" filled with "NA"
Tardigrada.COI_2Orders["duplicate"] <- NA

#https://stackoverflow.com/questions/28961442/create-duplicate-rows-based-on-conditions-in-r
#This creates duplicate rows containing Echiniscoidea Order and asigns a value of "2" to the "duplicate" coloumn for the duplicate rows
Tardigrada.COI_2Orders <- rbind(Tardigrada.COI_2Orders,
      Tardigrada.COI_2Orders %>% 
        filter(order_name == "Echiniscoidea") %>% 
        mutate(duplicate = 2,))

table(Tardigrada.COI_2Orders$order_name)
#Now there are 518 Echiniscoidea and 575 Parachaela


######----- Applying the coil pipeline again WITH SEQUENCES CONTAINING ONLY 2 ORDERS----

Tardigrada.COI_coi5p_pipe_2Orders <- lapply(1:length(Tardigrada.COI_2Orders$processid), function(i){
  coi5p_pipe(Tardigrada.COI_2Orders$nucleotides[i],
             name = Tardigrada.COI_2Orders$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(Tardigrada.COI_coi5p_pipe_2Orders)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Tardigrada_coi5p_df_2Orders <- flatten_coi5p(Tardigrada.COI_coi5p_pipe_2Orders)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Tardigrada_coi5p_df_2Orders$framed),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES AGAIN----

#Here I extract the framed sequences into their own dataframe

Tardigrada_coi5p_framed_2Orders <- subset(Tardigrada_coi5p_df_2Orders, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Tardigrada_coi5p_framed_2Orders$framed <- toupper(Tardigrada_coi5p_framed_2Orders$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters of all the sequences even if they do not have a dash.
Tardigrada_coi5p_framed_2Orders$framed <- gsub("^.{0,50}", "", Tardigrada_coi5p_framed_2Orders$framed)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Tardigrada_coi5p_framed_2Orders$framed),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I remove sequences containing dashes 
Tardigrada_coi5p_framed_2Orders <- Tardigrada_coi5p_framed_2Orders %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Tardigrada_coi5p_framed_2Orders$framed),
     main = paste("Tardigrada COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


#Here I generate 1 FASTA file

#This creates an empty character vector, with length twice the length of the dataframe; then puts the values from column1 in every second position starting at 1, and the values of column2 in every second position starting at 2.

Tardigrad_fasta_2Orders <- character(nrow(Tardigrada_coi5p_framed_2Orders) * 2)
Tardigrad_fasta_2Orders[c(TRUE, FALSE)] <- paste0(">", Tardigrada_coi5p_framed_2Orders$name)
Tardigrad_fasta_2Orders[c(FALSE, TRUE)] <- Tardigrada_coi5p_framed_2Orders$framed

#Then write using writeLines:

writeLines(Tardigrad_fasta_2Orders, "Tardigrada_coi5p_2Orders.fasta") 

#Here I make 2 FASTA files. One with 10% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#It will split your original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Tardigrada_coi5p_framed_2Orders)/2)
df1 = Tardigrada_coi5p_framed_2Orders[1:index,]
df2 = Tardigrada_coi5p_framed_2Orders[(index +1) : nrow(Tardigrada_coi5p_framed_2Orders),]

#In this case df1 has 496 observations and df2 has 497 observations

##df1 to FASTA (this wil be the file I add errors to):

Tardigrad_fasta_1 <- character(nrow(df1) * 2)
Tardigrad_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Tardigrad_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Tardigrad_fasta_1, "Tardigrada_1.fasta") 

##df2 to FASTA (this wil be the error free file):

Tardigrad_fasta_2 <- character(nrow(df2) * 2)
Tardigrad_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Tardigrad_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Tardigrad_fasta_2, "Tardigrada_clean.fasta") 



######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Tardigrada_error.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Tardigrada_clean.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#We will replace ACGT with 1234 in both dataframes

df1_error$sequence <- gsub('A', '1', df1_error$sequence)
df1_error$sequence <- gsub('C', '2', df1_error$sequence)
df1_error$sequence <- gsub('G', '3', df1_error$sequence)
df1_error$sequence <- gsub('T', '4', df1_error$sequence)

df2_clean$sequence <- gsub('A', '1', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '2', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '3', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '4', df2_clean$sequence)

#Combine the 2 dataframes

df3 <- rbind(df1_error, df2_clean)


#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled <- df3[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled <- df3_shuffled[ , !(names(df3_shuffled) %in% drop)]


#I want to split the "sequence" column into multiple columns.

#First, let's find the maximum and minimum sequence length in the sequence column

max(nchar(df3_shuffled$sequence)) #615
#615/3 = 205 (we can split into 3 characters per column)

min(nchar(df3_shuffled$sequence)) #553
#553/7 = 79 (we can split into 7 characters per column)

#This will split based on the size of the minimum sequence length

tmp <- read.fwf(
  textConnection(df3_shuffled$sequence),
  widths = rep(7, ceiling(max(nchar(df3_shuffled$sequence) / 7))),
  stringsAsFactors = FALSE)

df3_shuffled_split <- cbind(df3_shuffled, tmp)

#Remove 'sequence' label 

df3_shuffled_split$sequence <- NULL

#Ok, so it seems like this did not split based off the minimum sequence length (553). But we were still able to split 7 characters per column. So for now, just manually erase all columns after 79 so that there is no column with NA.

df3_shuffled_split <- df3_shuffled_split[, -c(81:89)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split, file = "Tardigrada_2Orders.csv", row.names = FALSE)






#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_1 <- cbind(df3_shuffled, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_1$sequence <- NULL

#Ok, so it seems like this did not split based off the minimum sequence length (553). But we were still able to split 1 character per column. So for now, just manually erase all columns after 554 so that there is no column with NA.

df3_shuffled_split_1 <- df3_shuffled_split_1[, -c(555:616)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_1, file = "Tardigrada_2Orders_1.csv", row.names = FALSE)





#Now I try one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('1', '0001', df1_error$sequence)
df1_error$sequence <- gsub('2', '0010', df1_error$sequence)
df1_error$sequence <- gsub('3', '0100', df1_error$sequence)
df1_error$sequence <- gsub('4', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('1', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('2', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('3', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('4', '1000', df2_clean$sequence)

#Combine the 2 dataframes

df3_oh <- rbind(df1_error, df2_clean)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh <- df3_oh[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled_oh <- df3_shuffled_oh[ , !(names(df3_shuffled_oh) %in% drop)]

#This will split based on the size of the minimum sequence length

tmp_oh <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp_oh)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#Ok, so it seems like this did not split based off the minimum sequence length (2212). But we were still able to split 1 character per column. So for now, just manually erase all columns after 2213 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2214:2461)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Tardigrada_2Orders_oh.csv", row.names = FALSE)








#Here I make 2 FASTA files. One with 5% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#It will split your original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Tardigrada_coi5p_framed_2Orders)/2)
df1 = Tardigrada_coi5p_framed_2Orders[1:index,]
df2 = Tardigrada_coi5p_framed_2Orders[(index +1) : nrow(Tardigrada_coi5p_framed_2Orders),]

#In this case df1 has 496 observations and df2 has 497 observations

##df1 to FASTA (this wil be the file I add errors to):

Tardigrad_fasta_1 <- character(nrow(df1) * 2)
Tardigrad_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Tardigrad_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Tardigrad_fasta_1, "Tardigrada_1_5per.fasta") 

##df2 to FASTA (this wil be the error free file):

Tardigrad_fasta_2 <- character(nrow(df2) * 2)
Tardigrad_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Tardigrad_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Tardigrad_fasta_2, "Tardigrada_clean_5per.fasta") 



######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Tardigrada_error_5per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Tardigrada_clean_5per.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#We will replace ACGT with 1234 in both dataframes

df1_error$sequence <- gsub('A', '1', df1_error$sequence)
df1_error$sequence <- gsub('C', '2', df1_error$sequence)
df1_error$sequence <- gsub('G', '3', df1_error$sequence)
df1_error$sequence <- gsub('T', '4', df1_error$sequence)

df2_clean$sequence <- gsub('A', '1', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '2', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '3', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '4', df2_clean$sequence)

#Combine the 2 dataframes

df3 <- rbind(df1_error, df2_clean)


#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled <- df3[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled <- df3_shuffled[ , !(names(df3_shuffled) %in% drop)]



#Now I do one-hot encoding

#We will replace 1234 with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('1', '0001', df1_error$sequence)
df1_error$sequence <- gsub('2', '0010', df1_error$sequence)
df1_error$sequence <- gsub('3', '0100', df1_error$sequence)
df1_error$sequence <- gsub('4', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('1', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('2', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('3', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('4', '1000', df2_clean$sequence)

#Combine the 2 dataframes

df3_oh <- rbind(df1_error, df2_clean)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh <- df3_oh[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled_oh <- df3_shuffled_oh[ , !(names(df3_shuffled_oh) %in% drop)]

#This will split sequences to one character per column

tmp_oh <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp_oh)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#Ok, so it seems like this did not split based off the minimum sequence length (2212). But we were still able to split 1 character per column. So for now, just manually erase all columns after 2213 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2214:2461)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Tardigrada_2Orders_oh_5%.csv", row.names = FALSE)










######----- COMBINING THE 2 FASTA FILES AFTER 2% ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Tardigrada_error_2per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Tardigrada_clean.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#We will replace ACGT with 1234 in both dataframes

df1_error$sequence <- gsub('A', '1', df1_error$sequence)
df1_error$sequence <- gsub('C', '2', df1_error$sequence)
df1_error$sequence <- gsub('G', '3', df1_error$sequence)
df1_error$sequence <- gsub('T', '4', df1_error$sequence)

df2_clean$sequence <- gsub('A', '1', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '2', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '3', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '4', df2_clean$sequence)

#Combine the 2 dataframes

df3 <- rbind(df1_error, df2_clean)


#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled <- df3[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled <- df3_shuffled[ , !(names(df3_shuffled) %in% drop)]



#Now I do one-hot encoding

#We will replace 1234 with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('1', '0001', df1_error$sequence)
df1_error$sequence <- gsub('2', '0010', df1_error$sequence)
df1_error$sequence <- gsub('3', '0100', df1_error$sequence)
df1_error$sequence <- gsub('4', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('1', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('2', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('3', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('4', '1000', df2_clean$sequence)

#Combine the 2 dataframes

df3_oh <- rbind(df1_error, df2_clean)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh <- df3_oh[rows, ]


#Remove the seq_name column
drop <- c("seq_name")
df3_shuffled_oh <- df3_shuffled_oh[ , !(names(df3_shuffled_oh) %in% drop)]

#This will split sequences to one character per column

tmp_oh <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp_oh)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#Ok, so it seems like this did not split based off the minimum sequence length (2212). But we were still able to split 1 character per column. So for now, just manually erase all columns after 2213 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2214:2461)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Tardigrada_2Orders_oh_2%.csv", row.names = FALSE)





