######----- INSTALLING AND LOADING NEEDED PACKAGES---- 

install.packages("tidyverse")
library(tidyverse)
library(dplyr)

#install.packages("coil")
library(coil)

#install.packages("seqinr")
library(seqinr)

#Install Biostrings
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library("Biostrings")

#install.packages("data.table")
library(data.table)


######----- ACQUIRING ANNELIDA DATA FROM BOLD---- 

Annelida <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Annelida&format=tsv")
Annelida
#Read on January 26, 2022


######----- FILTERING ANNELIDA COI-5P SEQUENCES FROM BOLD---- 

#Criteria:
#1)Folmer region of the sequence must be greater that 600 base pairs
#2)The taxonomy of the sequence must be known at least to the genus level
#3)There are no missing base pairs in the sequence


Annelida.COI <- Annelida %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[BDEFHIJKLMNOPQRSUVWXYZ]", TRUE)) %>% #Criteria 3 + removing any nucleotide sequences that contain any letters besides ACGT
  filter(!is.na(genus_name)) %>% #Criteria 2
  filter(nchar(nucleotides) > 600) #Criteria 1 (Maybe, as this is the whole COI gene not simply the Folmer region)


#Histogram of nucleotide sequence lengths
hist(nchar(Annelida.COI$nucleotides),
     main = paste("Annelida COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- BAR PLOTS FOR # OF DIFFERENT TAXANOMIC GROUPS----

#Class

barplot(table(Annelida.COI$class_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Classes of Annelida")
#Order

barplot(table(Annelida.COI$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Annelida")

#Family

barplot(table(Annelida.COI$family_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Families of Annelida",
        las=2,
        cex.names = 0.6,
        srt = 60)

#Better looking ggplot barplots (these contain N/A taxanomic groups):

#Class

ggplot(data = Annelida.COI) + 
  geom_bar(mapping = aes(x = class_name), stat = "count", fill = "turquoise") +
  labs(title = "Classes of Annelida", x = "Taxanomic Group", y = "Counts")

#Order

ggplot(data = Annelida.COI) + 
  geom_bar(mapping = aes(x = order_name), stat = "count", fill = "turquoise") +
  labs(title = "Orders of Annelida from BOLD", x = "Taxonomic Group", y = "Counts") + coord_flip() +
  theme(axis.text=element_text(size=15),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(angle = 0))

#Family

ggplot(data = Annelida.COI) + 
  geom_bar(mapping = aes(x = family_name), stat = "count", fill = "turquoise") +
  labs(title = "Families of Annelida", x = "Taxanomic Group", y = "Counts") +
  theme(axis.text.x = element_text(angle = 90))


######----- REMOVING ALL ORDERS EXCEPT THE 2 LARGEST----

Annelida_Order_Counts <- as.data.frame(table(Annelida.COI$order_name))

Annelida_Order_Counts <- Annelida_Order_Counts[order(-Annelida_Order_Counts$Freq), ]
Annelida_Order_Counts

#The two largest orders are:
#Haplotaxida 11307 & Phyllodocida  5789

Annelida_2Orders <- Annelida.COI %>%
  filter(str_detect(order_name, "Crassiclitellata", TRUE)) %>%
  filter(str_detect(order_name, "Terebellida", TRUE)) %>%
  filter(str_detect(order_name, "Sabellida", TRUE)) %>%
  filter(str_detect(order_name, "Arhynchobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Rhynchobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Eunicida", TRUE)) %>%
  filter(str_detect(order_name, "Spionida", TRUE)) %>%
  filter(str_detect(order_name, "Enchytraeida", TRUE)) %>%
  filter(str_detect(order_name, "Amphinomida", TRUE)) %>%
  filter(str_detect(order_name, "Lumbriculida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta_incertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Opheliida", TRUE)) %>%
  filter(str_detect(order_name, "Canalipalpata", TRUE)) %>%
  filter(str_detect(order_name, "Echiuroidea", TRUE)) %>%
  filter(str_detect(order_name, "Capitellida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta__order_incertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Myzostomida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta_insertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Branchiobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Aciculata", TRUE)) %>%
  filter(str_detect(order_name, "Acanthobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Capilloventrida", TRUE)) 

table(Annelida_2Orders$order_name)

barplot(table(Annelida_2Orders$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Annelida")


###Here I duplicate the Phyllodocida order

#This creates the new column named "duplicate" filled with "NA"
Annelida_2Orders["duplicate"] <- NA

#https://stackoverflow.com/questions/28961442/create-duplicate-rows-based-on-conditions-in-r
#This creates duplicate rows containing Phyllodocida Order and asigns a value of "2" to the "duplicate" coloumn for the duplicate rows
Annelida_2Orders <- rbind(Annelida_2Orders,
                                Annelida_2Orders %>% 
                                  filter(order_name == "Phyllodocida") %>% 
                                  mutate(duplicate = 2,))

table(Annelida_2Orders$order_name)
#Now there are 11307 Haplotaxida and 11578 Phyllodocida


######----- REMOVING ALL ORDERS EXCEPT FOR 3 MINOR----

#The 3 orders I will be isolating are Crassiclitellata (2161), Terebellida (2101), and Rhynchobdellida (1081)

Annelida_3Orders <- Annelida.COI %>%
  filter(str_detect(order_name, "Haplotaxida", TRUE)) %>%
  filter(str_detect(order_name, "Phyllodocida", TRUE)) %>%
  filter(str_detect(order_name, "Sabellida", TRUE)) %>%
  filter(str_detect(order_name, "Arhynchobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Eunicida", TRUE)) %>%
  filter(str_detect(order_name, "Spionida", TRUE)) %>%
  filter(str_detect(order_name, "Enchytraeida", TRUE)) %>%
  filter(str_detect(order_name, "Amphinomida", TRUE)) %>%
  filter(str_detect(order_name, "Lumbriculida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta_incertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Opheliida", TRUE)) %>%
  filter(str_detect(order_name, "Canalipalpata", TRUE)) %>%
  filter(str_detect(order_name, "Echiuroidea", TRUE)) %>%
  filter(str_detect(order_name, "Capitellida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta__order_incertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Myzostomida", TRUE)) %>%
  filter(str_detect(order_name, "Polychaeta_insertae_sedis", TRUE)) %>%
  filter(str_detect(order_name, "Branchiobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Aciculata", TRUE)) %>%
  filter(str_detect(order_name, "Acanthobdellida", TRUE)) %>%
  filter(str_detect(order_name, "Capilloventrida", TRUE)) 

table(Annelida_3Orders$order_name)

barplot(table(Annelida_3Orders$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Annelida")


######----- REMOVING ALL ORDERS EXCEPT HAPLOTAXIDA---- (don't actually need this)

Annelida_Haplotaxida <- Annelida_2Orders %>%
  filter(str_detect(order_name, "Phyllodocida", TRUE))

table(Annelida_Haplotaxida$order_name)



######----- REMOVING ALL ORDERS EXCEPT PHYLLODOCIDA---- (don't actually need this)

Annelida_Phyllodocida <- Annelida_2Orders %>%
  filter(str_detect(order_name, "Haplotaxida", TRUE))

table(Annelida_Phyllodocida$order_name)



######----- Applying the coil pipeline WITH ANNELIDA DATA CONTAINING ONLY 2 ORDERS----

annelida_coi5p_pipe_2Orders <- lapply(1:length(Annelida_2Orders$processid), function(i){
  coi5p_pipe(Annelida_2Orders$nucleotides[i],
             name = Annelida_2Orders$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(annelida_coi5p_pipe_2Orders)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Annelida_coi5p_df_2Orders <- flatten_coi5p(annelida_coi5p_pipe_2Orders)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Annelida_coi5p_df_2Orders$framed),
     main = paste("Annelida COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")



######----- Applying the coil pipeline WITH HAPLOTAXIDA ORDER DATA---- (don't actually need this)

annelida_coi5p_pipe_haplotaxida <- lapply(1:length(Annelida_Haplotaxida$processid), function(i){
  coi5p_pipe(Annelida_Haplotaxida$nucleotides[i],
             name = Annelida_Haplotaxida$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(annelida_coi5p_pipe_haplotaxida)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Annelida_coi5p_df_haplotaxida <- flatten_coi5p(annelida_coi5p_pipe_haplotaxida)


######----- Applying the coil pipeline WITH PHYLLODOCIDA ORDER DATA----(don't actually need this)

annelida_coi5p_pipe_phyllodocida <- lapply(1:length(Annelida_Phyllodocida$processid), function(i){
  coi5p_pipe(Annelida_Phyllodocida$nucleotides[i],
             name = Annelida_Phyllodocida$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(annelida_coi5p_pipe_phyllodocida)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Annelida_coi5p_df_phyllodocida <- flatten_coi5p(annelida_coi5p_pipe_phyllodocida)



##So, we will use the Annedlida sequences containing the 2 orders because those are the most common orders in the Annelida Phylum. Later we will use 3 minor orders as test data for the ML algorithms. 



######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES ----

#Here I extract the framed sequences into their own dataframe
Annelida_coi5p_framed_2Orders <- subset(Annelida_coi5p_df_2Orders, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Annelida_coi5p_framed_2Orders$framed <- toupper(Annelida_coi5p_framed_2Orders$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters of all the sequences even if they do not have a dash.
Annelida_coi5p_framed_2Orders$framed <- gsub("^.{0,50}", "", Annelida_coi5p_framed_2Orders$framed)

#Here I remove sequences containing dashes 
Annelida_coi5p_framed_2Orders <- Annelida_coi5p_framed_2Orders %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Annelida_coi5p_framed_2Orders$framed),
     main = paste("Annelida COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I remove sequences greater than 1000bp

Annelida_coi5p_framed_2Orders <- Annelida_coi5p_framed_2Orders %>%
  filter(nchar(framed) < 1000)


###Here I make 2 FASTA files. One with 10% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_2Orders)/2)
df1 = Annelida_coi5p_framed_2Orders[1:index,]
df2 = Annelida_coi5p_framed_2Orders[(index +1) : nrow(Annelida_coi5p_framed_2Orders),]

#In this case df1 has 10230 observations and df2 has 10231 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_1 <- character(nrow(df1) * 2)
Annelida_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Annelida_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Annelida_fasta_1, "Annelida_1.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2 <- character(nrow(df2) * 2)
Annelida_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Annelida_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2, "Annelida_clean.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Annelida_error_10per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Annelida_clean.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('A', '0001', df1_error$sequence)
df1_error$sequence <- gsub('C', '0010', df1_error$sequence)
df1_error$sequence <- gsub('G', '0100', df1_error$sequence)
df1_error$sequence <- gsub('T', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('A', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '1000', df2_clean$sequence)

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

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh$sequence)) #2016

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2532

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2018:2533)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Annelida_2Orders_oh_10per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 5% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_2Orders)/2)
df1 = Annelida_coi5p_framed_2Orders[1:index,]
df2 = Annelida_coi5p_framed_2Orders[(index +1) : nrow(Annelida_coi5p_framed_2Orders),]

#In this case df1 has 10230 observations and df2 has 10231 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_1 <- character(nrow(df1) * 2)
Annelida_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Annelida_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Annelida_fasta_1, "Annelida_major5.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2 <- character(nrow(df2) * 2)
Annelida_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Annelida_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2, "Annelida_clean_5per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Annelida_error_5per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Annelida_clean_5per.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('A', '0001', df1_error$sequence)
df1_error$sequence <- gsub('C', '0010', df1_error$sequence)
df1_error$sequence <- gsub('G', '0100', df1_error$sequence)
df1_error$sequence <- gsub('T', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('A', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '1000', df2_clean$sequence)

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

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh$sequence)) #2016

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2532

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2018:2533)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Annelida_2Orders_oh_5per.csv", row.names = FALSE)




#Here I make 2 FASTA files. One with 2% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_2Orders)/2)
df1 = Annelida_coi5p_framed_2Orders[1:index,]
df2 = Annelida_coi5p_framed_2Orders[(index +1) : nrow(Annelida_coi5p_framed_2Orders),]

#In this case df1 has 10230 observations and df2 has 10231 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_1 <- character(nrow(df1) * 2)
Annelida_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Annelida_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Annelida_fasta_1, "Annelida_major2.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2 <- character(nrow(df2) * 2)
Annelida_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Annelida_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2, "Annelida_clean_2per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Annelida_error_2per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Annelida_clean_2per.fasta")
seq_name = names(fastaFile2)
sequence = paste(fastaFile2)
df2_clean <- data.frame(seq_name, sequence)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error['Label']='1'

df2_clean['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error$sequence <- gsub('A', '0001', df1_error$sequence)
df1_error$sequence <- gsub('C', '0010', df1_error$sequence)
df1_error$sequence <- gsub('G', '0100', df1_error$sequence)
df1_error$sequence <- gsub('T', '1000', df1_error$sequence)

df2_clean$sequence <- gsub('A', '0001', df2_clean$sequence)
df2_clean$sequence <- gsub('C', '0010', df2_clean$sequence)
df2_clean$sequence <- gsub('G', '0100', df2_clean$sequence)
df2_clean$sequence <- gsub('T', '1000', df2_clean$sequence)

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

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh$sequence)) #2016

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2532

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2018:2533)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Annelida_2MajOrder_oh_2per.csv", row.names = FALSE)



######----- Applying the coil pipeline WITH ANNELIDA DATA CONTAINING 3 MINOR ORDERS----

annelida_coi5p_pipe_3Orders <- lapply(1:length(Annelida_3Orders$processid), function(i){
  coi5p_pipe(Annelida_3Orders$nucleotides[i],
             name = Annelida_3Orders$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(annelida_coi5p_pipe_3Orders)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Annelida_coi5p_df_3Orders <- flatten_coi5p(annelida_coi5p_pipe_3Orders)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Annelida_coi5p_df_3Orders$framed),
     main = paste("Annelida COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES ----

#Here I extract the framed sequences into their own dataframe
Annelida_coi5p_framed_3Orders <- subset(Annelida_coi5p_df_3Orders, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Annelida_coi5p_framed_3Orders$framed <- toupper(Annelida_coi5p_framed_3Orders$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters of all the sequences even if they do not have a dash.
Annelida_coi5p_framed_3Orders$framed <- gsub("^.{0,50}", "", Annelida_coi5p_framed_3Orders$framed)

#Here I remove sequences containing dashes 
Annelida_coi5p_framed_3Orders <- Annelida_coi5p_framed_3Orders %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Annelida_coi5p_framed_3Orders$framed),
     main = paste("Annelida COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I remove sequences greater than 1000bp

Annelida_coi5p_framed_3Orders <- Annelida_coi5p_framed_3Orders %>%
  filter(nchar(framed) < 1000)


#Here I make 2 FASTA files. One with 10% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor10 and df2_minor10 in case of nrow(df) = even and df1_minor10 will have 1 row less than df2_minor10 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_3Orders)/2)
df1_minor10 = Annelida_coi5p_framed_3Orders[1:index,]
df2_minor10 = Annelida_coi5p_framed_3Orders[(index +1) : nrow(Annelida_coi5p_framed_3Orders),]

#In this case df1 has 2141 observations and df2 has 2141 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_minor10 <- character(nrow(df1_minor10) * 2)
Annelida_fasta_minor10[c(TRUE, FALSE)] <- paste0(">", df1_minor10$name)
Annelida_fasta_minor10[c(FALSE, TRUE)] <- df1_minor10$framed

#Then write using writeLines:

writeLines(Annelida_fasta_minor10, "Annelida_minor10.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2_minor10 <- character(nrow(df2_minor10) * 2)
Annelida_fasta_2_minor10[c(TRUE, FALSE)] <- paste0(">", df2_minor10$name)
Annelida_fasta_2_minor10[c(FALSE, TRUE)] <- df2_minor10$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2_minor10, "Annelida_clean_minor10.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn10 <- readDNAStringSet("Annelida_error10_minor10.fasta")
seq_name_mn10 = names(fastaFile1_mn10)
sequence_mn10 = paste(fastaFile1_mn10)
df1_error_mn10 <- data.frame(seq_name_mn10, sequence_mn10)

fastaFile2_mn10 <- readDNAStringSet("Annelida_clean_minor10.fasta")
seq_name_mn10 = names(fastaFile2_mn10)
sequence_mn10 = paste(fastaFile2_mn10)
df2_clean_mn10 <- data.frame(seq_name_mn10, sequence_mn10)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error_mn10['Label']='1'

df2_clean_mn10['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error_mn10$sequence <- gsub('A', '0001', df1_error_mn10$sequence)
df1_error_mn10$sequence <- gsub('C', '0010', df1_error_mn10$sequence)
df1_error_mn10$sequence <- gsub('G', '0100', df1_error_mn10$sequence)
df1_error_mn10$sequence <- gsub('T', '1000', df1_error_mn10$sequence)

df2_clean_mn10$sequence <- gsub('A', '0001', df2_clean_mn10$sequence)
df2_clean_mn10$sequence <- gsub('C', '0010', df2_clean_mn10$sequence)
df2_clean_mn10$sequence <- gsub('G', '0100', df2_clean_mn10$sequence)
df2_clean_mn10$sequence <- gsub('T', '1000', df2_clean_mn10$sequence)

#Combine the 2 dataframes

df3_oh_mn10 <- rbind(df1_error_mn10, df2_clean_mn10)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh_mn10))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh_mn10 <- df3_oh_mn10[rows, ]

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh_mn10$sequence)) #2004

#maximum sequence length
max(nchar(df3_shuffled_oh_mn10$sequence)) #2592

#Here I remove sequences less than 2016 (since that was the minimum for the 2 major orders)

df3_shuffled_oh_mn10 <- df3_shuffled_oh_mn10 %>%
  filter(nchar(sequence) > 2016) #This only removed one row (the row containing sequence length of 2004)

min(nchar(df3_shuffled_oh_mn10$sequence)) #now 2084

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh_mn10$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh_mn10$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh_mn10 <- cbind(df3_shuffled_oh_mn10, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh_mn10$sequence_mn10 <- NULL

#Remove the seq_name column

drop <- c("seq_name_mn10")
df3_shuffled_split_oh_mn10 <- df3_shuffled_split_oh_mn10[ , !(names(df3_shuffled_split_oh_mn10) %in% drop)]

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh_mn10 <- df3_shuffled_split_oh_mn10[, -c(2018:2593)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn10, file = "Annelida_3MinorOrders_oh_10per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 5% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor2 and df2_minor2 in case of nrow(df) = even and df1_minor2 will have 1 row less than df2_minor2 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_3Orders)/2)
df1_minor5 = Annelida_coi5p_framed_3Orders[1:index,]
df2_minor5 = Annelida_coi5p_framed_3Orders[(index +1) : nrow(Annelida_coi5p_framed_3Orders),]

#In this case df1 has 2141 observations and df2 has 2141 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_minor5 <- character(nrow(df1_minor5) * 2)
Annelida_fasta_minor5[c(TRUE, FALSE)] <- paste0(">", df1_minor5$name)
Annelida_fasta_minor5[c(FALSE, TRUE)] <- df1_minor5$framed

#Then write using writeLines:

writeLines(Annelida_fasta_minor5, "Annelida_minor5.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2_minor5 <- character(nrow(df2_minor5) * 2)
Annelida_fasta_2_minor5[c(TRUE, FALSE)] <- paste0(">", df2_minor5$name)
Annelida_fasta_2_minor5[c(FALSE, TRUE)] <- df2_minor5$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2_minor5, "Annelida_clean_minor5.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn5 <- readDNAStringSet("Annelida_error_minor5.fasta")
seq_name_mn5 = names(fastaFile1_mn5)
sequence_mn5 = paste(fastaFile1_mn5)
df1_error_mn5 <- data.frame(seq_name_mn5, sequence_mn5)

fastaFile2_mn5 <- readDNAStringSet("Annelida_clean_minor5.fasta")
seq_name_mn5 = names(fastaFile2_mn5)
sequence_mn5 = paste(fastaFile2_mn5)
df2_clean_mn5 <- data.frame(seq_name_mn5, sequence_mn5)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error_mn5['Label']='1'

df2_clean_mn5['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error_mn5$sequence <- gsub('A', '0001', df1_error_mn5$sequence)
df1_error_mn5$sequence <- gsub('C', '0010', df1_error_mn5$sequence)
df1_error_mn5$sequence <- gsub('G', '0100', df1_error_mn5$sequence)
df1_error_mn5$sequence <- gsub('T', '1000', df1_error_mn5$sequence)

df2_clean_mn5$sequence <- gsub('A', '0001', df2_clean_mn5$sequence)
df2_clean_mn5$sequence <- gsub('C', '0010', df2_clean_mn5$sequence)
df2_clean_mn5$sequence <- gsub('G', '0100', df2_clean_mn5$sequence)
df2_clean_mn5$sequence <- gsub('T', '1000', df2_clean_mn5$sequence)

#Combine the 2 dataframes

df3_oh_mn5 <- rbind(df1_error_mn5, df2_clean_mn5)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh_mn5))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh_mn5 <- df3_oh_mn5[rows, ]

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh_mn5$sequence)) #2004

#maximum sequence length
max(nchar(df3_shuffled_oh_mn5$sequence)) #2592

#Here I remove sequences less than 2016 (since that was the minimum for the 2 major orders)

df3_shuffled_oh_mn5 <- df3_shuffled_oh_mn5 %>%
  filter(nchar(sequence) > 2016) #This only removed one row (the row containing sequence length of 2004)

min(nchar(df3_shuffled_oh_mn5$sequence)) #now 2084

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh_mn5$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh_mn5$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh_mn5 <- cbind(df3_shuffled_oh_mn5, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh_mn5$sequence_mn5 <- NULL
df3_shuffled_split_oh_mn5$sequence <- NULL

#Remove the seq_name column

drop <- c("seq_name_mn5")
df3_shuffled_split_oh_mn5 <- df3_shuffled_split_oh_mn5[ , !(names(df3_shuffled_split_oh_mn5) %in% drop)]

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh_mn5 <- df3_shuffled_split_oh_mn5[, -c(2018:2593)] 


#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn5, file = "Annelida_3MinorOrders_oh_5per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 2% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor2 and df2_minor2 in case of nrow(df) = even and df1_minor2 will have 1 row less than df2_minor2 in case of nrow(df) = odd

index = floor(nrow(Annelida_coi5p_framed_3Orders)/2)
df1_minor2 = Annelida_coi5p_framed_3Orders[1:index,]
df2_minor2 = Annelida_coi5p_framed_3Orders[(index +1) : nrow(Annelida_coi5p_framed_3Orders),]

#In this case df1 has 2141 observations and df2 has 2141 observations  

##df1 to FASTA (this wil be the file I add errors to):

Annelida_fasta_minor2 <- character(nrow(df1_minor2) * 2)
Annelida_fasta_minor2[c(TRUE, FALSE)] <- paste0(">", df1_minor2$name)
Annelida_fasta_minor2[c(FALSE, TRUE)] <- df1_minor2$framed

#Then write using writeLines:

writeLines(Annelida_fasta_minor2, "Annelida_minor2.fasta") 

##df2 to FASTA (this wil be the error free file):

Annelida_fasta_2_minor2 <- character(nrow(df2_minor2) * 2)
Annelida_fasta_2_minor2[c(TRUE, FALSE)] <- paste0(">", df2_minor2$name)
Annelida_fasta_2_minor2[c(FALSE, TRUE)] <- df2_minor2$framed

#Then write using writeLines:

writeLines(Annelida_fasta_2_minor2, "Annelida_clean_minor2.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn2 <- readDNAStringSet("Annelida_error_minor2.fasta")
seq_name_mn2 = names(fastaFile1_mn2)
sequence_mn2 = paste(fastaFile1_mn2)
df1_error_mn2 <- data.frame(seq_name_mn2, sequence_mn2)

fastaFile2_mn2 <- readDNAStringSet("Annelida_clean_minor2.fasta")
seq_name_mn2 = names(fastaFile2_mn2)
sequence_mn2 = paste(fastaFile2_mn2)
df2_clean_mn2 <- data.frame(seq_name_mn2, sequence_mn2)

#We will add a column to both dataframes with a label. 1 for error and 2 for clean.

df1_error_mn2['Label']='1'

df2_clean_mn2['Label']='2'

#Now we do one-hot encoding

#We will replace ACGT with on-hot encoded values in both dataframes

df1_error_mn2$sequence <- gsub('A', '0001', df1_error_mn2$sequence)
df1_error_mn2$sequence <- gsub('C', '0010', df1_error_mn2$sequence)
df1_error_mn2$sequence <- gsub('G', '0100', df1_error_mn2$sequence)
df1_error_mn2$sequence <- gsub('T', '1000', df1_error_mn2$sequence)

df2_clean_mn2$sequence <- gsub('A', '0001', df2_clean_mn2$sequence)
df2_clean_mn2$sequence <- gsub('C', '0010', df2_clean_mn2$sequence)
df2_clean_mn2$sequence <- gsub('G', '0100', df2_clean_mn2$sequence)
df2_clean_mn2$sequence <- gsub('T', '1000', df2_clean_mn2$sequence)

#Combine the 2 dataframes

df3_oh_mn2 <- rbind(df1_error_mn2, df2_clean_mn2)

#Shuffle the dataframe by rows:

#First, you set a random seed so that your work is reproducible and you get the same random split each time you run your script
set.seed(42)

#Next, you use the sample() function to shuffle the row indices of the dataframe(df3). You can later use these indices to reorder the dataset.
rows <- sample(nrow(df3_oh_mn2))

#Finally, you can use this random vector to reorder the dataframe:
df3_shuffled_oh_mn2 <- df3_oh_mn2[rows, ]

#Now we should determine the minimum (and maximum) sequence lengths to split based on the size of the minimum sequence length (i.e. remove all end characters of the other sequences that are longer than the minimum sequence)

#minimum sequence length
min(nchar(df3_shuffled_oh_mn2$sequence)) #2004

#maximum sequence length
max(nchar(df3_shuffled_oh_mn2$sequence)) #2592

#Here I remove sequences less than 2016 (since that was the minimum for the 2 major orders)

df3_shuffled_oh_mn2 <- df3_shuffled_oh_mn2 %>%
  filter(nchar(sequence) > 2016) #This only removed one row (the row containing sequence length of 2004)

min(nchar(df3_shuffled_oh_mn2$sequence)) #now 2084

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh_mn2$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh_mn2$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh_mn2 <- cbind(df3_shuffled_oh_mn2, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh_mn2$sequence_mn2 <- NULL
df3_shuffled_split_oh_mn2$sequence <- NULL

#Remove the seq_name column

drop <- c("seq_name_mn2")
df3_shuffled_split_oh_mn2 <- df3_shuffled_split_oh_mn2[ , !(names(df3_shuffled_split_oh_mn2) %in% drop)]

#The minimum sequence length was 2016. So now, we just manually erase all columns after 2017 so that there is no column with NA.

df3_shuffled_split_oh_mn2 <- df3_shuffled_split_oh_mn2[, -c(2018:2593)] 


#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn2, file = "Annelida_3MinorOrder_oh_2per.csv", row.names = FALSE)
