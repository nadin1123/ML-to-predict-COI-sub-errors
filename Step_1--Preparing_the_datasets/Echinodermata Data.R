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



######----- ACQUIRING ECHINODERMATA DATA FROM BOLD---- 

Echinodermata <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Echinodermata&format=tsv")
#Result: 35268 observations with 80 variables
#Read on March 9, 2022

Echinodermata <- read.delim("Echinodermata Data from BOLD.txt")


######----- FILTERING ECHINODERMATA COI-5P SEQUENCES FROM BOLD---- 

#Criteria:
#1)Folmer region of the sequence must be greater that 600 base pairs
#2)The taxonomy of the sequence must be known at least to the genus level
#3)There are no missing base pairs in the sequence


Echinodermata.COI <- Echinodermata %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[BDEFHIJKLMNOPQRSUVWXYZ]", TRUE)) %>% #Criteria 3 + removing any nucleotide sequences that contain any letters besides ACGT
  filter(!is.na(genus_name)) %>% #Criteria 2
  filter(nchar(nucleotides) > 600) #Criteria 1 (Maybe, as this is the whole COI gene not simply the Folmer region)


#Histogram of nucleotide sequence lengths
hist(nchar(Echinodermata.COI$nucleotides),
     main = paste("Echinodermata COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- BAR PLOTS FOR # OF DIFFERENT TAXANOMIC GROUPS----

#Class

barplot(table(Echinodermata.COI$class_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Classes of Echinodermata")
#Order

barplot(table(Echinodermata.COI$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Echinodermata")

#Family

barplot(table(Echinodermata.COI$family_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Families of Echinodermata",
        las=2,
        cex.names = 0.6,
        srt = 60)

#Better looking ggplot barplots (these contain N/A taxanomic groups):

#Class

ggplot(data = Echinodermata.COI) + 
  geom_bar(mapping = aes(x = class_name), stat = "count", fill = "turquoise") +
  labs(title = "Classes of Echinodermata", x = "Taxanomic Group", y = "Counts")

#Order

#https://jacintak.github.io/post/2021-08-01-rphylopic/#:~:text=Use%20add_phylopic%20to%20add%20the,alpha%20to%20control%20the%20transparency.
#https://cran.r-project.org/web/packages/rphylopic/readme/README.html
echin <- name_search(text = "Echinodermata", options = "namebankID")[[1]] # find names
echin_id_all <- name_images(uuid = echin$uid[1])
echin_id <- name_images(uuid = echin$uid[1])$same[[1]]$uid
echin_pic <- image_data(echin_id, size = 512)[[1]]

ggplot(data = Echinodermata.COI) + 
  geom_bar(mapping = aes(x = order_name), stat = "count", fill = "turquoise") +
  labs(title = "Orders of Echinodermata from BOLD", x = "Taxonomic Group", y = "Counts") + coord_flip() +
  theme(axis.text=element_text(size=15),
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title=element_text(size=16,face="bold")) + theme_bw() + add_phylopic(echin_pic, ysize = 10, x = 3000, y = 20)

#Family

ggplot(data = Echinodermata.COI) + 
  geom_bar(mapping = aes(x = family_name), stat = "count", fill = "turquoise") +
  labs(title = "Families of Echinodermata", x = "Taxanomic Group", y = "Counts") +
  theme(axis.text.x = element_text(angle = 90))



######----- REMOVING ALL ORDERS EXCEPT THE 3 LARGEST----

Echinodermata_Order_Counts <- as.data.frame(table(Echinodermata.COI$order_name))

Echinodermata_Order_Counts <- Echinodermata_Order_Counts[order(-Echinodermata_Order_Counts$Freq), ]
Echinodermata_Order_Counts

#The three largest orders are:
#Valvatida 3062, Comatulida 2459 & Amphilepidida 2280

Echinodermata_3MajOrders <- Echinodermata.COI %>%
  filter(str_detect(order_name, "Forcipulatida", TRUE)) %>%
  filter(str_detect(order_name, "Camarodonta", TRUE)) %>%
  filter(str_detect(order_name, "Ophiacanthida", TRUE)) %>%
  filter(str_detect(order_name, "Dendrochirotida", TRUE)) %>%
  filter(str_detect(order_name, "Ophiurida", TRUE)) %>%
  filter(str_detect(order_name, "Holothuriida", TRUE)) %>%
  filter(str_detect(order_name, "Clypeasteroida", TRUE)) %>%
  filter(str_detect(order_name, "Paxillosida", TRUE)) %>%
  filter(str_detect(order_name, "Synallactida", TRUE)) %>%
  filter(str_detect(order_name, "Spinulosida", TRUE)) %>%
  filter(str_detect(order_name, "Arbacioida", TRUE)) %>%
  filter(str_detect(order_name, "Euryalida", TRUE)) %>%
  filter(str_detect(order_name, "Spatangoida", TRUE)) %>%
  filter(str_detect(order_name, "Cidaroida", TRUE)) %>%
  filter(str_detect(order_name, "Elasipodida", TRUE)) %>%
  filter(str_detect(order_name, "Velatida", TRUE)) %>%
  filter(str_detect(order_name, "Echinothurioida", TRUE)) %>%
  filter(str_detect(order_name, "Notomyotida", TRUE)) %>%
  filter(str_detect(order_name, "Apodida", TRUE)) %>%
  filter(str_detect(order_name, "Diadematoida", TRUE)) %>%
  filter(str_detect(order_name, "Molpadida", TRUE)) %>%
  filter(str_detect(order_name, "Pedinoida", TRUE)) %>%
  filter(str_detect(order_name, "Brisingida", TRUE)) %>%
  filter(str_detect(order_name, "Isocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Holasteroida", TRUE)) %>%
  filter(str_detect(order_name, "Hyocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Cyrtocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Aspidodiadematoida", TRUE)) %>%
  filter(str_detect(order_name, "Cassiduloida", TRUE)) %>%
  filter(str_detect(order_name, "Echinoneoida", TRUE)) %>%
  filter(str_detect(order_name, "Stomopneustoida", TRUE)) %>%
  filter(str_detect(order_name, "Echinolampadoida", TRUE)) %>%
  filter(str_detect(order_name, "Micropygoida", TRUE)) %>%
  filter(str_detect(order_name, "Persiculida", TRUE)) %>%
  filter(str_detect(order_name, "Phymosomatoida", TRUE)) 
  

table(Echinodermata_3MajOrders$order_name)

barplot(table(Echinodermata_3MajOrders$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Echinodermata")



######----- REMOVING ALL ORDERS EXCEPT FOR 3 MINOR----

#The 3 orders I will be isolating are Forcipulatida 1537, Camarodonta 1432, and Ophiacanthida 1303

Echinodermata_3MinOrders <- Echinodermata.COI %>%
  filter(str_detect(order_name, "Valvatida", TRUE)) %>%
  filter(str_detect(order_name, "Comatulida", TRUE)) %>%
  filter(str_detect(order_name, "Amphilepidida", TRUE)) %>%
  filter(str_detect(order_name, "Dendrochirotida", TRUE)) %>%
  filter(str_detect(order_name, "Ophiurida", TRUE)) %>%
  filter(str_detect(order_name, "Holothuriida", TRUE)) %>%
  filter(str_detect(order_name, "Clypeasteroida", TRUE)) %>%
  filter(str_detect(order_name, "Paxillosida", TRUE)) %>%
  filter(str_detect(order_name, "Synallactida", TRUE)) %>%
  filter(str_detect(order_name, "Spinulosida", TRUE)) %>%
  filter(str_detect(order_name, "Arbacioida", TRUE)) %>%
  filter(str_detect(order_name, "Euryalida", TRUE)) %>%
  filter(str_detect(order_name, "Spatangoida", TRUE)) %>%
  filter(str_detect(order_name, "Cidaroida", TRUE)) %>%
  filter(str_detect(order_name, "Elasipodida", TRUE)) %>%
  filter(str_detect(order_name, "Velatida", TRUE)) %>%
  filter(str_detect(order_name, "Echinothurioida", TRUE)) %>%
  filter(str_detect(order_name, "Notomyotida", TRUE)) %>%
  filter(str_detect(order_name, "Apodida", TRUE)) %>%
  filter(str_detect(order_name, "Diadematoida", TRUE)) %>%
  filter(str_detect(order_name, "Molpadida", TRUE)) %>%
  filter(str_detect(order_name, "Pedinoida", TRUE)) %>%
  filter(str_detect(order_name, "Brisingida", TRUE)) %>%
  filter(str_detect(order_name, "Isocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Holasteroida", TRUE)) %>%
  filter(str_detect(order_name, "Hyocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Cyrtocrinida", TRUE)) %>%
  filter(str_detect(order_name, "Aspidodiadematoida", TRUE)) %>%
  filter(str_detect(order_name, "Cassiduloida", TRUE)) %>%
  filter(str_detect(order_name, "Echinoneoida", TRUE)) %>%
  filter(str_detect(order_name, "Stomopneustoida", TRUE)) %>%
  filter(str_detect(order_name, "Echinolampadoida", TRUE)) %>%
  filter(str_detect(order_name, "Micropygoida", TRUE)) %>%
  filter(str_detect(order_name, "Persiculida", TRUE)) %>%
  filter(str_detect(order_name, "Phymosomatoida", TRUE))  

table(Echinodermata_3MinOrders$order_name)

barplot(table(Echinodermata_3MinOrders$order_name),
        xlab = "Taxanomic Group",
        ylab = "Counts",
        main = "Orders of Echinodermata")


######----- Applying the coil pipeline WITH ECHINODERMATA DATA CONTAINING 3 MAJOR ORDERS----

echinodermata_coi5p_pipe_3MajOrders <- lapply(1:length(Echinodermata_3MajOrders$processid), function(i){
  coi5p_pipe(Echinodermata_3MajOrders$nucleotides[i],
             name = Echinodermata_3MajOrders$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(echinodermata_coi5p_pipe_3MajOrders)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Echinodermata_coi5p_df_3MajOrders <- flatten_coi5p(echinodermata_coi5p_pipe_3MajOrders)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Echinodermata_coi5p_df_3MajOrders$framed),
     main = paste("Echinodermata COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- Applying the coil pipeline WITH ECHINODERMATA DATA CONTAINING 3 MINOR ORDERS----

echinodermata_coi5p_pipe_3MinOrders <- lapply(1:length(Echinodermata_3MinOrders$processid), function(i){
  coi5p_pipe(Echinodermata_3MinOrders$nucleotides[i],
             name = Echinodermata_3MinOrders$processid[i],
             trans_table = 5)
})

#The coi5p objects are in a list
class(echinodermata_coi5p_pipe_3MinOrders)

#This next function (flatten_coi5p) flattens our list of coi5p output objects into a dataframe
Echinodermata_coi5p_df_3MinOrders <- flatten_coi5p(echinodermata_coi5p_pipe_3MinOrders)

#Histogram of framed nucleotide sequence lengths
hist(nchar(Echinodermata_coi5p_df_3MinOrders$framed),
     main = paste("Echinodermata COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")


######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES FOR 3 MAJOR ORDERS ----

#Here I extract the framed sequences into their own dataframe
Echinodermata_coi5p_framed_MajOrders <- subset(Echinodermata_coi5p_df_3MajOrders, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Echinodermata_coi5p_framed_MajOrders$framed <- toupper(Echinodermata_coi5p_framed_MajOrders$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters of all the sequences even if they do not have a dash.
Echinodermata_coi5p_framed_MajOrders$framed <- gsub("^.{0,50}", "", Echinodermata_coi5p_framed_MajOrders$framed)

#Here I remove sequences containing dashes 
Echinodermata_coi5p_framed_MajOrders <- Echinodermata_coi5p_framed_MajOrders %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Echinodermata_coi5p_framed_MajOrders$framed),
     main = paste("Echinodermata COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I remove sequences greater than 1000bp

Echinodermata_coi5p_framed_MajOrders <- Echinodermata_coi5p_framed_MajOrders %>%
  filter(nchar(framed) < 1000)


###Here I make 2 FASTA files. One with 10% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MajOrders)/2)
df1 = Echinodermata_coi5p_framed_MajOrders[1:index,]
df2 = Echinodermata_coi5p_framed_MajOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MajOrders),]

#In this case df1 has 1779 observations and df2 has 1779 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_1 <- character(nrow(df1) * 2)
Echinodermata_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Echinodermata_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_1, "Echinodermata_1.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermatafasta_2 <- character(nrow(df2) * 2)
Echinodermata_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Echinodermata_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2, "Echinodermata_clean_major_10per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Echinodermata_error_major_10per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Echinodermata_clean_major_10per.fasta")
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
min(nchar(df3_shuffled_oh$sequence)) #2080

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2828

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2082:2829)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Echinodermata_3MajOrders_oh_10per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 5% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MajOrders)/2)
df1 = Echinodermata_coi5p_framed_MajOrders[1:index,]
df2 = Echinodermata_coi5p_framed_MajOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MajOrders),]

#In this case df1 has 1779 observations and df2 has 1779 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_1 <- character(nrow(df1) * 2)
Echinodermata_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Echinodermata_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_1, "Echinodermata_major5.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermata_fasta_2 <- character(nrow(df2) * 2)
Echinodermata_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Echinodermata_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2, "Echinodermata_clean_major_5per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Echinodermata_error_major_5per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Echinodermata_clean_major_5per.fasta")
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
min(nchar(df3_shuffled_oh$sequence)) #2080

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2828

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2082:2829)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Echinodermata_3MajOrders_oh_5per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 2% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1 and df2 in case of nrow(df) = even and df1 will have 1 row less than df2 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MajOrders)/2)
df1 = Echinodermata_coi5p_framed_MajOrders[1:index,]
df2 = Echinodermata_coi5p_framed_MajOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MajOrders),]

#In this case df1 has 1779 observations and df2 has 1779 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_1 <- character(nrow(df1) * 2)
Echinodermata_fasta_1[c(TRUE, FALSE)] <- paste0(">", df1$name)
Echinodermata_fasta_1[c(FALSE, TRUE)] <- df1$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_1, "Echinodermata_major2.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermata_fasta_2 <- character(nrow(df2) * 2)
Echinodermata_fasta_2[c(TRUE, FALSE)] <- paste0(">", df2$name)
Echinodermata_fasta_2[c(FALSE, TRUE)] <- df2$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2, "Echinodermata_clean_major_2per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1 <- readDNAStringSet("Echinodermata_error_major_2per.fasta")
seq_name = names(fastaFile1)
sequence = paste(fastaFile1)
df1_error <- data.frame(seq_name, sequence)

fastaFile2 <- readDNAStringSet("Echinodermata_clean_major_2per.fasta")
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
min(nchar(df3_shuffled_oh$sequence)) #2080

#maximum sequence length
max(nchar(df3_shuffled_oh$sequence)) #2828

#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh <- cbind(df3_shuffled_oh, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh$sequence <- NULL

#The minimum sequence length was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh <- df3_shuffled_split_oh[, -c(2082:2829)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh, file = "Echinodermata_3MajOrders_oh_2per.csv", row.names = FALSE)





######----- CREATING A FASTA FILE OF THE NUCLEOTIDE SEQUENCES FOR 3 MINOR ORDERS ----

#Here I extract the framed sequences into their own dataframe
Echinodermata_coi5p_framed_MinOrders <- subset(Echinodermata_coi5p_df_3MinOrders, select=c("name", "framed"))

#FASTA file DNA sequences have to be uppercase. So here I convert all sequences from lowercase to uppercase.
Echinodermata_coi5p_framed_MinOrders$framed <- toupper(Echinodermata_coi5p_framed_MinOrders$framed)

#Here I trim the sequences so that there are no dashes (-) at the start.I remove the first 50 characters of all the sequences even if they do not have a dash.
Echinodermata_coi5p_framed_MinOrders$framed <- gsub("^.{0,50}", "", Echinodermata_coi5p_framed_MinOrders$framed)

#Here I remove sequences containing dashes 
Echinodermata_coi5p_framed_MinOrders <- Echinodermata_coi5p_framed_MinOrders %>%
  filter(str_detect(framed, "[-]", TRUE))

#Histogram of framed nucleotide sequence lengths
hist(nchar(Echinodermata_coi5p_framed_MinOrders$framed),
     main = paste("Echinodermata COI-5P Sequence Lengths from BOLD"),
     xlab = "Sequence Lengths (bp)")

#Here I remove sequences greater than 1000bp

Echinodermata_coi5p_framed_MinOrders <- Echinodermata_coi5p_framed_MinOrders %>%
  filter(nchar(framed) < 1000)



#Here I make 2 FASTA files. One with 10% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor10 and df2_minor10 in case of nrow(df) = even and df1_minor10 will have 1 row less than df2_minor10 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MinOrders)/2)
df1_minor10 = Echinodermata_coi5p_framed_MinOrders[1:index,]
df2_minor10 = Echinodermata_coi5p_framed_MinOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MinOrders),]

#In this case df1 has 1294 observations and df2 has 1294 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_minor10 <- character(nrow(df1_minor10) * 2)
Echinodermata_fasta_minor10[c(TRUE, FALSE)] <- paste0(">", df1_minor10$name)
Echinodermata_fasta_minor10[c(FALSE, TRUE)] <- df1_minor10$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_minor10, "Echinodermata_minor10.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermata_fasta_2_minor10 <- character(nrow(df2_minor10) * 2)
Echinodermata_fasta_2_minor10[c(TRUE, FALSE)] <- paste0(">", df2_minor10$name)
Echinodermata_fasta_2_minor10[c(FALSE, TRUE)] <- df2_minor10$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2_minor10, "Echinodermata_clean_minor_10per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn10 <- readDNAStringSet("Echinodermata_error_minor_10per.fasta")
seq_name_mn10 = names(fastaFile1_mn10)
sequence_mn10 = paste(fastaFile1_mn10)
df1_error_mn10 <- data.frame(seq_name_mn10, sequence_mn10)

fastaFile2_mn10 <- readDNAStringSet("Echinodermata_clean_minor_10per.fasta")
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
min(nchar(df3_shuffled_oh_mn10$sequence)) #2208

#maximum sequence length
max(nchar(df3_shuffled_oh_mn10$sequence)) #2540


#I will divide each base into 1 column

tmp1 <- read.fwf(
  textConnection(df3_shuffled_oh_mn10$sequence),
  widths = rep(1, ceiling(max(nchar(df3_shuffled_oh_mn10$sequence) / 1))),
  stringsAsFactors = FALSE)

df3_shuffled_split_oh_mn10 <- cbind(df3_shuffled_oh_mn10, tmp1)

#Remove 'sequence' label 

df3_shuffled_split_oh_mn10$sequence_mn10 <- NULL
df3_shuffled_split_oh_mn10$sequence <- NULL

#Remove the seq_name column

drop <- c("seq_name_mn10")
df3_shuffled_split_oh_mn10 <- df3_shuffled_split_oh_mn10[ , !(names(df3_shuffled_split_oh_mn10) %in% drop)]

#The minimum sequence length is 2208. The minimum sequence length for the major orders was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh_mn10 <- df3_shuffled_split_oh_mn10[, -c(2082:2541)]  

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn10, file = "Echinodermata_3MinOrders_oh_10per.csv", row.names = FALSE)



#Here I make 2 FASTA files. One with 5% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor2 and df2_minor2 in case of nrow(df) = even and df1_minor2 will have 1 row less than df2_minor2 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MinOrders)/2)
df1_minor5 = Echinodermata_coi5p_framed_MinOrders[1:index,]
df2_minor5 = Echinodermata_coi5p_framed_MinOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MinOrders),]

#In this case df1 has 1294 observations and df2 has 1294 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_minor5 <- character(nrow(df1_minor5) * 2)
Echinodermata_fasta_minor5[c(TRUE, FALSE)] <- paste0(">", df1_minor5$name)
Echinodermata_fasta_minor5[c(FALSE, TRUE)] <- df1_minor5$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_minor5, "Echinodermata_minor5.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermata_fasta_2_minor5 <- character(nrow(df2_minor5) * 2)
Echinodermata_fasta_2_minor5[c(TRUE, FALSE)] <- paste0(">", df2_minor5$name)
Echinodermata_fasta_2_minor5[c(FALSE, TRUE)] <- df2_minor5$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2_minor5, "Echinodermata_clean_minor_5per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn5 <- readDNAStringSet("Echinodermata_error_minor_5per.fasta")
seq_name_mn5 = names(fastaFile1_mn5)
sequence_mn5 = paste(fastaFile1_mn5)
df1_error_mn5 <- data.frame(seq_name_mn5, sequence_mn5)

fastaFile2_mn5 <- readDNAStringSet("Echinodermata_clean_minor_5per.fasta")
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
min(nchar(df3_shuffled_oh_mn5$sequence)) #2208

#maximum sequence length
max(nchar(df3_shuffled_oh_mn5$sequence)) #2540

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

#The minimum sequence length is 2208. The minimum sequence length for the major orders was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh_mn5 <- df3_shuffled_split_oh_mn5[, -c(2082:2541)] 

#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn5, file = "Echinodermata_3MinOrders_oh_5per.csv", row.names = FALSE)


#Here I make 2 FASTA files. One with 2% substitution errors to be added (using error generator) and another with no substitution errors added

#First I divide the dataframe in half
#I will split the original df into two equal dataframe df1_minor2 and df2_minor2 in case of nrow(df) = even and df1_minor2 will have 1 row less than df2_minor2 in case of nrow(df) = odd

index = floor(nrow(Echinodermata_coi5p_framed_MinOrders)/2)
df1_minor2 = Echinodermata_coi5p_framed_MinOrders[1:index,]
df2_minor2 = Echinodermata_coi5p_framed_MinOrders[(index +1) : nrow(Echinodermata_coi5p_framed_MinOrders),]

#In this case df1 has 1294 observations and df2 has 1294 observations  

##df1 to FASTA (this wil be the file I add errors to):

Echinodermata_fasta_minor2 <- character(nrow(df1_minor2) * 2)
Echinodermata_fasta_minor2[c(TRUE, FALSE)] <- paste0(">", df1_minor2$name)
Echinodermata_fasta_minor2[c(FALSE, TRUE)] <- df1_minor2$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_minor2, "Echinodermata_minor2.fasta") 

##df2 to FASTA (this wil be the error free file):

Echinodermata_fasta_2_minor2 <- character(nrow(df2_minor2) * 2)
Echinodermata_fasta_2_minor2[c(TRUE, FALSE)] <- paste0(">", df2_minor2$name)
Echinodermata_fasta_2_minor2[c(FALSE, TRUE)] <- df2_minor2$framed

#Then write using writeLines:

writeLines(Echinodermata_fasta_2_minor2, "Echinodermata_clean_minor_2per.fasta")


######----- COMBINING THE 2 FASTA FILES AFTER ERROR GENERATION TO 1 OF THEM INTO A DATAFRAME & CSV FILE WITH LABELS----

#First we read the 2 FASTA files into a dataframe

fastaFile1_mn2 <- readDNAStringSet("Echinodermata_error_minor_2per.fasta")
seq_name_mn2 = names(fastaFile1_mn2)
sequence_mn2 = paste(fastaFile1_mn2)
df1_error_mn2 <- data.frame(seq_name_mn2, sequence_mn2)

fastaFile2_mn2 <- readDNAStringSet("Echinodermata_clean_minor_2per.fasta")
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
min(nchar(df3_shuffled_oh_mn2$sequence)) #2208

#maximum sequence length
max(nchar(df3_shuffled_oh_mn2$sequence)) #2540

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

#The minimum sequence length is 2208. The minimum sequence length for the major orders was 2080. So now, we just manually erase all columns after 2081 so that there is no column with NA.

df3_shuffled_split_oh_mn2 <- df3_shuffled_split_oh_mn2[, -c(2082:2541)] 


#Export DataFrame to CSV

write.csv(df3_shuffled_split_oh_mn2, file = "Echinodermata_3MinOrders_oh_2per.csv", row.names = FALSE)

