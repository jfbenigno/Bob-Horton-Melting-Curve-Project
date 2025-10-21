if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER", force = TRUE)
BiocManager::install("Biostrings", force = TRUE)

library(Biostrings)
library(DECIPHER)
mySeq <- ("gagttttatcgcttccatgacgcagaagttaacactttcggatatttctgatgagtcgaaaaattatcttgataaagcaggaattactactgcttgtttacgaattaaatcgaagtggactgctggcggaaaatgagaaaattcgacctatccttgcgcagctcgagaagctcttactttgcgacctttcgccatcaactaacgattctgtcaaaaactgacgcgttggatgaggagaagtggcttaatatgcttggcacgttcgtcaaggactggtttagatatgagtcacattttgttcatggtagagattctcttgttgacattttaaaagagcgtggattactatctgagtccgatgctgttcaaccactaataggtaagaaatcatgagtcaagttactgaacaatccgtacgtttccagaccgctttggcctctattaagctcattcaggcttctgccgttttggatttaaccgaagatgatttcgattttctgacgagtaacaaagtttggattgctactgaccgctctcgtgctcgtcgctgcgttgaggcttgcgtttatggtacgctggactttgtgggataccctcgctttcctgctcctgttgagtttattgctgccgtcattgcttatttatgttcatcccgtcaacattcaaacggcctgtctcatcatggaaggcgctgaatttacggaaaacattattaatggcgtcgagcgtccggttaaagccgctgaattgttcgcgtttaccttgcgtgtacgcgcaggaaacactgacgttcttactgacgcagaagaaaacgtgcgtcaaaaattacgtgcggaaggagtgatgtaatgtctaaaggtaaaaaacgttctggcgctcgccctggtcgtccgcagccgttgcgaggtactaaaggcaagcgtaaaggcgctcgtctttggtatgtaggtggtcaacaattttaattgcaggggcttcggccccttacttgaggataaattatgtctaatattcaaactggcgccgagcgtatgccgcatgacctttcccatcttggcttccttgctggtcagattggtcgtcttattaccatttcaactactccggttatcgctggcgactccttcgagatggacgccgttggcgctctccgtctttctccattgcgtcgtggccttgctattgactctactgtagacatttttactttttatgtccctcatcgtcacgtttatggtgaacagtggattaagttcatgaaggatggtgttaatgccactcctctcccgactgttaacactactggttatattgaccatgccgcttttcttggcacgattaaccctgataccaataaaatccctaagcatttgtttcagggttatttgaatatctataacaactattttaaagcgccgtggatgcctgaccgtaccgaggctaaccctaatgagcttaatcaagatgatgctcgttatggtttccgttgctgccatctcaaaacatttggactgctccgcttcctcctgagactgagctttctcgccaaatgacgacttctaccacatctattgacattatgggtctgcaagctgcttatgctaatttgctactgaccaagaacgtgattacttcatgcagcgttaccatgatgttatttcttcatttggaggtaaaacctcttatgacgctgacaaccgtcctttacttgtcatgcgctctaatctctgggcatctggctatgatgttgatggaactgaccaaacgtcgttaggccagttttctggtcgtgttcaacagacctataaacattctgtgccgcgtttctttgttcctgagcatggcactatgtttactcttgcgcttgttcgttttccgcctactgcgactaaagagatttcagtaccttaacgctaaaggtgctttgacttataccgatattgctggcgaccctgttttgtatggcaacttgccgccgcgtgaaatttctatgaaggatgttttccgttctggtgattcgtctaagaagtttaagattgctgagggtcagtggtatcgttatgcgccttcgtatgtttctcctgcttatcaccttcttgaaggcttcccat tcattcaggaaccgccttctggtgatttgcaagaacgcgtacttattcgccaccatgattatgaccagtgtttccagtccgttcagttgttgcagtggaatagtcaggttaaatttaatgtgaccgtttatcgcaatctgccgaccactcgcgattcaatcatgacttcgtgataaaagattgagtgtgaggttataacgccgaagcggtaaaaattttaatttttgccgctgaggggttgaccaagcgaagcgcggtaggttttctgcttaggagtttaatcatgtttcagacttttatttctcgccataattcaaactttttttctgataagctggttctcacttctgttactccagcttcttcggcacctgtttta")
 
#clean sequence automatically 
mySeq <- gsub("[^ACGTNacgtn]", "", mySeq)
#fas <- system.file("extdata", "IDH2.fas", package="DECIPHER")
dna_obj <- DNAString(mySeq)
#names(dna) <- "mySeq"
dna_rc <-reverseComplement(dna_obj)
####dna <- DNAStringSet(mySeq)
dna <-DNAStringSet(list(dna_obj, dna_rc))
names(dna) <- c("plus", "minus")

# plot the melt curve for the two alleles
temps <- seq(69, 75, 0.2)
####m <- MeltDNA(dna,
             #type="melt", temps=temps, ions=0.0165)
####matplot(temps, m,
        #type="l", xlab="Temperature (\u00B0C)", ylab="Average Theta")

m <- MeltDNA(dna, type="derivative", temps=temps -5, ions=0.0165)
matplot(temps, m, type="l", xlab="Temperature (°C)", ylab="-d(Theta)/dT")



legend("topright", names(dna), lty=seq_along(dna), col=seq_along(dna))



# plot the negative derivative curve for a subsequence of the two alleles
temps <- seq(80, 95, 0.25)
m <- MeltDNA(subseq(dna, 492, 542),
             type="derivative", temps=temps)
matplot(temps, m,
        type="l", xlab="Temperature (\u00B0C)", ylab="-d(Theta)/dTemp")
legend("topright", names(dna), lty=seq_along(dna), col=seq_along(dna))



# plot the positional helicity profile for the IDH2 allele
temps <- seq(90.1, 90.5, 0.1)
m <- MeltDNA(dna[1],
             type="position", temps=temps, ions=0.5)
matplot(seq_len(dim(m[[1]])[2]), t(m[[1]]),
        type="l", xlab="Nucleotide Position", ylab="Theta")
temps <- formatC(temps, digits=1, format="f")
legend("topright", legend=paste(temps, "\u00B0C", sep=""),
       col=seq_along(temps), lty=seq_along(temps), bg="white")

