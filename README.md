# Bioinformatics Practice Package
Practicing writing basic bioinformatics algorithms to create a bioinformatics package.

### Classes

**FASTA_seq**
> Initialize with the sequence name and sequence. These can be pulled from a .fasta file with names being in the first line following a ">" and the corresponding sequence in the line bellow.

**FASTQ_seq**
> Initialize with the sequence "name", sequence, signs, and quality scores. These are the 4 lines for each seq in a .fastq file.

### Functions

**FASTA_open(FASTA_filepath)**
> Takes the argument of a filepath to a FASTA formatted .txt file. This fasta file can also be zipped (a .gz file) and this function will still be able to open the file. There is no need to unzip the file before using this function. Utilizes the FASTA_readlines(file) in order to read the lines of the file and create the FASTA_seq's. This will create new FASTA_seqs and return them as a list if there are multiple sequences within the file and will return the single FASTA_seq if there is only one sequence so it will not need to be unpacked later.
>>**FASTA_readlines(file)
>> This is a helper function to FASTA_open() and takes the file which must be inputted as a .txt and reads the lines of the file. It creates new FASTA_seq's and returns them as a list.


