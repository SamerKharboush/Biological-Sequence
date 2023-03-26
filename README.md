# About Dataset

Drosophila Melanogaster
Drosophila Melanogaster, the common fruit fly, is a model organism which has been extensively used in entymological research. It is one of the most studied organisms in biological research, particularly in genetics and developmental biology.

When its not being used for scientific research, D. melanogaster is a common pest in homes, restaurants, and anywhere else that serves food. They are not to be confused with Tephritidae flys (also known as fruit flys).

# First jupytor notebook
Folder [Biopython]

Getting Started with Biopython
Myles O'Neill - based on http://biopython.org/DIST/docs/tutorial/Tutorial.html

Biopython is a library we can use to analyze bioinformatic data. Lets have a look at some of the things we can do with it!

This script will go through a selection of tutorial excercises from the Biopython cookbook. I've tweaked them a bit to work with the Drosophila dataset we are using


# Seconed Jupytor notebook
Folder [notebooks]

![image](https://user-images.githubusercontent.com/88907921/227752967-5695840e-fa19-4eeb-8d9b-969450ef986d.png)

Whilst there are wonderful libraries like BioPython that allow us to work with biological sequences
It's often quite benefitial to understand the workings of the code & expand on the code if a need arises, in this notebook, we use a very simple class structure
This is more time consuming, nevertheless it's definitely more interesting & allows us to incorporate different libraries (eg.bokeh)

In this notebook, emphasis will be placed on replicating something similar to two BioPython modules:
from Bio import SeqIO ( Class for readng sequences )
from Bio.Seq import Seq ( Class for Sequence Operations )

At the very end of this notebook, we'll also create a whl package, which is uploaded & updated in the dataset: bioseq, which includes classes used here
whl is the format of packages used on pip & these files can be uploaded to datasets, installed via pip & called as you would an existing pip package on Kaggle
If you are interested in bioinformatics & I'd definitely recommend getting to know bioconductor as well, I've made a simple introductory notebook bioconductor - If you are interested in bioinformatics & I'd definitely recommend getting to know bioconductor as well, I've made a simple introductory notebook bioconductor - bioinformatics basics as well
1 | Background
1.1 | Section Overview

In this section, we'll introduce some key terminology associated with what we will look into in this notebook:

1.2 - We'll introduce information about biological cells
1.3 - We'll introduce the two commonly used alphabets (neucleotides & amino acids)
1.4 - DNA is made up of two strands, we'll mention a couple of words about them
1.5 - In this section we'll introduce two main sythesis methods (RNA & Protein)
1.6 - We'll talk about what we'll be making in this notebook, in order to do some basic biological seqence operations
1.2 | Cell Information

Cells contain molecules which will reference to througout the notebook: amino acids, nucleotides, proteins

A cell is mostly composed of water:

a bacteria cell has a weight composition of roughly 70% water & 30% chemical origin,
7% -> small molecules Inc. amino acids and nucleotides
23% -> macro molecules (proteins,lipids,polysaccharides)
According to their internal structure, they can be divided into to major categories:

Prokaryotic cells : have no nucleus or internal membranes
Eukaryotic cells : which have a defined nucleus, internal membranes and functional elements called organelles.
At a structural level, all cells are surrounded by a structure called cell membrane or plasma membrane.
This membrane is permeable to molecules that cells need to absorb from or excrete to the outside medium.
Within the cell we find the cytoplasm (largely composed of water), which serves as the medium for the cell
1.3 | Sequence Alphabets

ABC (I/II) Nucleic Acids

Among molecules with a biological role, we can find nucleic acids
Nucleic acids encode and express the genetic code that is kept within the cell
There are two major types of nucleic acids:
Deoxyribo Nucleic Acid (DNA)
Ribonucleic Acid (RNA) (Obtainable via transcription)
DNA contains the information necessary to build a cell, and keep it functioning.
In eukaryotic cells, DNA will be found in the nucleus, whilst in the prokaryotic cells, it will be found in the cytoplasm.
IUPAC defines the full list of nucleotides as shown in the table below, with A,T,G,C being the main four:
Another type of nucleotide list often used is IUB Ambiguity Codes, which we use later in the notebook as well
ABC (II/II) Amino Acids

Amino acids:

The building blocks of proteins, which are macromolecules that perform most of the functions inside a cell

Proteins have a broad range of functions, spanning from catalytic to structural functions:

Enzymes : Type of abundant proteins that promote chemical reactions & convert some molecules into other types of molecules required for the functioning of the cell
Carbohydrates : Serve as energy storage, both for immediate and long term energy demands
Lipids : Part of the plasma membrane, doing signaling and energy storage
The cell also contains other components of varying complexity. Of importance:

Mitochondria & the Chloroplasts : Organelles involved in the production of energy.
Ribosomes : Large and complex molecules composed of a mixture of genetic material, req. to assemble proteins and play a central role in the flow of genetic information
1.4 | DNA Strands

COMPLEMENTARY STRANDS IN DNA
DNA is a molecule composed of two complementary strands that form and stick together due to the connections established between the nucleotides in both strands
This is made possible by due to the chemical phenomenon:
Where Adenine (A) bonds only with Thymine (T) nucleotides; as a result of two hydrogen connections
Similarly, Guanine (G) bonds only with Cytosine (C) nucleotides by three hydrogen connections
REVERSE COMPLEMENT
This results in two complementary and anti-parallel strands (connected in opposite directions), if we know the nucleotide sequence in one of the strands, we can get the sequence in the opposite strand by taking the complement of its nucleotides, which are also read backwards, thus we have the reverse complement of the other strand.
It has become a standard to describe the DNA though only one of the strands, due to this complementarity using [A,T,G,C].
The existence of these two strands is essential in order to pass on genetic information to new cells and produce proteins.
1.5 | Central Dogma of Molecular & Cell Biology

The process talked about in the next section is quite complex & not all aspects are competely understood, a general and simplistic picture is provided below
DNA, RNA & Proteins are the central elements of the flow of genetic information that occurs in two steps:
STEP 1 RNA SYNTHESIS : TRANSCRIPTION

Preliminary step required to produce a protein
The nucleotide sequence of a gene from one of the DNA strands is transcribed ( copied into a complementary molecule of RNA )
The complementarity of the genetic code allows recovering the information encoded in the original DNA sequence, a process performed by the enzyme, RNA polymerase
Additional steps of RNA processing, including stabilising elements at the end of the molecule, are performed by different protein complexes
After these steps, which occur within the nucleus of the cell, an RNA molecule; mature messenger RNA (mRNA) is obtained.
The mRNA is then transported to the cytoplasm, where it will be used by the cellular machine to guide the production of a protein
STEP 2 PROTEIN SYNTHESIS : TRANSLATION

Process in which the nucleotide sequence of the mRNA is transcribed into a chain of amino acids, forming a polypeptide.
Proteins are cellular entities that have either:
(a) Structural function - Participating in the physical definition of a cell.
(b) Chemical function - Being involved in chemical reactions occuring in the cell.
In order to function as expected, a protein needs to acquire the appropriate structure & this structure is often decomposed at different complexity levels
The primary structure is defined by the chain of amino acids; polypeptide, consisting in part or completely of protein.

This process is performed by the ribosomes that attach and scan the mRNA from one end to the other, in groups of nucleotide triplets / Codons.
CODONS

In each position of the triplet, we have 1/4 nucleotides, ie. there are 4x4x4 (64) possible triplets/Codons.
For each codon in the mRNA sequence, we have a corresponding amino acid in the polypeptide chain
START AND STOP CODON

Some of these codons represent significant signals that indicate the initiation or the termination of the translation process.
Once the ribosome detects an initiation codon, it starts the formation of the amino acid chain,
When it scans the stop codon, it stops the translation and detaches from the mRNA molecule
There are 20 types of amino acids used to form polypeptides (for IUPAC) & less than the 64 possible codons, therefore we have more than one codon corresponds to a type of amino acid
During the translation process:
A type of small RNA molecule, transfer RNA (tRNAs) will bring to the ribosome, the amino acids of the corresponding type, which will be complementary to the mRNA codon that is currently being scanned.
Each mRNA molecule can be scanned multiple times by different ribosomes, giving rise to multiple copies of the polypeptide.
With its redundency, where more than one codon encodes an amino acid, the genetic code encloses a very efficient code-correction mechanism that minimises the impact of errors in the nucleotide sequence occuring in DNA replication.
OPEN READING FRAMES

During translation process:

Parsing of the mRNA sequence by the ribosome may start at different nucleotides.
Given that a codon is composed of three nucleotides, the mRNA sequence may have 3 possible interpretations.
These three ways of parsing the sequence are called reading frames.
2 | Sequence Operations Class
We define a custom class, SQ, which contains basic sequence related operations
The sequence class SQ can be coupled with other classes that incorporate more sophisticated sequence based operations; eg. sequence alignment
CLASS | Sequence Class

ATTRIBUTES

seq : String format of sequence
seq_type : Sequence Character Type
id : Sequence Identifier
name : used to allocate a name to the sequence
description : used to show more information about the sequence
annotation : used to store annotation information about our sequence
CLASS METHODS

info : Show the sequence & type information only
abc : Show the alphabet characters of the the sequence type
validate : Determine if the current sequence contains valid subsets of charaters for its type
freq : Calculate & plot the frequencey of each character in the alphabet of the define sequence type
count_purines : Count purines A&G sum & pyrimidines C&T
groupfreq : Cound dinucleotides & trincucleotide groupings in sequence
gc : Calculate the GC content of the sequence
reverse_comp : Rearrange the DNA sequnce to show its reverse complement sequence pair
transcription : Convert DNA sequence into an RNA sequence
get_protein : Derive all the amino acid chains in all reading frames and store as one sequence
find_pattern : Find standard/string patterns in the sequence & more specialised/specific patterns
cut_pattern : Cut sequence in specific parts based on a pattern
MAPPING DICTIONARY

3 | Basic Sequence Operations
3.1 | Section Overview

In this section, we'll look at basic operations associated with our classes, focusing on just using the availble commands assoiated with blocks
3.2 - Define a sequence from a python string
3.3 - A string must contain characters from either the nucleotide or animo acid alphabet, let's check if our sequence has any missprints
3.4 - All DNA strings come in pairs, we always input one of its strands, let's obtain the reverse of the complementary strand
3.5 - sequences in SQ only need to be initialised with the sequence itself, however we can also add mode information about out sequence
3.6 - Each sequence may contain subsequences or regions of interest, we can add annotations that will specify the region location & a description about this subsequence

3.2 | Defining Basic Sequences from Strings

We can create an instance of object SQ using the following ( including the sequence type )
The sequence and sequence type summary can be view by calling info method
DNA sequences are also colour coded, making it a little easier to distinguish the different characters
We can also view the sequence by calling the view method, which accepts an argument split_id (indicating after how many characters a pause should be made when visualising it)

3.3 | Check Validity of Strings

A sequence is valid if its alphabet corresponds to a specific set of code, class SQ uses the IUPAC code standard
We can check whether a sequence is valid or not by using the validate function, returning a logical output
This can be handy when our sequece is very large, and we can check the entire sequence for errors

3.4 | Switch to its Reverse Complement

DNA has two complementary strands
Due to the complementarity of the DNA strands, usually only one of the strands is provided in a sequence files obtained from databases
The second strand to the input DNA sequence can be obtained by calling the reverse_comp function

3.5 | Defining a Detailed Sequence

SQ can also be used to store more detailed information about our sequence; id and description
Later we'll need to look at a sequence & jot down some information about parts of the sequence after we have identified them as well

3.6 | Making Annotations

Say we searched a database to find if our sequence contains something interesting
We have two hits, 'TTTT','CCTCA', which we want to save in annotation
Usually, when we search databases, we'll have a sequence identifier & a description/comment about the sequence, so we'll be able to save this info
The search result sequences can be extracted from the alignment format & saved in SQ format

4 | Counting Alphabet Characters
4.1 | Nucleotides % in a Sequence

We may want to obtain the count for each each nucleotide or amino acid in the given sequence that has the same amount of bases
We can use the freq method to show:
The count/percentage of each of the alphabet compare two sets of sequences with the addition of compare

4.2 | Counting Amino Acids in a Sequence

We may want to compare sequences of different length, so in that case it's best to use percentage (default)
show_id='count' is good for when we just want to see the character count of sequence

4.3 | Di/Trinucleotides

We may want to know the different combinations of double or tripple nucleotides that is present within the sequence
One example would be to count all the codons present in a sequence

4.4 | GC Content

GC-Content | WIKIPEDIA & SCIENCEDIRECT

In molecular biology and genetics, GC-content (or guanine-cytosine content) is the percentage of nitrogenous bases in a DNA or RNA molecule that are either guanine (G) or cytosine (C). This measure indicates the proportion of G and C bases out of an implied four total bases, also including adenine and thymine in DNA and adenine and uracil in RNA.
GC content is strongly correlated with biological features of genome organization such as distribution of various classes of repeated elements, gene density, level and tissue-specificity of transcription, and mutation rate.

4.5 | Purines & Pyrimidines

Purines - (adenine and guanine) have a two-ringed structure consisting of a nine-membered molecule with four nitrogen atoms
Pyrimidines - (cytosine, uracil, and thymine) only have one single ring, which has just six members and two nitrogen atoms

5 | Decoding DNA for Protein Generation
5.1 | Transcription

The mRNA is created as the two complementary strands are split, and is complement to one of the strands.
Transcription of DNA sequence can be obtained using the .transcription() function, where the T is replaced by U in the sequence string.
rna_seq = sq_n.decode.transcription()
rna_seq.info()
SEQ: AUGACGGAUCAGCCGCAAGCGGAAUUGGCGUUUACGUACGAUGCGCCGUAA TYPE: rna
5.2 | Translation

The background to what happens during the translation process was already outlined in 1.5, however it makes sense to expand on theory a little more.
Proteins are synthesised by creating chains of aminoacids, according to information contained in the messenger RNA (mRNA) in a process called translation.
START & STOP CODONS
Translation of a protein:
Always begins with a specific codon ATG -> M (Methionine), which is always the first amino in the protein.
Translation process terminates when a stop codon is found; TAA,TAG,TGA -> _.
An example sequence where we know exactly where the start and termination codons are:
ATGACGGATCAGCCGCAAGCGGAATTGGCGTTTACGTACGATGCGCCGTAA
We can note it starts with ATG & ends with one of the three stop codons.
As the sequence follows the rule of start and ending codon, in such a case we can use the staticmethod defined in SQ; .translate(p0=0) directly.
Some interesting discussions about the starting amino acid (ATG is usually the starting codon)
OPEN READING FRAME
![image](https://user-images.githubusercontent.com/88907921/227753445-d5b3401c-f12b-4d5b-a0b0-8f43d6dfa884.png)

A reading frame is a way of dividing the DNA sequence into a set of consecutive, non-overlapping triplet nucleotides (possible codons) (using dictionary mapping).
A given sequence has 3 possible reading frames, first, second and third nucleotide positions. In addition, considering there is another complementary strand, we should compute the only 3 frames corresponding to the reverse compliment.
In many cases, given a DNA sequence, we don't know in advance where the coding regions are, especially when dealing with complete sequences.
In such cases, we need to scan the DNA sequence for the coding region. First, we need to divide and compute these reading frames (6 in total). The frames function stores the converted 6 converted amino acid strings in a list.
OPEN READING FRAME | genomove.gov

An open reading frame is a portion of a DNA molecule that, when translated into amino acids, contains no stop codons
The genetic code reads DNA sequences in groups of three base pairs, which means that a double-stranded DNA molecule can read in any of six possible reading frames: three in the forward direction and three in the reverse.
A long open reading frame is likely part of a gene.
Open Reading Frames (ORF)[1]
Drawing
IDEAL CASE (P=0)
Let's try one of the reading frames at the start of the sequence (p0=0), which we know follows the correct rules observed in life.
Let's also consider all other possible ORFS & get the final protein found in the DNA sequence.

Some interesting information about how many amino acids do we need to actually classify the decoded result as a protein

proteins & amino acids (Both animal and plant proteins are made up of about 20 common amino acids)
We will mention those with less than about 20 amino acids below when we'll look at an example
6 | Loading Files Containing Sequences
REAL SEQUENCES
Most realistic application of bioinformatics certainly involve working with sequences that are too long
Often we need to work with databases, so it's convenient to have an effective I/O system in place
One of such formats is FASTA, its made to be a very flexible format, allowing us to work with single, multiple or even alignment sequences.
THE FASTA FORMAT
FASTA | FASTA format:

In bioinformatics and biochemistry, the FASTA format is a text-based format for representing either nucleotide sequences or amino acid (protein) sequences, in which nucleotides or amino acids are represented using single-letter codes. The format also allows for sequence names and comments to precede the sequences. The format originates from the FASTA software package, but has now become a near universal standard in the field of bioinformatics.

Commonly used to store nucleotide or protein sequences. It's less detailed & usually containing only the sequence, and name/header only.
File extension | Typically changed based on the sequence type content of the file, defined below in the dictionary FASTA_dic.
Multisequences | The format can contain any number of sequences, each starting with the symbol: >.
The name/header, defined after the symbol > contains an origin identifier (NCBI identifiers), defined in the the dictionary identifiers_dic.
The NCBI defined a standard for the unique identifier used for the sequence (SeqID) in the header line.
This allows a sequence that was obtained from a database to be labelled with a reference to its database record.
The database identifier format is understood by the NCBI tools like makeblastdb and table2asn.
The following list describes the NCBI FASTA defined format for sequence identifiers
CLASS | read_seq

READING SEQUENCE FILES CLASS
We define a custom class, read_seq(), which methods used for reading the FASTA format
The class automatically should detect the class type and call the corresponding class, storing the sequence upon instantiation
Sequences that are read are stored in SQ objects, which allows us to also to also store additional information about the sequences
CONSTRUCTOR

name : stores the pathway to the file
format : stores the FASTA format identifier
CLASS METHODS

read_FASTA : used to read FASTA formats, storing sequence data in lst_seq & descriptions in SQ objects
get_sq : used to store the return a list or a sequence of SQ class objects

Some interesting information about how many amino acids do we need to actually classify the decoded result as a protein

proteins & amino acids (Both animal and plant proteins are made up of about 20 common amino acids)
We will mention those with less than about 20 amino acids below when we'll look at an example
6 | Loading Files Containing Sequences
REAL SEQUENCES
Most realistic application of bioinformatics certainly involve working with sequences that are too long
Often we need to work with databases, so it's convenient to have an effective I/O system in place
One of such formats is FASTA, its made to be a very flexible format, allowing us to work with single, multiple or even alignment sequences.
THE FASTA FORMAT
FASTA | FASTA format:

In bioinformatics and biochemistry, the FASTA format is a text-based format for representing either nucleotide sequences or amino acid (protein) sequences, in which nucleotides or amino acids are represented using single-letter codes. The format also allows for sequence names and comments to precede the sequences. The format originates from the FASTA software package, but has now become a near universal standard in the field of bioinformatics.

Commonly used to store nucleotide or protein sequences. It's less detailed & usually containing only the sequence, and name/header only.
File extension | Typically changed based on the sequence type content of the file, defined below in the dictionary FASTA_dic.
Multisequences | The format can contain any number of sequences, each starting with the symbol: >.
The name/header, defined after the symbol > contains an origin identifier (NCBI identifiers), defined in the the dictionary identifiers_dic.
The NCBI defined a standard for the unique identifier used for the sequence (SeqID) in the header line.
This allows a sequence that was obtained from a database to be labelled with a reference to its database record.
The database identifier format is understood by the NCBI tools like makeblastdb and table2asn.
The following list describes the NCBI FASTA defined format for sequence identifiers
CLASS | read_seq

READING SEQUENCE FILES CLASS
We define a custom class, read_seq(), which methods used for reading the FASTA format
The class automatically should detect the class type and call the corresponding class, storing the sequence upon instantiation
Sequences that are read are stored in SQ objects, which allows us to also to also store additional information about the sequences
CONSTRUCTOR

name : stores the pathway to the file
format : stores the FASTA format identifier
CLASS METHODS

read_FASTA : used to read FASTA formats, storing sequence data in lst_seq & descriptions in SQ objects
get_sq : used to store the return a list or a sequence of SQ class objects

7 | Recognising Patterns
Another important part relating to sequences are pattern recognition within a sequence
There are specific patterns of nucleotides we can try to find, in order to understand our sequence content a little more
Having looked at the decoding of proteins in Section 4, we might want to know more about our amino acid chains
For example, an amino acid chain containing the pattern Zinc finger RING-type | Prosite Database

Some proteins known to include a RING finger include Mammalian breast cancer type 1 susceptibility protein (BRCA1)
Let's load the Breast cancer type 1 susceptibility protein BRCA1_HUMAN in FASTA format (from the UniProtKB database)
Using the loader from Section 5 to see if we have a match as well, since they should be related
Other ZnF types of prosite information @github
PROSITE PATTERN FORMAT
Looking at the example, it should be quite clear what the below format rules mean,
Last point not currently added so you'd need to be specific, which case you're actualy searching for eg. x(2),x(3),x(4) instead of x(2,4)
Standard IUPAC amino acid used to as bases in pattern, separated by -
x -> any amino acid acceptable
[] -> ambiguity represented by list, any aa in that list acceptable
{} -> ambiguity represented by list, any aa other than in {} accepted
Repetition of pattern element shown below: x(3) -> x-x-x, x(2,4) -> to x-x or x-x-x or x-x-x-x

8 | Creating Subsequences
8.1 | Creating Subsequences

Our entire sequence may be very long & we might have want cut the current sequence for a variety of purposes
Example applications:
Where a full sequence was cut into non overlapping unitigs - Identifying Antibiotic Resistant Bacteria
Where a full sequence was cut into a predefined sets of nucleotide bases & OHE was used - Transcription Factor Binding Location Prediction
Restriction Enzymes can cut DNA sequences at particular segments of the DNA, defined by a particular pattern (we'll look at this here)
8.2 | Dividing Sequences into Parts

CREATING LISTS OF SUBSEQUENCES
Whatever the mechanism, the function cut_pattern is general, requiring us to specify where in the sequence a cut occurs using "|"
# Arbitrary cut pattern at |
n_seq = SQ('AAGATTCGAGCATGCAAACCGGATACA',seq_type='dna')
n_seq.cut.cut_pattern('AGC|AT')
['AAGATTCGAGC', 'ATGCAAACCGGATACA']
8.3 | Generating OHE DNA Features

We can cut the sequence at a starting location & number of characters that must be present in each sequence
characters which aren't sufficient to create another segment; the leftover characters are discarded

9 | Covid-19 : Protein Identification
9.1 | Section Overview

Having defined a class that can read, store (read_seq) & derive proteins (SQ) from a sequence file in the FASTA format:
We can read a sequence from a database and store it in our sequence class & derive proteins/amino acid chains from the DNA sequence
Once we have derived our amino acid chains, we will want to identify them, so we need another class that will be able to search a database and match our derived sequence to the ones found inside the database, which we will be able to do using BLASTwww
BLAST - DATABASE SEARCHING
Steps we'll be taking:

(1) Get the proteins encoded in the genome ( using .translate )
(2) Identify the particular protein present & pick some that interest us
(3) Using a local sequence alignmnent tool called BLAST; (using class BLASTwww)
(3) Visualise the alignment, so we can see which parts of the sequence match to that of the database
We will look into Local Sequence Alignment in another notebook Biological Sequence Alignment, which BLAST uses.
However for the time being keeping things simple; a protein sequence is compared to others already found in a database and returns very simular sequences
9.2 | Reading Sequence from Databses

USING THE NCBI API (via Biopython)
To use the class function BLASTwww, we'll first need to:

(1) decode the proteins from the DNA sequence
(2) select a specific protein & search for matching results in the database
BLASTwww - Uses SearchIO.read in Biopython, which enables us to visualise the alignment using view

9.3 | Protein Identification

READ SEQUENCE FROM DATABASE
Let's use the coronavirus genomoe/sequence, already uploaded to Kaggle & original source
9.4 | Manually Searching Other Databases

There are quite a few databses, we only checked the pdb database & not all are accessible via the NCBI API
We might not actually find a match for the entire sequence all the time, especially if it's as long as the one we have, naturally it will likely contain different subsequences
We can use this string; copy & paste this protein & search for it using the PSI BLAST in order to find proteins already in a predefined database.
For the protein sequence that we searched for:
The results indicate that a 100% match was found with an already existing sequence in the database: sp|P0C6U8|R1A_SARS Replicase polyprotein 1a

What is a polyprotein | NCBI

Polyproteins are chains of covalently conjoined smaller proteins that occur in nature as versatile means to organize the proteome of viruses including HIV.

10 | Summary
CUSTOM CLASSES
In this notebook, we have defined a class that stores basic information about a specific sequence of interest.
Another class was also defined in order to read a commonly encountered format type FASTA.
OVERVIEW
We also looked at common sequence related operations such as transcription & translation, which are quite critical concepts.
We also looked at some basic visualisation of the nucleotide & amino acid distributions, which can be helpful when comparing not only the count, but also the percentage of each alphabet for the corresponding sequence type.
WHERE TO FROM HERE
As with the BioPython module, we have defined classes which store sequence data & which can be expanded further.
The class, SQ will be used as a basis for other sequence related operations, such as pairwise & multiple sequence alignment, which due to their content size alone, can have their own dedicated classes, similar to the class READ_SEQ, which is dedicated to reading and passing on the sequence information to class SQ.
The notebook will of course be continuously updated and futher functions will be added in the future.
