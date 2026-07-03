# Comparative Genomics: Study and Analysis of Human Accelerated Regions (HARs) in Humans(hg38) and Chimps(PanTro6) 


### Project Overview

Human Accelerated Regions are regions in our genome that have been widely conserved in other primates but have mutated quickly in the human genome, and these regions are believed to be part of what makes Humans "special" from other primates. This project aims to analyze mutations that occured in Human Accelerated Regions primarily between Human sequences and Chimp sequences. The main question is to visualize the difference in DNA sequence between the two species as well as characterize the mutations in order to better understand if they contribute to human evolution or if they are due to non-adaptive biased GC conversions. 

### The Biological Problem

During DNA repair, there is a natural bias leaning towards G/C alleles being replaced rather than A/T, this is known as GC Biased Gene Conversion (gBGC). In evolutionary genomics, being able to identify these mutations and separate them from adaptive mutations that arose from selection will help us better pinpoint exactly where and how the human genome differs from that of other primates. 

### Metric 1: GC Content


Since DNA naturally undergoes deamination (C to T mutation,) regions in DNA that go against this trend and exhibit GC richness are an indication of gBGC. So measuring the GC content of a sequence helps us identify potential regions of gBGC.  

### Metric 2: Bias Ratio


Calculating a "bias ratio" helps determine which mutations are due to gBGC and which ones are potential adaptive mutations involved in human evolution, effectively helping us answer "What makes Humans special?" by eliminating false positives. The ratio divides the number of Weak to Strong (A/T to G/C) mutations by Strong to Weak (G/C to A/T) mutations.


### Technical Implementation


I organized this project using a modular design to ensure the code is reusable, scalable, and easy to document across multiple themes.



1) src/utils.py: I developed a comprehensive library of bioinformatics functions to handle the full data lifecycle:
- Data Acquisition: Implemented fetch_dna using the UCSC DAS API and Regular Expressions to automate sequence retrieval directly from genomic coordinates.
- Sequence Parsing: Includes custom functions like read_fasta and count_bases to handle standard genomics file formats and nucleotide distribution analysis.
- Evolutionary Metrics: Developed core functions for sequence comparison, including get_percent_identity, find_conserved_positions, and count_deaminations (identifying chemical decay).
- Mutation Analysis: The flagship analyze_mutations tool classifies substitutions into Transitions vs. Transversions and implements a W to S Bias Ratio to detect non-adaptive evolutionary signatures.



2) notebooks/: I utilized Jupyter Notebooks as a medium to document the step-by-step logic of the analysis. This includes coordinate normalization (handling the chimp 2A/2B fusion) and final data visualization of mutational patterns.



3) The Pipeline: I implemented an automated iteration engine that batch-processes orthologous sequence JSON files. The pipeline uses dictionary-merging logic (Splat operator) to join original genomic metadata with calculated percentages (GC content and Bias Ratios) for 281 regions.




### Key Results & Inference

Through batch analysis of the mutation patterns in the human lineage, I identified a specific region that exemplifies the "False Positive" problem in accelerate genomics.



- Candidate ID: ZOOHAR.1

- This region exhibited 9 total base substitutions with a Bias Ratio of 10.00.

- Conclusion: Every single mutation in ZOOHAR.1 is a Weak to Strong (A/T to G/C) transition. Because the directionality is 100% biased toward G/C alleles, I can infer that this region's "acceleration" is a product of non-adaptive gBGC. This confirms that high mutation rates in HARs can be driven by recombination glitches rather than adaptive natural selection.


### How to Run

Ensure all dependencies are installed using pip install requests.

Open notebooks/har_analysis.ipynb and run all cells to reproduce the mutation results from the raw sequence data.
