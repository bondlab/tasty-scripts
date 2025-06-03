# tasty scripts

Scripts for annotating and poking around cytochromes

The HemeHunter 3000 will identify multiheme c-type cytochromes from a GenBank file or directory of files.

This script scans translated CDS features and selects those with (default) at least 3,
or the threshold set with -m, to find that many 'CXXCH' motifs. before looking for non-cannonical motifs.

It is careful to only use -CH dipeptides once, keeping track of motif locations, to avoid
common errors of double-counting motifs

It calculates total motif counts, molecular weight, and motif density, and outputs a summary TSV.
Use the -f flag to also save a FASTA file of the cytochrome amino acid sequences.

If you need translations, or want to look for missing genes, use the option -t 
which will automatically translate any untranslated regions, including truncated genes.
It will then use the new _translated.gbk version of the file to call all cytochromes

It will autodetect if you provide a directory or single file, and process appropriately

Results are ready to be layered onto another annotation file, such as with featurebot

-----


The Featurebot is a multipurpose tool to map annotations onto genomes. 

It will update a GenBank file with output from SignalP, TMHMM beta barrel, or multiheme cytochrome analyses. 
It can run in single mode (use -gff and -gbk for feature and genome file), or batch (specify -gbkdir and -feturedir). 

In batch mode, it will match feature files to genome files if acession number is filename prefix, 
like GCF_000816575.1_regionoutput.gff3. Use -proteinid to use protein_id instead of locus_tag to match CDS features to annotations

----

Other scripts are specific for moving files around in post-processing, likely are not useful