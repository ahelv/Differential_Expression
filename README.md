# Differential_Expression

The code provided here can be used to analyse differential gene expression in R using Fisher's Exact Test. 

Required Inputs:
  1. bedfile of called peaks (i.e. from MACS2) with only chr, start, stop columns
  2. RNA-Seq data with 4 columns: ext_gene, logFC, pvalue, padj
  3. List of genes with 4 columns: gene name, chromosome, txStart, and txStop (each gene appearing only once)
  4. Flanking gene list (same genes but including flanking regions of 100kb)

Using the function: 
1. RNA is the RNA-seq file you would like to use here 
2. export_name1 is the name of the bedfile with peaks that aren't dysregulated (i.e. "genes_no_diff.bed")
3. peak_file is a file with the MACS2 called peaks of interest
4. export_name2 is the name of the bedfile of peaks for export  
5. peak_sort_name is the name of the sorted bedfile of peaks for bedtools 
6. nearest_gene_name is the name of the file with the nearest gene info from bedtools Closest
7. nearest_no_diff_name is the name of the file with the nearest gene info from bedtools Closest for the non-dysregulated genes
8. distance is the distance within a nearest gene will be searched (default is 100kb)

