## Methyl-SNP-seq
Scripts used for Methyl-SNP-seq related analysis. <br>
Developed and maintained by Bo Yan (New England Biolabs, yan@neb.com). <br>

---------------
### **Read_Processing** <br>
This directory contains scripts used for Methyl-SNP-seq data processing, including: <br>
adapter trim using TrimRead.py <br>
Reference-free Read Read Deconvolution using DeconvolutionConversion_v2.py <br>
Reference-dependent Read Deconvolution using DeconvolutionWithCalibration directory <br>
Correction of undeconvoluted Read using DeconvolutionUnmatchCorrect.py <br>
Removal of multiple mapping reads using MarkUniread.py <br>
Removal of duplicates using MarkDup.py <br>
Addition of XM tag for methylation extraction using AddXMtag.py <br>

Please see **README.pdf** for details. <br>

---------------

### **Methylation_Motif_Calling** <br>
This directory contains scripts used for methylation motif identificaion, including: <br>

### **motifExtraction.py** <br>
Used to extract a certain size of sequence (e.g. 8bp) containing a methylated C or C at a given position based on Methyl-SNP-seq reads. <br>
The extracted sequences can be used for methylation motif calling using identifyMotif.py and clusterMotif.py. <br>

Requirement: Python 3, regex module

Logic: <br>
If --contigs is provided: <br>
&emsp;&emsp;Use the reads mapped to the target contig(s) for sequenc extraction. <br>
Based on the methylation report, extract the nearby sequence of methylated C (methylC.sequence) or unmethylated C (unmethylC.sequence) or all C (allC.sequence). <br>

Usage: <br>
<pre>
<b>python motifExtraction.py --input Ecoli_ABK_Deconvolution_R1.fq --report Ecoli_ABK.Deconvolution.5mC -left 1 -right 4 --sam Ecoli_ABK.sam --name Ecoli_ABK_Node-27-70 --contigs NODE_27_length_74030_cov_17.136755 NODE_70_length_61531_cov_18.801319</b>
</pre>

Parameters: <br>
--input: Methyl-SNP-seq Deconvolution Read fastq file <br>

--report: corresponding Deconvolution methylation report <br>

--name: name for all output files, Default TestMotif <br>
Generate Output files: <br>
name.Deconvolution_R1.fq: containing reads used for seq extraction if --contigs is on. <br>
name.Deconvolution.5mC: containing methylation information for reads in name.Deconvolution_R1.fq. <br>
name_motif_methylC.sequence, name_motif_unmethylC.sequence and name_motif_allC.sequence, e.g. <br>
e.g. <br>
CCAGGC <br>
CCAGGA <br>
.. <br>
Note: <br>
name_motif_methylC.sequence and name_motif_allC.sequence can be used as input files for identifyMotif.py. <br>

--left/-l INT default 1, --right/-r INT default 4: <br>
-l and -r define the number of bases upstream of downstream of a methylated or unmethylated C. <br>
With default setting, seq is x5mCxxxx in methylC.sequence; xCxxxx in unmethylC.sequence; xC/5mCxxxx in allC.sequence. <br>

--contigs: chr(s), nargs="+",  Not required <br>
If provided, only the reads mapped to the given chr(s)/contig(s) are used to extract sequences. <br>
These given chr(s) must be present in the sam file otherwise no read will be used. <br>

--sam: Not required <br>
Deconvolution Reads mapped to genome assembly, only required if --contigs is provided. <br>

### **identifyMotif.py** <br>
Use to identify the significantly methylated motif by comparing the counts in sample and in reference based on Binomial test and Bonferroni correction. <br>
null hypothesis: the methylation level of motif in sample is not significantly higher than Pi_0 (defined by --mode). <br>

Requirement: Python 3, scipy>=1.6.0 <br>

Usage: <br>
<pre>
<b>python identifyMotif.py --sample Ecoli_motif_methylC.sequence --reference Ecoli_motif_allC.sequence --output Ecoli_SignificantMotif.txt</b>
</pre>

Parameters: <br>
--sample, --reference: <br>
motif sequences served as positive (seq containing methylated C) and background (seq containing all C) for binomial test. <br>
e.g. <br>
CCAGGC <br>
CCAGGA <br>
GCAACA <br>
... <br>

Test Files provided: --sample Ecoli_motif_methylC.sequence --reference Ecoli_motif_allC.sequence <br>
In the sample file, the C at position 4 is methylated C; in the reference file, the C at position 4 is methylated or unmethylated C. <br>

--output: the first line is header <br>
Name&emsp;Sample&emsp;Reference&emsp;pvalue&emsp;corrected_alpha&emsp;Significance <br>
GCCCAGGT&emsp;157&emsp;285&emsp;0.0&emsp;6.105378838756945e-09&emsp;True <br>
CACCAGGC&emsp;282&emsp;429&emsp;0.0&emsp;6.105378838756945e-09&emsp;True <br>

--alpha: Default 0.0001 <br>
alpha value used to determine the corrected pvalue for significance. <br>

--cutoff: Default 2 <br>
To avoid calling significance of motif seq showing accidentally in sample, only the motif seq having count>=cutoff in the reference will be called as significant if adjusted pvalue below sigificance value. <br>

--mode: options [average, top, a float number] <br>
controling the Pi_0 used for binomial test <br>
null hypothesis: the methylation level of seq in sample is not significantly higher than Pi_0 <br>
&emsp;&emsp;average: Default mode, Pi_0 = (number of methylated seq in sample) / (number of all seq in reference) <br>
&emsp;&emsp;top: Pi_0 is the 95 percentile of the methylation level of all the seq in sample <br>
&emsp;&emsp;a float value: giving a custom value betweeen 0 and 1 for Pi_0, e.g. 0.67 <br>

### **clusterMotif.py** <br>
Use to cluster the significant motif seqs defined by identifyMotif.py. <br>

Requirement: Python 3, scipy>=1.6.0, matplotlib, seaborn, pandas <br>

Logic: <br>
Calculate the distance (1bp mismatch equals to distance=1) between each pair of the motif seq; <br>
based on distance matrix, cluster the motif seqs; <br>
determine the number of clusters (p_number) to cut the dengram tree; <br>
report the motif seqs in each cluster and count the base composition at teach position, <br>
the returing motif seqs can be used for motif logo calling using weblogo. <br>

Usage: <br>
<pre>
<b>python clusterMotif.py --input Ecoli_SignificantMotif.txt --name Ecoli --number 2 > report.txt </b>
</pre>

Parameters: <br>
--input: File generated by identifyMotify.py containing the significant motif seq, e.g. <br>
Name&emsp;Sample&emsp;Reference&emsp;pvalue&emsp;corrected_alpha&emsp;Significance <br>
GCCCAGGT&emsp;157&emsp;285&emsp;0.0&emsp;6.105378838756945e-09&emsp;True <br>
CACCAGGC&emsp;282&emsp;429&emsp;0.0&emsp;6.105378838756945e-09&emsp;True <br>
.. <br>

Test Files provided: --input Ecoli_SignificantMotif.txt  <br>
Ecoli_SignificantMotif.txt is generated by identifyMotif.py using provided test files with default mode. <br>

--name: Name for the heatmap and dengram images and cluster seq files, Not Required <br>
name_Motif_cluster_heatmap.png, this heatmap helps the determination of p_number <br>
name_Motif_cluster_dengram.png, cut tree returning a dengram in which the number of clusters is determined by p_number <br>
name_Motifcluster1, name_Motifcluster2 ... These files have the motif seq that can be used to call motif logo. <br>

--number: p for hierarchy.dendrogram, determining the number of clusters generated by cut tree <br>

Note: <br>
Returning is redirected and saved in report.txt: <br>
The number of motif seqs in this cluster is: 64 <br>

|   | Pos1 | Pos2 | Pos3 | Pos4 | Pos5 | Pos6 | Pos7 | Pos8 |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| A  | 16  | 16  | 0  | 0  | 64  | 0  | 0  | 16  |
| T  | 16  | 16  | 0  | 0  | 0  | 0  | 0  | 16  |
| C  | 16  | 16  | 64  | 64 | 0  | 0  | 0  | 16  |
| G  | 16  | 16  | 0  | 0  | 0  | 64  | 64  | 16  |

The number represents the number of seqs having A/T/C/G at this position. <br>
So the above motif is NNCCAGGN. <br>
The returning motif seqs can be used to call motif logo using weblogo. <br>

