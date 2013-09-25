RAD for genetic analysis of natural populations
===
Type II B RAD-seq processing pipeline: raw reads through population analysis

Type II B RAD uses special restriction enzymes to excise 36b DNA "tags" distributed evenly across the genome of (theoretically) any organism.  

Here, we used the enzyme BcgI, which recognizes the site ...∇(N10)CGA(N6)TGC(N12)∇..., to digest the genome of the scleractinian coral *Acropora millepora*.  


## Trimming Raw Reads
---
There are several schools of thought about trimming strategies.  Of course, you must trim to remove low-quality bases, but should you require each read to contain the restriction site or not?  Reads that do not contain the restriction site may be produced by either sequencing errors or enzyme fidelity errors.  Since we don't know much about this strange Type II B enzyme, it's within the realm of possibility that the enzyme may cut at slightly non-canonical sites.  But, do we then even care about these sites?  I can see reasons to either require or not require the site to be present in each read, but at the end of the day we decided to require it.  

 Jacob Malcom does NOT insist on each read having the restriction site: 
> "My guess is that your trimming is fine.  I don't think you want to require the recognition site 
	be present...I've found lots of good support that the enzymes are pretty darn high-fidelity, but 
	they're not perfect and well-supported variants occur in the 6-base recognition sequence.  I 
	filter out reads with > 20% of bases with QUAL < 20, but have found that Illumina data has very 
	few of such reads...so it might not be worth the time.  I remember reading a paper 6-12 months ago 
	suggesting that B-tails in the fastq qual string (ord(B) = 66, so if offset = 64 then QUAL = 2) 
	can be an issue.  And homopolymers can occur, but they don't like to map well anyway, so generally 
	not an issue."
	
Anyway, I used a script written by M. Matz (April 2013) to remove Illumina adaptors and require the restriction site for each read.  The output was then funneled into the fastxtoolkit to filter for low-qual reads:

```
module load fastx_toolkit
trim2bRAD.pl K208.fastq '.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}' 'AGATCGGAAGA' 0 | fastq_quality_filter -Q33 -q 20 -p 90 > K208.trim
```

## RAD Tag length
---
This enzyme leaves sticky ends when it cuts.  The question here is if we should analyze the full 36b length of the tag or if we should remove the two overhanging bases from each end of the read (32b reads).  I can imagine a world in which these ends have a higher error rate through ligation/degradation/etc than does the rest of the read, but I have not found any evidence of this.  N.B.:  Earlier, we did find evidence of many more SNPs located at positions 2 and 35 than in the rest of the reads, which indicated that we should just cut them off and move on.  However, I am convinced this is an artifact of using one of our homemade scripts to call genotypes and not a true signal.

Evidence: when I mapped 36b reads to the whole genome using Bowtie2 and then genotyped SNPs using the Genome Analysis Toolkit (GATK)'s UnifiedGenotyper, I collected **X** SNPs, which should be  X SNPs / 36 b = **Y** SNPs per base, IF they are evenly distributed.  
When I then used GATK's hard-clipping tool on the sorted.bam file to delete the first and last two bases, I collected **X** SNPs.  X SNPs / 32 b = **Y** SNPs per base.  These are the same numbers.  There is no SNP excess on the tails of the reads when you use robust genotypers and mappers.

## Reference to map reads to
---


### Population Statistics in general:
---

To filter by Minor Allele Frequency or not?  See: http://csb.scichina.com:8080/kxtbe/EN/abstract/abstract413250.shtml, "Effects of cutoff thresholds for minor allele frequencies on HapMap resolution: A real dataset-based evaluation of the Chinese Han and Tibetan populations"  Basically you should include all the SNPs.

Also see: https://www.biomedcentral.com/content/pdf/1753-6561-3-S7-S41.pdf ;"The effect of minor allele frequency on the likelihood of obtaining false positives"; namely, "These results suggest that removal of low MAF SNPs from analysis due to concerns about inflated false-positive results may not be appropriate"

Also, from an email exchange with Mark Kirkpatrick: 
> I'm working with linkage diseq and Fst at the moment.  The BayeScan manual (Fst) just notes that including low-MAF loci is bad and will skew the results, but doesn't elaborate on why.  I've asked the developer but he's pretty bad about email responses.  And for LD, r^2 is influenced by allele frequencies, so would that matter when comparing pairs of sites with super rare alleles?
"Yes.  Throwing out real data will bias your statistics."


## Linkage Disequilibrium
---
See Association Mapping of Kernel Size and Milling Quality in Wheat (Triticum aestivum L.) Cultivars, Flavio Breseghello*,† and Mark E. Sorrells 2006 : http://www.genetics.org/content/172/2/1165.full.pdf

Basically, you can either draw an arbitrary cutoff value at r^2 = 0.1, or you can mark the 95th percentile of the distribution of unlinked loci and intersect that line with the loess regression.

Also see : http://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/


### Population Connectivity:
---
## STRUCTURE:
Try to use as many SNPs as possible.  To filter for MAF or not?  Maybe not such a good idea...  But I can't find any decent documentation or "best practices".  Sarah says to not throw out any data; private alleles with a MAF of 0 are still good data!  So, maybe leave them all.  Run "5pops all SNPs" x 2, 3, 4k and "6pops all SNPs" x 2, 3, 4k.  

But, do filter out the outlier loci as identified by BayeScan.  

1.  Convert vcf into genepop  (vcf2genepop.pl)
2.  Delete loci that are implicated as outliers by BayeScan (remember, we did BayeScan several times: once per pairwise comparison, once for the five pops and once for all six pops, so consider/remove all outlying loci) 
	'''
	remove_loci_from_genepop.pl
	'''
3.  Split the genepop into random pieces (note: *SNPs are randomly selected, no attention paid to LD*): "genepop_split.pl six_pops.neutral.genepop (number of chunks) r"
3.  convert .genepop into .structure: java -jar /work/01408/rlc2489/RoxyTools/PGDSpider_2.0.4.0/PGDSpider2-cli.jar -inputfile five_pops.neutral.genepop -inputformat GENEPOP -outputfile five_pops.neutral.structure -outputformat STRUCTURE -spid genepop2structure.spid
4.  make mainparams, extraparams files
5.  run STRUCTURE
	**Time for running STRUCTURE, k=2
		30k SNPs:  1000 runs/2hrs = 1M runs/83days
		10k SNPs:  3000 runs/6mins = 30,000 runs/hr = 1M runs/33hrs
		 5k SNPs:
		 1k SNPs:
	       0.5k SNPs:

		Split it into smaller pieces: , try 10k SNPs/structure run




