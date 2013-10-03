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

## Reference to map reads to: 
---
1.  Extracted _in silico_ digested fragments?  Allowing a single mismatch in the restriction site?  Forcing perfect matches?
2.  The whole genome?

## Mapping algorithm:
---
1.  BWA
2.  Bowtie1 (not indel-aware/doesn't support gapped alignments)
3.  Bowtie2 
4.  SHRiMP
5.  Stampy

## Mapping settings/switches:
---

## Genotype callers:
1.  Samtools mpileup
2.  in-house haplotype caller
3.  GATK UnifiedGenotyper
4.  GATK HaplotypeCaller (can't get it to work)

## VCF filtering
The "Ultimate vcf" file (26 Sept 2013): 
	
1.  filtered for retaining 90% of "true variants" (as determined by the GATK UnifiedGenotyper VQSR stage)
2.  has only loci that are genotyped in >80% of individuals
3.  has only individuals in which >60% of all loci are genotyped.  
4.  Remove Reef 21-121, those other-species guys from Night (Ni15, Ni17, Ni18, O166)  (`delete_inds_from_vcf.pl`)
5.  Remove clones (high relatedness; `find_clones.pl`,`delete_inds_from_vcf.pl`) The following pairs are at least 95% similar:
 	* K210 = K212 = K213 = K216 (keep K212, has the most genotypic information of the four)
	* A395 = A397 (not using A reef any more)
	* K219 = K211 (keep K211)
	* A412 = Ni16, (delete Ni16, likely a different species)
	* K4 = O5, (delete both; likely some labeling screw-up)
	* M17 = M16 (delete M16; M17 has more genotypic information)
6.  Remove all Reef 21-121 (super high Fst between those and all other reefs - 0.2, are likely a different species),  
7.  Remove those individuals with low heterozygosity (`vcftools --vcf ultimate.vcf --het`) 
  	* delete those with F values > 0.4 (only Ni15, Ni17, Ni18, which are already deleted)
	* delete those with F values > 0.3 (Ni2, Ni10, Ni3, Ni5, O3, Ni16)  (whatever, it's okay - Night is just poorly sequenced/genotyped in the first place so we'll probably end up not even using that population)




### Population Statistics in general:
---

To filter by Minor Allele Frequency when calculating statistics or not?  See: http://csb.scichina.com:8080/kxtbe/EN/abstract/abstract413250.shtml, "Effects of cutoff thresholds for minor allele frequencies on HapMap resolution: A real dataset-based evaluation of the Chinese Han and Tibetan populations"  Basically you should include all the SNPs.

Also see: https://www.biomedcentral.com/content/pdf/1753-6561-3-S7-S41.pdf ;"The effect of minor allele frequency on the likelihood of obtaining false positives"; namely, "These results suggest that removal of low MAF SNPs from analysis due to concerns about inflated false-positive results may not be appropriate"

Also see: http://www.biomedcentral.com/content/pdf/1471-2148-12-94.pdf ; "Uninformative polymorphisms bias genome scans for signatures of selection"; "Markers with a very low minor allele frequency therefore lack the adequate sensitivity to capture the historical signatures of drift and hitchhiking, the key processes in genome scans." 

Also, from an email exchange with Mark Kirkpatrick: 
> I'm working with linkage diseq and Fst at the moment.  The BayeScan manual (Fst) just notes that including low-MAF loci is bad and will skew the results, but doesn't elaborate on why.  I've asked the developer but he's pretty bad about email responses.  And for LD, r^2 is influenced by allele frequencies, so would that matter when comparing pairs of sites with super rare alleles?
"Yes.  Throwing out real data will bias your statistics."




## BayeScan Fst analysis
---
Preparing data:
####

--- **Pairwise comparisons, or all populations considered together?**
Do you consider loci that are outliers when comparing two pops pairwise, or do you apply BayeScan to all five or six or twenty pops together?  For example, if you have southern, northern, eastern and western pops, do you compare S-N, S-W, S-E, N-W, N-E and E-W, or S-N-E-W at once?  Do you lose power for having multiple comparisons?  Must you correct for multiple comparisons?  Are the true signatures of local adaptaion that may be unique to the southern population lost when you compare all four populations instead of just S-N?

This is what Matthieu Foll (BayeScan) says:
> It is true that you lose power by only doing pairwise comparisons. Also, once you will have your 15 pairwise results, what do you do? They are not independent, as populations share some history. So how are you going to assign significance to a given locus? If one outlier pops up in let's say 3 or 4 pairwise comparisons you will be happy? What is the probability that this happens by chance? I know some papers are published with this approach but I personally dislike like.
> It is also true that if you have only a few populations leaving in the specific environment you are interested in, the signal will be diluted when mixed with several others... So no easy solution.

So, indeed, no easy answer, but because our approach is to simply identify "interesting" contigs in the genome, we decided to perform BayeScan on pairs of pops AND five pops together while simply keeping decent records for which BayeScan run the "interesting" label came from.

--- **To filter on MAF or not:**  The BayeScan manual says this:
> ... However, in the extreme case of totally uninformative data, the posterior odds will simply be equal to the prior odds. This is for example the case for uninformative loci such as monomorphic markers or markers with a very low minor allele frequency. We advise people to exclude these markers from the analysis. 

So, here, it seems appropriate to remove MAF < 0.05, at least for BayeScan Fst calcs

---  **MAF in general:**  Definition and considerations

_Definition:_ MAF can either be considered as the frequency of the real minor allele, or simply as the frequency of the ALTERNATIVE/non-reference allele.  I think the latter definition is commonly used, but the alt allele is defined simply through comparison to a genome, which is usually derived from a single individual.  That's biologically and statistically meaningless, if a conventional definition.  The informative/uninformative SNP label I'm trying to apply to my data set is relevant to the actual freq instead, not some arbitrary 'ref'/'alt' status.

Therefore, my script calculates the freq of the TRUE minor allele and only prints out SNPs that are above whatever supplied threshold in at least one pop.  It also requires each locus to be genotyped in at least five individuals per population.  

_Considerations:_

if you are extracting pairs of pops, do you remove a locus if the MAF is 0.03 in one pop BUT is 0.2 in the other pop? (no.) Does it matter if the minor allele is a different allele in each pop?  (yes.)  What is the cutoff threshold for MAF?  How rare does an allele have to be to be considered an uninformative SNP?

Consider this:

pop1 allele "A" | pop2 allele "A" | MAF pop1 | MAF pop2 | notes
| --- | --- | --- | --- | --- |
95% | 97% | 0.05 | 0.03 | This locus is uninformative
95% | 3% | 0.05 | 0.03 | But what about this one?  
80% | 98% | 0.20 | 0.02 | One pop has a too-low MAF but the other doesn't

**my decisions:**  
	*  _Exclude a locus when MAF in both pops is too low_;  "low" may have to be empirically determined, as in the Roesti paper, i.e. "try a bunch of cutoffs and stop when the Fst plots quit changing."
	*  _Include the locus as long as the MAF in ONE pop is high enough_, even if it's below threshold in the other 
	*  _Include the locus even if MAF is below threshold, AS LONG AS the MAF is for a DIFFERENT allele_ (i.e., almost fixed but for opposite/different alleles)
	*  _Exclude a locus is there are fewer than 5 individuals genotyped PER pop_.  Note: one of the criteria we applied to making the multi-pop vcf in the first place was that each SNP had to be genotyped in at least 80% of individuals.  I can imagine a circumstance where a locus is poorly represented or maybe missing entirely in one pop, but genotyped in 100% of individuals from all other pops.  I'd like to exclude these cases, so that's how I coded it.


--- **MAF threshold for filtering:** 
How low is "very low"?  Some lit discards MAF lt 0.25 (Bradbury et al 2010, Atlantic cod paper)!  (seems pretty high to me).  Do you discard MAF lt 0.05 (for 50 inds * 2 chromosomes = 100 alleles; MAF = 0.05 = 5 alleles)?  MAF lt 0.01 (a single heterozygote among 50 individuals)?  Roesti et al 2012 pretty much says that you should look at your own data closely and see when the Fst graphs stop changing, and use that cutoff.  They also use an "n" threshold instead of a MAF percentage, as in, the minor allele must be seen at least 1 or 4 or 10 (etc) times.  

Or, it seems more common to simply discard alleles that are only seen once in the data set.  See: http://hpc.ilri.cgiar.org/training/embo/course_material/genotyping/NGS_and_Genotyping_Rounsley.pdf - "Analyze at the population level, not sample level.  Filter SNPs based on frequency in population; discard SNPs found only in single samples."  Also see: http://rspb.royalsocietypublishing.org/content/277/1701/3725.short - "Parallel adaptive evolution of Atlantic cod on both sides of the Atlantic Ocean in response to temperature" Bradbury et al 2010

This is actually a very good idea, especially for Type II B RAD: we can't remove PCR duplicates due to the nature of the digestion, so it's quite possible that an allele seen only once is indeed genotyped due to an artificial increase in read depth via PCR duplication.  Though, do keep in mind that many rare alleles could still be quite real; see Roesti et al 2012 for a paragraph about how lots of rare variants are totally expected.

--- **What to do with low-MAF loci:** Delete the minor allele, or delete the locus?
In Stacks, for example, if a locus has a too-low MAF, that minor allele gets deleted such that the locus now looks homozygous.  My preference is to simply delete the entire locus.  I can see arguments for either method and am open to a conversation about it, but it seems that the minor allele is (not always but usually) a true biological signal and shouldn't simply be ignored (see Roesti et al 2012).  If that locus is uninformative, then you should just not consider that locus, rather than squint at it a little bit and massage the data such that the minor allele straight-up disappears.  

--- **Do you filter for MAF within or across populations?** 
It seems reasonable to make sure that a particular locus' MAF is gt threshold in at least ONE population of the two compared.  For example, if the threshold is "two alleles"/"must be seen at least twice", do you require those two alleles to be found in a single population, or is it okay if there is one allele in one pop and one allele in the other?  Is it okay if the MAF in pop1 is 0.5 but the MAF in pop2 is 0.03?  (yes, that's a true signal)
I think it's probably best to consider MAF as a within-pop thing.  



Test this shit out:
#### BayeScan - Testing for MAF effects.  Run BayeScan on 30k randomly sampled SNPs from KxO population comparison using different MAF cutoffs: 
1.  All SNPs equally likely to be selected as part of the 30k (no MAF filtering at all)
2.  Exclude singleton SNPs only (as in, if an allele is only present once in EITHER population/a single alternate-allele-having heterozygote exists in any one pop)
3.  Exclude SNPs that are only seen twice (as in, a single diploid alt-allele-having homozygote across populations
4.  Require SNPs to have MAF gt 0.1 in at least one pop
5.  Require SNPs to have MAF gt 0.25 in BOTH pops



1.  Extract pairwise pops based on MAF cutoffs.  
	* `vcf_extract_pops_by_MAF.pl` v. 2 Oct 2013 R. Capper to extract pairwise populations and delete loci based on MAF percentage threshold.
	* ` ` v. 2 Oct 2013 R. Capper to extract pairwise populations and delete loci based on number of times an allele is encountered (remove singleton alleles, alleles seen only twice, etc)
3.  vcf -> genepop
4.  genepop -> bayescan
5.  run bayescan


#### BayeScan - decisions, data prep, the runs, and processing
Two steps:

1.  Extract pairs of pops from ultimate.vcf
2.  Apply MAF > 0.05 filtering to all pairwise pop vcfs, AND to all_five_pops.vcf 

How many SNPs result per file?  Can BayeScan actually run on that number of SNPs?  Seems like BayeScan really can only handle ~30k SNPs per run.  If you have more than 30k SNPs, I suggest splitting them into pieces and running BayeScan independently on each set of SNPs, then combining later.

```

```


## Pi
---
Can you actually calculate pi from RAD data?  Likely not: we are NOT looking at every variant site in the genome, only a subsample.  Unless you had a breakdown of known invariant/variant sites then this won't work; i.e., you'd need to the lengths of each tag and where the SNPs are within those tags, not simply a pretty .vcf file.


## Linkage Disequilibrium
---
See Association Mapping of Kernel Size and Milling Quality in Wheat (Triticum aestivum L.) Cultivars, Flavio Breseghello*,† and Mark E. Sorrells 2006 : http://www.genetics.org/content/172/2/1165.full.pdf

Basically, you can either draw an arbitrary cutoff value at r^2 = 0.1, or you can mark the 95th percentile of the distribution of unlinked loci and intersect that line with the loess regression.

Also see : http://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/

## Tajima's D
---
`vcftools --vcf KxO.vcf --TajimaD 10000`
`mv out.Tajima.D KxO.TajimaD`


### Population Connectivity:
---
## STRUCTURE:
Try to use as many SNPs as possible.  To filter for MAF or not?  Maybe not such a good idea...  But I can't find any decent documentation or "best practices".  Sarah says to not throw out any data; private alleles with a MAF of 0 are still good data!  So, maybe leave them all.  Run "5pops all SNPs" x 2, 3, 4k and "6pops all SNPs" x 2, 3, 4k.  

But, do filter out the outlier loci as identified by BayeScan.  

1.  Convert vcf into genepop  (vcf2genepop.pl)
2.  Delete loci that are implicated as outliers by BayeScan (remember, we did BayeScan several times: once per pairwise comparison, once for the five pops and once for all six pops, so consider/remove all outlying loci) `remove_loci_from_genepop.pl`
3.  Split the genepop into random pieces (note: *SNPs are randomly selected, no attention paid to LD*): `genepop_split.pl six_pops.neutral.genepop (number of chunks) r`
4.  convert .genepop into .structure: `java -jar /work/01408/rlc2489/RoxyTools/PGDSpider_2.0.4.0/PGDSpider2-cli.jar \` `-inputfile five_pops.neutral.genepop -inputformat GENEPOP \` `-outputfile five_pops.neutral.structure -outputformat STRUCTURE \ -spid genepop2structure.spid`
5.  make mainparams, extraparams files using the desktop version of STRUCTURE
6.  run STRUCTURE:  Note, here, try running to 500,000 runs (50,000 burn-in) and check for convergence.  1M runs will likely not finish within 24hrs.


#### STRUCTURE Stats:
---

##### Time for running five_pops STRUCTURE, k=2
number of SNPs | raw rate | per hr | per 1M |
| ------------ | ------------- | ------ | ---- |
30k SNPs: | 1000 runs/2hrs | 500 runs/hr | 1M runs/83days |
10k SNPs: | 3000 runs/6mins | 30,000 runs/hr | 1M runs/33hrs |
5k SNPs: | 4500 runs/5mins | 54,000 runs/hr | 1M runs/18.5hrs |
1k SNPs: | 28400 runs/6mins |  284,000 runs/hr | 1M runs/3.5hrs  |
0.5k SNPs: | 140,200 runs/15mins | 560,800 runs/hr | 1M runs/1.8hrs |
0.1k SNPs: | 218,000 runs/5mins | 2,616,000 runs/hr | 1M runs/0.4hrs |
		
##### Time for running FIVE_pops STRUCTURE, 5k SNPs, various k's
| K | raw rate | per hr | per 1M |
| --- | --- | --- | --- |
k=2: | 13,000 runs/14mins | 55,714 runs/hr | 1M runs/17.9hrs
k=3: |  5,000 runs/7mins | 42,857 runs/hr | 1M runs/23.3hrs
k=4: |  4,800 runs/7mins | 41,143 runs/hr | 1M runs/24.3hrs
k=5: |  5,600 runs/10mins | 33,600 runs/hr | 1M runs/29.7hrs
		
##### Time for running SIX_pops STRUCTURE, 5k SNPs, various k's
| K | raw rate | per hr | per 1M |
| --- | --- | --- | --- |
k=2: | 6,400 runs/10mins | 38,400 runs/hr | 1M runs/26.0hrs |
k=3: | 5,600 runs/10mins | 33,600 runs/hr | 1M runs/29.8hrs |
k=4: | 6,500 runs/13mins | 30,000 runs/hr | 1M runs/33.3hrs |
k=5: | 8,100 runs/17mins | 28,588 runs/hr | 1M runs/35.0hrs |



