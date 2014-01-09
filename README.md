RAD for genetic analysis of natural populations
===
Type II B RAD-seq processing pipeline: raw reads through population analysis

Type II B RAD uses special restriction enzymes to excise 36b DNA "tags" distributed evenly across the genome of (theoretically) any organism.  

Here, we used the enzyme BcgI, which recognizes the site ...∇(N10)CGA(N6)TGC(N12)∇..., to digest the genome of the scleractinian coral *Acropora millepora*.  

Note: BcgI is methylation-**sensitive**; we're not sure what kind of bias, if any, this introduces, but it is definitely something to consider.


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

---  **MAF in general:**  False positives, definition and considerations

_Effect on false positives:_  See: http://www.nature.com.ezproxy.lib.utexas.edu/nrg/journal/v6/n2/pdf/nrg1521.pdf, "GENOME-WIDE ASSOCIATION STUDIES FOR COMMON DISEASES AND COMPLEX TRAITS"; Box 4.

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
+ **Exclude a locus when MAF in both pops is too low**;  "low" may have to be empirically determined, as in the Roesti paper, i.e. "try a bunch of cutoffs and stop when the Fst plots quit changing."  Note: make sure that the MAFs are relative to the same allele here.  It would be a shame to delete a locus at which pop1 is 95% 'A' and pop2 is 5% 'A'.
+  **Include the locus as long as the MAF in ONE pop is high enough**, even if it's below threshold in the other 
+  **Include the locus even if MAF is below threshold, AS LONG AS the MAF is for a DIFFERENT allele** (i.e., almost fixed but for opposite/different alleles)
+  **Exclude a locus is there are fewer than 5 individuals genotyped PER pop**.  Note: one of the criteria we applied to making the multi-pop vcf in the first place was that each SNP had to be genotyped in at least 80% of individuals.  I can imagine a circumstance where a locus is poorly represented or maybe missing entirely in one pop, but genotyped in 100% of individuals from all other pops.  I'd like to exclude these cases, so that's how I coded it.



--- **MAF threshold for filtering:** 
How low is "very low"?  Some lit discards MAF lt 0.25 (Bradbury et al 2010, Atlantic cod paper)!  (seems pretty high to me).  Do you discard MAF lt 0.05 (for 50 inds * 2 chromosomes = 100 alleles; MAF = 0.05 = 5 alleles)?  MAF lt 0.01 (a single heterozygote among 50 individuals)?  Roesti et al 2012 pretty much says that you should look at your own data closely and see when the Fst graphs stop changing, and use that cutoff.  They also use an "n" threshold instead of a MAF percentage, as in, the minor allele must be seen at least 1 or 4 or 10 (etc) times.  

Or, it seems more common to simply discard alleles that are only seen once in the data set.  See: http://hpc.ilri.cgiar.org/training/embo/course_material/genotyping/NGS_and_Genotyping_Rounsley.pdf - "Analyze at the population level, not sample level.  Filter SNPs based on frequency in population; discard SNPs found only in single samples."  Also see: http://rspb.royalsocietypublishing.org/content/277/1701/3725.short - "Parallel adaptive evolution of Atlantic cod on both sides of the Atlantic Ocean in response to temperature" Bradbury et al 2010

This is actually a very good idea, especially for Type II B RAD: we can't remove PCR duplicates due to the nature of the digestion, so it's quite possible that an allele seen only once is indeed genotyped due to an artificial increase in read depth via PCR duplication.  Though, do keep in mind that many rare alleles could still be quite real; see Roesti et al 2012 for a paragraph about how lots of rare variants are totally expected.

--- **What to do with low-MAF loci:** Delete the minor allele, or delete the locus?
In Stacks, for example, if a locus has a too-low MAF, that minor allele gets deleted such that the locus now looks homozygous.  My preference is to simply delete the entire locus.  I can see arguments for either method and am open to a conversation about it, but it seems that the minor allele is (not always but usually) a true biological signal and shouldn't simply be ignored (see Roesti et al 2012).  If that locus is uninformative, then you should just not consider that locus, rather than squint at it a little bit and massage the data such that the minor allele straight-up disappears.  

--- **Do you filter for MAF within or across populations?** 
It seems reasonable to make sure that a particular locus' MAF is gt threshold in at least ONE population of the two compared.  For example, if the threshold is "two alleles"/"must be seen at least twice", do you require those two alleles to be found in a single population, or is it okay if there is one allele in one pop and one allele in the other?  Is it okay if the MAF in pop1 is 0.5 but the MAF in pop2 is 0.03?  (yes, that's a true signal)
I think it's probably best to consider MAF as a within-pop thing.  

--- **Coding issues and solutions**

There are three (four) methods for calculating MAF:

1.  By raw counts, across populations: `vcf_extract_by_MAF_counts.pl`, v. 9 Oct 2013.  This is the easiest way to do it for sure, but counting rare alleles does not account for changing number of genotyped individuals at a particular locus.  For example, maybe you have 20 dudes but only 10 are genotyped at a certain locus.  Is it better to insist on seeing 2 alleles no matter if you have 20 dudes (MAF = 0.05) or if you have 10 (MAF = 0.1)?  
2.  By raw counts, within populations.  This is tricky to script, at least for me, and likely not necessary.  This introduces a whole lot of conditions that may be stupid.  For example, 
	+ monomorphs get thrown out.  
 	+ if the rare allele is below threshold in BOTH pops, clearly, exclude it.  
 	+ If it's below in one pop but above in the other, then maybe keep that locus.
 	+ If it's below threshold in BOTH pops but for DIFFERENT alleles (i.e., just about fixed for opp alleles) then keep that.
 	+ If it's totally fixed in one but fixed for the opposite allele in the other, then definitely keep it even though the MAF of one pop is 0.  
3.  By frequency/pecentage, across populations.  I like this method better than the by-counts way because you can make MAF = some threshold regardless of how many individuals are genotyped, but I don't have a strong opinion about it.
4.  By frequency/percentage, within populations: `vcf_extract_by_MAF_percentage.pl`.  This is how I did it first.  However, at this point, I can't remember why I thought within-pop MAF had any advantage over across-pop MAF.  At best, it's probably going to give you the same results as across-pop calculations, but with more room for screwing up the script.  For example, you may have 90 'A' alleles and 10 'a' alleles across 2 pops.  If your MAF cutoff is 0.06, you would keep that locus if you look across pops (10/100 = 0.1, > 0.06).  But if you look within pops, the 'A' allele could be distributed equally (MAF pop1 = 0.1, MAF pop2 = 0.1).  But, maybe all 'a' alleles are in one pop (MAF pop1 = 0, MAF pop2 = 0.2).   I would still keep that site.  Is there a situation where the 'a' allele could be distributed unequally such that I would discard that site when looking within-pops instead of across-pops?  As long as one pop is above threshold, then I tend to keep that locus.  Maybe there are some cases where different numbers of individuals/contributed alleles would affect the cutoff such that the distribution of minor alleles is below threshold in both pops but above overall but I can't even come up with an illustration of this at the moment.  Whatever, I already did the work, and here was my logic behind it:

example if threshold = 0.1

pop1 allele A | pop2 allele A | MAF pop1 | MAF pop2 | notes
| ---- | ---- | ---- | ----- | ----- |
    95%       |      97%      |   0.05   |   0.03   | Exclude; too low in both pops
    95%       |       3%      |   0.05   |   0.03   | Include; low MAF but different alleles, i.e. nearly-private alleles
    80%       |      98%      |   0.20   |   0.02   | Include; low MAF but only in a single pop
   100%       |     100%      |   1.00   |   1.00   | Exclude; monomorphic
   100%       |       0%      |   1.00   |   0.00   | Include; fixed for different alleles


> I scripted `vcf_extract_by_MAF_percentage.pl` to calculate MAF within pops instead of across them.  This flexibility means that it's possible for the minor alleles to be different between pop1 and pop2; for example, pop1 may be homozygous reference (allele 0, frequency = 1) while pop2 is nearly-but-not-quite homozygous for the alternate allele (allele 0, frequency = 0.05).  This is actually the kind of signal we're looking for and we'd be remiss to throw it out.  However, a word of caution about coding this type of thing:  

> My script calculates MAF by counting the number of times each allele is seen **within** a pop, finding which allele is rarer, then dividing that rare allele's count by twice the number of individuals with genotyped data at that SNP.  So, if you're a homozygote, there is only one allele seen in the data, which is then considered the 'rare' allele, resulting in the MAF being calculated as 1 and not as 0.  If pop2 is het and the rare allele is the SAME as the homozyg allele in pop1, then everything is fine.  However, if pop2 is het and the rare allele is DIFFERENT than the homozyg allele in pop1, then the true MAF of pop1 is 0 and not 1.

> Imagine that pop1 has the following allelic composition: {0, 0, 0, 0}.  The frequency of allele 0 is 1; the freq of allele 1 is 0.  MAF gets calculated as 1.
Now imagine that pop2 looks like this: {0, 0, 0, 1}.  Freq{0} = 0.75, freq{1} = 0.25.  MAF gets calculated as 0.25.
It would be erroneous to compare pop1's MAF=1 to pop2's MAF=0.25.  The MAF of pop1, relative to the MAF of pop2, is actually 0.  

    






Test out how MAF affects Fst using BayeScan AND vcftools fst.
---------
#### BayeScan - Testing for MAF effects.  Run BayeScan on SNPs from KxO population comparison (59 individuals = 118 alleles) using different MAF cutoffs (`vcf_extract_by_MAF_counts.pl`).  Question:  How does varying MAF cutoff affect Fst distribution?  

1.  All non-monomorphic SNPs (51910 SNPs)
1.  Exclude singleton SNPs only (as in, if an allele is only present once in EITHER population/a single alternate-allele-having heterozygote exists in any one pop) (MAF ~= 0.01)
3.  Exclude SNPs that are seen 3 or fewer times (MAF ~= 0.025)
4.  Require SNPs to have MAF gt 0.1 in at least one pop (see each allele 12 or more times across both pops)
5.  Require SNPs to have MAF gt 0.25 in BOTH pops (see each allele 31 times or more across both pops)

---

1.  vcf -> genepop: `vcf2genepop.pl vcf="KxO.maf_gt_0.01.vcf" pops=K,O > KxO.0.01.genepop`
4.  genepop -> bayescan: `java -jar /work/01408/rlc2489/RoxyTools/PGDSpider_2.0.4.0/PGDSpider2-cli.jar\` `-inputfile AxK.genepop -inputformat GENEPOP -outputfile AxK.bayescan\` `-outputformat GEST_BAYE_SCAN -spid genepop2bayescan.spid`
5.  run bayescan with threads option: `bayescan_2.1 KxO.0.01.bayescan -threads 24`

#### vcftools weir and cockerham Fst:  How does varying MAF cutoff affect Fst distribution/profiles?

Answer:  It doesn't really.  Each locus' Fst is calculated independently from the others, so the values themselves don't change.  What does change is the NUMBER of SNPs included.  If a lot of your low-MAF SNPs are the same freq across pops (like, pop1 = 0.01, pop2 = 0.00), then there's no variation!  And therefore Fst is low!  So you're cutting out all those guys that bring the baseline down, which usually would bring the baseline up and/or resolve signals and peaks better, or at least make them higher/clearer.

For example: see fig 3 in Roesti et al.  I can't figure out how to post pictures in markdown, but basically the figure is Fst vs chr position, showing that the Fst profile is much lower and has fewer peaks/lower signal when INCLUDING low-MAF SNPs, than when EXCLUDING them.  It makes sense that you might see a lower profile when including more, lower SNPs especially if you use sliding window/smoothing or whatever, as you probably should.

#### vcftools weir and cockerham Fst:  How does varying MAF cutoff affect global Fst values?

1.  Extract SNPs by MAF cutoff using `vcf_extract_by_MAF_counts.pl` as above, five conditions: 
	+  no MAF filter (51910 SNPs), 
	+  MAF gt 0.01 (40906 SNPs), 
	+  MAF gt 0.025 (30458 SNPs), 
	+  MAF gt 0.10 (16438 SNPs), 
	+  MAF gt 0.25 (6500 SNPs)
2.  run vcftools `vcftools --vcf KxO.all_maf.vcf --weir-fst-pop K.pop --weir-fst-pop O.pop`
3.  output = weighted Fst; Graph Fst vs. MAF; Graph Fst vs. number of SNPs
	+  however, this may very well be confounded by number of SNPs!  So, to disentangle this (just for fun), extract 50k, 41k, 30k, 16k and 6.5k SNPs randomly without regard to MAF from the KxO non-monomorphic SNP set, run vcftools again, graph.
4.  Extract random SNPs:








#### BayeScan - decisions, data prep, the runs, and processing
Two steps:

1.  Extract pairs of pops from ultimate.vcf
2.  Apply MAF > 0.05 filtering to all pairwise pop vcfs, AND to all_five_pops.vcf 

How many SNPs result per file?  Can BayeScan actually run on that number of SNPs?  Seems like BayeScan really can only handle ~30k SNPs per run.  If you have more than 30k SNPs, I suggest splitting them into pieces and running BayeScan independently on each set of SNPs, then combining later.



Fst Plotting:
---
+  Local polynomial fitting?
+  sliding windows?

```

```


## Pi
---
Can you actually calculate pi from RAD data?  We are NOT looking at every variant site in the genome, only a subsample.  Unless you had a breakdown of known invariant/variant sites then this won't work; i.e., you'd need to the lengths of each tag and where the SNPs are within those tags, not simply a pretty .vcf file.

However, if you get GATK to output a list of EVERY gENOTYPED BASE including invariants, then yes, you CAN calculate pi.  See my script vcf2pi.pl.

## Linkage Disequilibrium
---
See Association Mapping of Kernel Size and Milling Quality in Wheat (Triticum aestivum L.) Cultivars, Flavio Breseghello*,† and Mark E. Sorrells 2006 : http://www.genetics.org/content/172/2/1165.full.pdf

Basically, you can either draw an arbitrary cutoff value at r^2 = 0.1, or you can mark the 95th percentile of the distribution of unlinked loci and intersect that line with the loess regression.

Also see : http://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/

Also see: "fig3: Plots of linkage disequilibrium (LD) values (r2) against genetic distance (cM) between pairs of markers in multiple classes of cultivated tomato. All possible pair-wise combinations of markers on the same chromosome were plotted to visualize LD decay within chromosomes over the entire genome. The r2 values were calculated separately for processing and fresh market cultivars (B and C, respectively) as well as processing, fresh market, and vintage cultivar classes combined (A). Curves were fit for each plot by second-degree LOESS. The horizontal dotted lines indicate the baseline r2 values based on the 95th percentile of the distribution of unlinked r2 values (black) and the fixed r2 value of 0.1 (grey).

Mentions: The decay of LD over genetic distance was investigated by plotting pair-wise r2 values against the distance (cM) between markers on the same chromosome (Fig. 3). A smooth line was fit to the data using second-degree locally weighted scatterplot smoothing (LOESS; Breseghello and Sorrells, 2006) as implemented in SAS. To describe the relationship between LD decay and genetic distance, two methods of establishing baseline r2 values were investigated. Critical values of r2 were based on a fixed value of 0.1 (Nordborg et al., 2002; Palaisa et al., 2003; Remington et al., 2001) and from the parametric 95th percentile of the distribution of the unlinked markers (Breseghello and Sorrells, 2006). The relationship between these baseline r2 values and genetic distance was determined using the LOESS curve and a 1 cM moving means approach. For the LOESS estimation of LD decay, genetic distance was estimated as the point where the LOESS curve first crosses the baseline r2 value. For the moving means approach, the distance between linked markers was used to divide marker pairs into bins of 1 cM. Markers separated by 0–0.9 cM were placed in the first bin, marker distances from 1–1.9 were in the second bin, etc. The mean of the r2 values within each bin was calculated and LD decay was estimated as the first bin where the baseline r2 value was lower than the bin mean."  from : http://openi.nlm.nih.gov/detailedresult.php?img=3060673_jexboterq367f03_ht&req=4


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



