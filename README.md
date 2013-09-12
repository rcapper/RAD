RAD for genetic analysis of natural populations
===
Type II B RAD-seq processing pipeline: raw reads through population analysis

Type II B RAD uses special restriction enzymes to excise 36b DNA "tags" distributed evenly across the genome of (theoretically) any organism.  

Here, we used the enzyme BcgI, which recognizes the site ...∇(N10)CGA(N6)TGC(N12)∇..., to digest the genome of the scleractinian coral *Acropora millepora*.  


## Trimming Raw Reads
---
There are several schools of thought about trimming strategies.  Of course, you must trim to remove low-quality bases, but should you require each read to contain the restriction site or not?  Reads that do not contain the restriction site may be produced by either sequencing errors or enzyme fidelity errors.  Since we don't know much about this strange Type II B enzyme, it's within the realm of possibility that the enzyme may cut at slightly non-canonical sites.  But, do we then even care about these sites?  I can see reasons to either require or not require the site to be present in each read, but at the end of the day we decided to require it.  





