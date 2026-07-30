# Manuscript outline and Methods

**Working title options**

1. Tissue-specific regulatory architecture distinguishes repressed from induced genes in the renal cortex of cats with chronic kidney disease
2. Cis-regulatory annotation of the feline renal cortex and its application to a spontaneous chronic kidney disease transcriptome
3. The cortical transcriptomic signature of feline chronic kidney disease reflects a shift from kidney-specific to ubiquitous regulatory programs

Option 3 states the principal finding and is preferred. Option 2 is the safer fallback if reviewers judge the finding descriptive.

---

## Section plan

**Abstract.** Structured or single paragraph depending on journal. Order: feline CKD is common and its cortical transcriptome has been described but not interpreted against regulatory architecture; a cis-regulatory map of seven feline tissues exists but has never been applied to a disease transcriptome; we integrated the two; the repressed and induced programs differ systematically in the tissue breadth of their enhancers; three independent lines indicate that the dominant axis in bulk cortex is cellular composition; we report an assembly-harmonization hazard specific to these datasets. State explicitly that the study is exploratory and hypothesis-generating.

**Introduction.** Three paragraphs. (1) Feline CKD, prevalence in aged cats, and the existing transcriptomic and multi-omic description. (2) What bulk tissue transcriptomes cannot resolve, and the fact that prior work has raised cellular composition as a caveat without measuring it. (3) The feline cis-regulatory atlas as an unused resource, and the specific questions addressed here.

**Results.** Six subsections following the analysis order: regulatory annotation and its validation; differential expression and reproduction of the source study; coexpression structure is not reproducible; the dominant axis correlates with tubular content; motif enrichment and the ERR observation with its architecture controls; tissue breadth of regulatory elements distinguishes repressed from induced genes.

**Discussion.** Principal finding; relation to the two prior analyses of this dataset; why the composition interpretation is favoured; what the ERR observation does and does not support; the assembly hazard as a practical warning; limitations; implications for cell-resolved study designs.

**Conclusions.** Brief, on the theme rather than on the study.

---

## Materials and Methods

### Data sources

Renal transcriptomic data were obtained from the multi-omic study of spontaneous feline CKD deposited under GEO accession GSE303653 (Li et al., *Communications Biology*, 2025). Gene-level expression matrices for renal cortex and medulla were retrieved from Supplementary Data 3 and 4 of that publication. Although the source text describes these files as raw read counts, the deposited values are log2-transformed counts per million normalised by the trimmed mean of M-values, as evidenced by non-integer and negative entries; all downstream modelling was chosen accordingly. The cortical matrix comprised 15,149 genes across 21 samples (6 control, 8 CKD IRIS stage 1/2, 7 CKD IRIS stage 3/4) and the medullary matrix 15,267 genes across 18 samples (6, 7 and 5 respectively). Sample group and sex were read from the header rows of the same files. Age was not used as a covariate because the deposited age row contained values incompatible with the cohort description in the source publication, indicating a formatting error.

Chromatin data were obtained from GEO accession GSE182952, comprising replicated peak sets for H3K4me3, H3K27ac and CTCF in seven tissues from adult male domestic shorthair cats, with two biological replicates per tissue. Sample metadata identify the renal samples as renal cortex. Peaks had been called with MACS2 following ENCODE3 specifications for histone modifications, retaining peaks from the pooled replicate set that overlapped by at least half their length with peaks from both biological replicates. Peak sets are not filtered against a blacklist because none exists for *Felis catus*.

Genome annotation and sequence were obtained from the NCBI RefSeq assembly GCF_018350175.1 (F.catus_Fca126_mat1.0), NCBI *Felis catus* Annotation Release 105, together with the corresponding assembly report.

### Assembly harmonisation

The chromatin peaks in GSE182952 were called against GCA_018350175.1, whereas the transcriptomic data were quantified against Felis_catus_9.0. Because the two assemblies use identical chromosome names, naive harmonisation by chromosome name retains 99.96% of peaks without warning while producing coordinate assignments that are effectively random. In our data this manifested as promoter overlap at 1.55-fold the permuted expectation for H3K4me3, where a genuine promoter mark is expected to exceed tenfold. All coordinate-dependent analyses were therefore performed against F.catus_Fca126_mat1.0. Peak chromosome names were translated in two steps: the `chr` prefix was removed, and the resulting name was mapped to its RefSeq accession using the assembly report, matching against the Sequence-Name, Assigned-Molecule and UCSC-style-name fields. Analyses were restricted to chromosome-level sequences, defined as those whose maximum annotated gene coordinate exceeded 10 Mb.

Expression data were joined to the annotation by gene symbol rather than by coordinate. This is valid because the expression matrix is indexed by symbol and coordinates are used only to position genes relative to peaks; it means the two layers derive from different assemblies, which we state explicitly as a design choice. Symbol correspondence was 99.87% for cortex and 99.88% for medulla.

We recommend that any reuse of these two datasets include a permutation-based positive control, as described below, since the failure mode is silent.

### Differential expression

Because the deposited values are log2 counts per million rather than integer counts, negative-binomial frameworks are not applicable and mean-variance modelling by `voom` cannot be used. Differential expression was performed with `limma` using empirical Bayes moderation with an intensity-dependent trend and robust hyperparameter estimation. The design was `~ 0 + group + sex`; sex was retained because the compartments are not sex-balanced. Contrasts were CKD 1/2 versus control, CKD 3/4 versus control, and CKD 3/4 versus CKD 1/2, each fitted separately within compartment. Multiple testing was controlled by the Benjamini–Hochberg procedure. The primary reporting threshold was a false discovery rate below 0.10, matching the threshold used in the source publication.

Counts of differentially expressed genes were compared against the source study as a reproduction check. At a false discovery rate below 0.10 we obtained 3, 5,087, 1,908 and 2,821 genes for cortical early, cortical late, medullary early and medullary late contrasts respectively, against 6, 4,616, 1,995 and 2,484 reported by the source study using edgeR. The CKD 3/4 versus CKD 1/2 contrast yielded no genes in cortex and eight in medulla and was not carried forward; this decision was made after inspecting the result and is reported as such.

### Regulatory annotation

Gene models were taken from Annotation Release 105, restricted to features with an assigned gene symbol. Promoter windows were defined as 2,000 bp upstream and 500 bp downstream of the annotated transcription start site, oriented by strand. Because no blacklist exists for this species, the widest 0.5% of peaks for each mark were discarded as a conservative proxy for mapping artefacts.

Regulatory domains were constructed for each gene as the interval bounded by the nearest CTCF peak midpoint upstream and downstream of the transcription start site, following the approach used by the Functional Annotation of Animal Genomes consortium, and capped at 250 kb on each side to avoid pathological domains in CTCF-sparse regions. The resulting domains had a median width of 39 kb.

Genes were assigned to one of four mutually exclusive classes in hierarchical order: active promoter, where the promoter window overlapped both an H3K4me3 and an H3K27ac peak; primed promoter, where it overlapped H3K4me3 only; distal enhancer only, where no promoter mark was present but at least one H3K27ac peak fell within the regulatory domain outside any promoter window; and no active mark. Per-mark peak counts within promoter windows and within domains were retained as continuous variables so that classes can be redefined without reprocessing.

### Validation of the regulatory annotation

Coordinate correspondence between the chromatin and annotation layers was verified by permutation. Peaks were shuffled within chromosome, preserving peak number and width, over 100 iterations, and the number of genes with at least one promoter peak was recorded. Observed counts were divided by the mean of the permuted distribution. The analysis pipeline aborts if H3K4me3 promoter enrichment falls below fivefold. Observed enrichment was 8.59-fold for H3K4me3, 5.09-fold for CTCF and 3.92-fold for H3K27ac, with 49.2%, 40.5% and 48.7% of genes carrying a promoter peak of the respective mark.

### Coexpression structure

Differentially expressed protein-coding genes were partitioned by direction and clustered separately, since combining directions groups genes by sign rather than by co-regulation. Gene-gene distances were computed as one minus the Spearman correlation across samples, followed by average-linkage hierarchical clustering and adaptive branch cutting with a minimum module size of 30 genes.

Module reproducibility was assessed by repeatedly removing 20% of animals at random, reclustering, and computing the maximum Jaccard index between each original module and any module in the resampled solution, over 200 iterations. Modules with a mean Jaccard index of at least 0.6 were retained. Module eigengenes were defined as the first principal component of the scaled expression of member genes, with sign fixed to the mean loading so that the eigengene increases when member genes increase. Association with disease group was tested by Kruskal–Wallis.

Cell-composition scores were computed as the mean of per-gene z-scores across marker sets defined a priori and specific to compartment: proximal tubule markers for cortex, and thick ascending limb, collecting duct, intercalated cell and thin limb markers combined for medulla, together with shared immune and fibroblast marker sets. Marker set membership and the number of markers detected in each matrix are reported in the supplement.

### Motif enrichment

Regulatory regions were defined from peak summits, using the summit offset recorded in the narrowPeak files, and standardised to 300 bp centred on the summit so that sequence length does not differ between comparison sets. Foreground regions were the promoter-overlapping H3K4me3 peaks, or the domain-contained non-promoter H3K27ac peaks, of the genes in each set. Background regions were the corresponding elements of expressed protein-coding genes not differentially expressed in that contrast, sampled to at most three times the foreground and capped at 4,500 sequences; foreground was capped at 1,500. Sequences containing ambiguous bases were discarded and each sequence was given a unique coordinate-derived name.

Enrichment was computed with `monaLisa` using binned motif enrichment with the background drawn from the opposing bin and GC bins spanning 0 to 1. Position weight matrices were taken from the JASPAR2020 CORE vertebrate collection and restricted to motifs whose corresponding gene was present in the expression matrix of the compartment, leaving 401 of 633 matrices. Enriched motifs were cross-referenced against the differential expression result for the corresponding transcription factor, and a motif was considered directionally coherent only when the factor changed in the direction compatible with the regulated set.

Each analysis was repeated with cell-type marker genes removed from the foreground set as a sensitivity analysis. Because foreground sequences were subsampled independently in each run, the two runs are not paired and this sensitivity analysis should be interpreted with that limitation in mind.

### Architecture controls for the oestrogen-related receptor observation

Regions were scanned for individual matrices using `findMotifHits` with a minimum score of 10, and hits were mapped back to genes through the region-to-promoter assignment. CTCF and SP1 matrices were scanned in parallel as negative controls. Functional enrichment of motif-positive genes was computed with `g:Profiler` against human orthologs using the g:SCS correction, with the background restricted to the remaining genes of the same expression class rather than to the transcriptome, so that the test addresses specificity within the class.

Mitochondrial annotation was defined as membership of GO:0005739 including descendant terms, retrieved through `org.Hs.eg.db`. Two logistic models were fitted over all expressed protein-coding genes with at least one scanned promoter region. The first modelled mitochondrial annotation as a function of motif presence, expression state and their interaction, testing whether the motif-function association differs between states. The second modelled repression as a function of motif presence, mitochondrial annotation and the logarithm of the number of scanned regions per gene, the last included because genes with more promoter sequence scanned have a higher prior probability of containing any motif.

### Tissue breadth of regulatory elements

Replicated peak sets for all seven tissues were harmonised to the annotation as described above. For each renal element, presence was recorded in each tissue by overlap, giving a tissue breadth between one and seven. Elements were classified as kidney-specific (one tissue), restricted (two to three), shared (four to five) or ubiquitous (six to seven). This measure derives entirely from chromatin in healthy animals and is independent of the disease expression matrix.

Distal elements were defined as renal H3K27ac peaks not overlapping any annotated promoter window and were assigned to genes by overlap with the gene's CTCF-delimited domain. For each gene we computed the mean tissue breadth of its distal elements, the fraction of those elements that were kidney-specific, the mean fold-change signal of those elements, and their number.

Genes were classified as repressed, induced or unchanged from the cortical CKD 3/4 versus control contrast. Association between enhancer tissue specificity and expression state was tested by Kruskal–Wallis followed by pairwise Wilcoxon tests, and by logistic regression of expression state on the kidney-specific fraction, adjusted for the logarithm of the number of enhancers. Because peak detection sensitivity differs between tissues and weaker peaks are less likely to be detected elsewhere, the models were refitted with mean enhancer signal as an additional covariate; both adjusted and unadjusted estimates are reported.

### Study design and reporting

This is an exploratory, hypothesis-generating secondary analysis of two publicly deposited datasets. It was not conducted under a prespecified analysis plan, and several analytical decisions were made in response to intermediate results: the CKD 3/4 versus CKD 1/2 contrast was dropped after it returned almost no genes; the primary compartment was selected after examining differential expression in both; and coexpression modules were abandoned as an analytical unit after they failed reproducibility assessment. These decisions are stated here rather than presented as prespecified. Negative results are reported in the main text with the same detail as positive ones.

### Software

All analyses were performed in R version [fill from `sessionInfo()`]. Principal packages were `limma` [version], `rtracklayer` [version], `GenomicRanges` [version], `Rsamtools` [version], `Biostrings` [version], `dynamicTreeCut` [version], `monaLisa` version 1.18.0, `TFBSTools` [version], `JASPAR2020` [version], `org.Hs.eg.db` [version], `gprofiler2` [version], `ggplot2` [version] and `patchwork` [version]. Random seeds were fixed at the top of each module and are recorded in the deposited code.

### Data and code availability

All primary data are publicly available under GEO accessions GSE303653 and GSE182952 and NCBI assembly GCF_018350175.1. Analysis code is organised as sequentially numbered modules and is available at [repository URL], archived at Zenodo under DOI [DOI]. Intermediate outputs required to reproduce every figure and table are included in the archive.

### Use of artificial intelligence

Claude (Anthropic) was used throughout this work. Its role extended beyond language editing and is described in full in the interest of transparency.

The model was used to draft and debug the analysis code, to propose and critique analytical designs, to identify confounders and specify the control analyses addressing them, to interpret intermediate results, and to draft and revise the text of this manuscript. Several methodological decisions reported here originated in that exchange, including the permutation-based validation of coordinate correspondence, the resampling assessment of module reproducibility, the substitution of fixed distance windows by CTCF-delimited domains, and the tissue-breadth analysis that produced the principal finding.

All code was executed by the author, and all numerical results reported here were produced by that execution rather than by the model. Every analytical claim was checked by the author against the primary output. Errors introduced by the model during development were identified and corrected during the work, including an incorrect assumption about the reference assembly of the chromatin data, which the permutation control detected.

The model is not an author and does not meet authorship criteria, as it cannot take responsibility for the work. The author takes full responsibility for the content of this manuscript, including all sections in which artificial intelligence was used.

---

## Notes before drafting the remaining sections

Fill software versions from `sessionInfo()` rather than from memory.

Confirm the target journal's policy on artificial intelligence disclosure and adjust the placement of that statement accordingly; some journals require it in the cover letter, some in the Methods, some in a dedicated declaration after the Discussion.

The upper-bound caveat on the kidney-specific fraction belongs in the Discussion with the per-tissue peak counts alongside it: the renal H3K27ac set is the largest of the seven at 81,603 peaks against 35,524 for small intestine, so the 19.58% kidney-specific figure is an upper bound.

State in the Results that the odds ratio of 2.24 compares a gene whose distal elements are all kidney-specific against one with none, and that the median kidney-specific fraction is zero in all three groups, so the effect resides in a minority of genes.
