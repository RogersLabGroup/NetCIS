# Outline for NetCIS

## Tools

- cutadapt
- bowtie2
- bam2bed?
- python 3.10
- networkx
- cytoscape
- GO enrichment
- GSEA?

## Preprocess data

- cutadapt to trim adapter and SB tags
- map reads to reference with bowtie2
- only keep reads that start with "TA"
- output as bam files
- change bam files to bed files, one insertion site per line? How to do this
- strand orientation of insertions
  - positive strand
    - right read maps to positive
    - left read maps to negative
  - negative strand
    - right read maps to negative
    - left read maps to positive

## Read in data

- use pandas to read in tabular data of SB insertions
  - convert to numpy data structures if needed

- have a file that contains the regions that are always hit by the SB we used
  - be up front about what is being removed from the anaylsis
    - provide a way for others to add/change this

- optionally choose what chromosomes to exlude to prevent infalted local hopping insertion events

- database of targeted genes to compare to

## QC checks for input data

- BamQC for quality control of BAM files? <https://github.com/s-andrews/BamQC>

- To do first
- visual representations of a sliding window for the number of inseritons events per sample for the genome
  - try a variety of window sizes, or look at KC-RBM to see how they did it
  - yaxis is chromosomes/bp, and xaxis is the samples. Compare cases to controls for a sample to see what this looks.
- <https://transit.readthedocs.io/en/latest/tpp.html#tpp-statistics>
- <https://transit.readthedocs.io/en/latest/transit_features.html>
- per chromosome, then whole genome, then whole genome minus local hopping chromosomes

## Construct network (using pseudo code from graph framework article)

- be wary of going too far with this. keep it simple

- input connection threshold: default 50 kb
  - default 50 kilobases from graph based framework paper
  - defualt 50 kilobases to expand promoter region from IAS

- create insertion site network
  - Split insertion sites (IS) by chromosome
  - Order IS for each chromosom
  - for each IS(i)
    - add node(i) as an IS location into the network
      - for each node (j) in the network
        - if distance between node(i) and node(j) is less than threshold
          - add edge(ij) with a weight of the distance (or inverse?) to the network

## Analyze network

- Louvain algorithm, graph coarsening, or modularity metric for non-overlapping sub-graphs?

- export network into CytoScape and analyze there OR use networkx for analysis

- for each connected subgraph in the graph
  - if subgraph is not a random network
    - add subgraph to list of non-random CIS

- Explore non-random CIS

- Obtain statistics
  - spatial significance of CIS
    - Poisson or ZINB (zero inflated negative binomial) distribution?
    - Difference to a random distribution?
  - significance of difference between case and controls per CIS and per gene
  - use the controls to determine what the distributionn should look like at a CIS

- Compare to a SB database of known targeted genes
  - candidate cancer gene database from Tim Starr
  - Sleeping Beauty Cancer Driver Database
  - another unselected database from Adam or create one ourselves?

- predict CIS effects?
  - skew, kurtosis, transposon orientation

- Output results

## Thoughts to consider

- saturation of SB screens

- providing in depth statistical analysis of data before running our methods on it
  - look at QC methods from TRANSIT

- Is there a way to provide a case control analysis and one without?

- Use appropriate statistical methods and be very clear about what we are doing
  - rank based methods are preferred
    - Wilcoxon
    - Mann-Whiteney

- Somehow get a p-value of each CIS (look at dist of control)
  - Is there a way to tie-back to the QC for a CIS region to help determine significance?
  - What's the quality of the underlying data?

- Is there a way to incorporate insertion direction into graph?

- Add in functional enrichment analysis?

- Create a feature table for NetCIS and comparing to the other tools (5 tools total, 6 with NetCIS)
  - Rylee is assigned this
