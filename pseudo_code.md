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

## Construct network (using pseudo code from graph framework article)

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

- Is there a way to incorporate insertion direction into graph?

- Add in functional enrichment analysis?
