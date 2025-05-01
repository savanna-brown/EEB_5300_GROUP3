# **Comparative analysis of odorant receptor (OR) gene expansion in eusocial versus non-eusocial Hymenopteran and non-Hymenopteran insects**
## *Savanna Brown, Weijun Liang, Tyler Elias*<br><br>

# **Introduction**

In animals, eusociality is a complex and extreme social structure defined by having three characteristics: reproductive division of labor (castes), cooperative care of young, and overlapping generations (Wilson and Hölldobler, 2005). In these societies, some individuals forgo personal reproduction to help rear the offspring of genetic relatives, a scenario often explained by inclusive fitness or kin selection theory (Eberhard, 1975). Inclusive fitness posits that genes promoting behavior that increases the number of offspring of close genetic kins can spread (Eberhard, 1975). Kin selection, a kind of selection that favors a trait due to its positive effects on the reproductive success of kins, is therefore regarded as a primary mechanism enabling the evolution of eusociality (Eberhard, 1975).

Eusocial insects mainly occur in the orders Hymenoptera and Blattodea (Hölldobler and Wilson, 1990; Ross and Matthew, 2018; Thorne, 1997). Within Hymenoptera (bees, wasps, ants), eusociality has evolved multiple times: many bees and wasps are eusocial, and all ants are eusocial (Hölldobler and Wilson, 1990; Ross and Matthew, 2018). In Blattodea, termites are eusocial cockroaches and live in large colonies with a queen, king, sterile workers, and soldiers (Thorne, 1997). By contrast, typical cockroaches lack the three characteristics of eusociality (Bell et al., 2007). Therefore, eusocial Blattodea (termites) form true colonies with castes, whereas non-social Blattodea (cockroaches) do not.

Insect eusociality depends critically on chemical communication. Social insects use pheromones and other odorant cues to recognize nestmates, mark trails, regulate caste behavior, and coordinate many aspects of colony life (Zhou et al., 2015). In the sensilla of antennae, odorant receptors (ORs) are membrane proteins in the olfactory sensory neurons that detect volatile compounds (Steinbrecht, 2007; Zhou et al., 2015). Each OR is tuned to specific chemicals, and OR signals are processed in the brain’s antennal lobes, and notably, ORs are essential for sensing chemical cues from nestmates and castes (Steinbrecht, 2007; Zhou et al., 2015). Because chemosensation is so central to eusocial behavior, sophisticated OR repertoires have been hypothesized to underpin social life (Zhou et al., 2015).

Building on this idea, genomic studies found that some eusocial insects possess very large OR gene families. For example, sequenced ant and honeybee genomes revealed some of the largest OR repertoires known in insects. Zhou et al. (2015) reported widespread chemoreceptor gene expansions in ants and bees and suggested that these expansions likely facilitated the transition to eusociality. In particular, certain OR subfamilies (such as the 9-exon subfamily) are greatly expanded in ants and social bees and have been proposed to mediate kin recognition and pheromone communication in eusocial Hymenoptera (Zhou et al., 2015). This led to the expectation that eusocial species might in general maintain more OR genes than non-social relatives.

However, recent broad comparisons challenge the simple link between eusociality and OR gene count. Gautam et al. (2024) found no evidence to support that OR repertoires in Hymenoptera are linked to the evolution of eusociality. Instead, being wingless might shape the expansion of OR genes. Thus, the relationship between OR gene family size and eusocial lifestyle remains an open question. Furthermore, few studies have tackled the OR gene count in eusocial Blattodea, and it is unclear whether eusociality in this order is also associated with relatively large OR gene repertoires.<br><br>

# Research Question, Hypothesis, and Objective

**Research question: Is there a correlation between the number of odorant receptor (OR) genes and eusociality in Hymenoptera and Blattodea?**

To answer our research question, we compare OR gene counts across Hymenoptera and Blattodea. 

We hypothesize that eusocial taxa will exhibit larger OR repertoires across 19 OR subfamilies than non-eusocial taxa. <br><br>

## Study System
To explore the relationship between odorant receptor (OR) gene repertoires and eusociality, we selected representative species from four focal categories:
- **Eusocial Hymenoptera**
- **Non-eusocial Hymenoptera**
- **Eusocial Blattodea (termites, formerly order Isoptera)**
- **Non-eusocial Blattodea (cockroaches)**

An outgroup species from Odonata was also included for comparative purposes. Species were chosen based on the availability of high-quality annotated proteomes.<br>

| **Category**               | **Species**                                                                 |
|---------------------------|------------------------------------------------------------------------------|
| Eusocial Hymenoptera      | *Apis mellifera* (Western honeybee)<br>*Harpegnathos saltator* (jumping ant)<br>*Vespa crabro* (European hornet) |
| Non-eusocial Hymenoptera  | *Nasonia vitripennis* (jewel wasp)<br>*Microplitis demolitor* (braconid parasitoid wasp)<br>*Orussus abietinus* (parasitic wood wasp) |
| Eusocial Blattodea        | *Zootermopsis nevadensis* (dampwood termite)<br>*Cryptotermes secundus* (drywood termite)<br>*Coptotermes formosanus* (Formosan subterranean termite) |
| Non-eusocial Blattodea    | *Periplaneta americana* (American cockroach)<br>*Diploptera punctata* (Pacific beetle cockroach) |
| Outgroup                  | *Ischnura elegans* (blue-tailed damselfly) <br>                                 



All eusocial species have clear reproductive and non-reproductive caste systems, whereas their listed non-social relatives lack such division of labor. By surveying OR gene family sizes in each of these taxa, **we aim to determine whether a consistent expansion of OR genes accompanies the evolution of eusociality in both orders.** <br><br>


# **Motivation**

We started this project because we are interested in insects and gene family evolution. While exploring project ideas, we found several studies on the expansion of odorant receptors in eusocial insects. Eusociality is a fascinating topic because it is a unique way of living that influences the evolution of eusocial species.

Eusocial insects are found in two major groups: Hymenoptera (which includes all ants, some bees, and some wasps) and Blattodea (which includes mostly termites, while non-eusocial species are mostly cockroaches). Previous research suggests a link between eusociality and the expansion of odorant receptors, but most studies focus on Hymenoptera. There has been less research on eusociality in Blattodea. To expand on existing studies, we took a different comparative approach by including a non-Hymenopteran lineage. We specifically focused on the 9-exon family, a well-documented group of odorant receptor genes in Hymenoptera.


<br><br><br><br><br><br><br><br>
#INSERT REFERENCES HERE# 
<br><br><br><br><br><br>



# Workflow <br>

## **Data Collection**

We searched NCBI for each target species to ensure that a predicted proteome in FASTA format was available for download. Each protein fasta was downloaded directly from NCBI  using a wget command and the FTP links. <br>


```bash
# wget download
wget –c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/059/185/GCA_964059185.1_ibEctPall1.hap1.1/GCA_964059185.1_ibEctPall1.hap1.1_genomic.fna.gz
```
<br><br>

## **Quality Control**

To assess completeness of selected proteomes, we ran **BUSCO** (Benchmarking Universal Single-Copy Orthologs) on each protein FASTA file using the `insecta_odb10` database.

```
module load busco

# Define directories and database
GENOME_DIR="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/pep_assemblies"
BUSCO_DIR="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/BUSCO"
LINEAGE="/isg/shared/databases/BUSCO/odb10/lineages/insecta_odb10"

# Create BUSCO results directory if it doesn't exist
mkdir -p ${BUSCO_DIR}

# Change to the BUSCO results directory
cd ${BUSCO_DIR}

# Loop through each genome file that ends with ".fna" and run BUSCO
for genome_file in ${GENOME_DIR}/*.faa; do
  base_name=$(basename ${genome_file})
  output_name=$(echo ${base_name} | cut -d'_' -f1,2)

  # Run BUSCO
  echo "Running BUSCO on ${genome_file}"
  busco -i ${genome_file} -o ${output_name}_busco -l ${LINEAGE} -m protein -c 8 -f

  echo "BUSCO analysis complete for ${genome_file}. Results saved to ${BUSCO_DIR}/${output_name}_busco"
done

```
<br>


 All selected proteomes had a BUSCO completeness score greater than 85% and were considered sufficient for downstream analysis. <br><br>
 

### BUSCO Summary

![BUSCO Summary](figures/BUSCO_summary.png)


<br><br>
### OR Gene Screening

We also used **DIAMOND** (v2.1.8; Buchfink et al., 2021) to perform sequence similarity searches of OR genes in our selected proteomes. A DIAMOND database was built from a curated set of *Polistes* odorant receptor (OR) proteins from Zhou et al. (2015). Each target species’ protein FASTA file was aligned against this database using the `blastp` mode with default parameters to identify putative OR genes. <br><br>


### Diamond databasing of curated Polistes OR protein set

```bash
module load diamond/2.1.8

diamond makedb --in polistes_OR_protein_sequences.fasta \
               --db Polistes_OR_Protein_Sequences.dmnd \
               --threads 16
```
<br>
### Similarity search - screen for OR genes
```bash
module load diamond/2.1.8

# Define paths
DB="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/Polistes_OR_Protein_Sequences.dmnd"  
INPUT_DIR="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/pep_assemblies"  
OUTPUT_DIR="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_SCREEN_RESULTS"

# Loop through all .faa files in the input directory and run DIAMOND BLASTP
for QUERY in $INPUT_DIR/*.faa; do
    # Extract the base name of the input file (e.g., file1.faa -> file1)
    BASENAME=$(basename $QUERY .faa)
    
    # Define the output file for this query
    OUT="$OUTPUT_DIR/${BASENAME}_vs_pep_database.tsv"
    
    # Run DIAMOND BLASTP
    echo "Running DIAMOND BLASTP for $QUERY..."
    diamond blastp \
        --db $DB \
        --query $QUERY \
        --out $OUT \
        --outfmt 6 qseqid sseqid pident length evalue bitscore qseq sseq \
        --threads 16
    echo "DIAMOND BLASTP completed for $QUERY. Output saved to $OUT."
done
```
<br><br>

All selected proteomes showed successful OR gene retrieval with this method, indicating that our chosen input data is sufficient for downstream analysis. 
<br><br>


## **OrthoFinder**

OrthoFinder is a comparative genomics tool that identifies orthogroups by combining sequence similarity searches, gene tree inference, and species tree construction (Emms & Kelly, 2017). It begins with an all-vs-all sequence search across protein FASTA files from all species, clusters homologous genes into orthogroups, and then reconstructs gene trees for each orthogroup. It also infers a rooted species tree and maps gene duplication events onto the tree by comparing the topology of gene trees and the species tree, ultimately giving key information about how genes are related both among and within species.

For this project, we ran OrthoFinder (v2.5.4) using predicted protein sequences from our selected hymenopteran and non-hymenopteran species to define orthogroups across lineages. <br>


```bash
module load OrthoFinder/2.5.4
module load mafft/7.471 FastME/2.1.5 RAxML/8.2.11 fasttree/2.1.10

# input directory
ORTHO_INPUT="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/ORTHOFINDER_INPUT"

# path to OF python script on cluster
ORTHOFINDER="/isg/shared/apps/OrthoFinder/2.5.4/OrthoFinder_source/orthofinder.py"

python $ORTHOFINDER -f $ORTHO_INPUT/ \
-t 16 \
-a 16 \
-M msa \
-A mafft \
-T fasttree \
-S diamond \
-n OrthoFinder_run1
```
<br><br>



## Diamond Similarity Search

These orthogroups served as the basis for identifying candidate odorant receptor (OR) gene families. To identify putative OR gene families, we used our DIAMOND database of curated *Polistes* OR proteins to perform a similarity search against the longest representative protein sequence of each orthogroup.
<br>

### Longest Representative Protein Pull

```bash
OUTPUT_FILE="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/longest_per_orthogroup.fasta"
INPUT_DIR="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/ORTHOFINDER_INPUT/OrthoFinder/Results_OrthoFinder_run1/Orthogroup_Sequences"

# loop through each orthogroup file
for FILE in "$INPUT_DIR"/*.fa; do
    BASENAME=$(basename "$FILE" .fa)  # Extract the orthogroup ID

    # extract the longest sequence
    awk -v og="$BASENAME" '
    /^>/ { 
        if (seq) { 
            if (length(seq) > max_len) { 
                max_len = length(seq); 
                longest_header = header; 
                longest_seq = seq; 
            } 
        } 
        header = $0; seq = ""; next 
    } 
    { seq = seq $0 } 
    END { 
        if (length(seq) > max_len) { 
            longest_header = header; 
            longest_seq = seq; 
        } 
        print ">" og "_" substr(longest_header, 2) "\n" longest_seq 
    }' "$FILE" >> $OUTPUT_FILE
done
```
<br>


### Diamond blastp

```bash
module load diamond/2.1.8 

DB="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_screen/Polistes_OR_Protein_Sequences.dmnd"  # db
INPUT="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/longest_per_orthogroup.fasta"  # input
OUT="/home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits.tsv"  # output

diamond blastp \
    --db $DB \
    --query $INPUT \
    --out $OUT \
    --outfmt 6 qseqid sseqid pident length evalue bitscore qseq sseq \
    --threads 8
```
<br><br>

Alignment results were filtered based on alignment quality thresholds of e-value less than 1e-10 and an alignment percent identity greater than or equal to 30%. These thresholds were chosen after exploring results and finding an ideal balance of data preservation while eliminating noise from lower quality hits that made 1:1 orthogroup:OR gene family calls unclear.
 
```bash
# get a clean summary of ther diamond results to see what families hit to what OGs
awk -F'\t' '{split($2, a, "_"); print $1 "\t" a[2]}' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits.tsv > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_all_summary.tsv


# filter the DIAMOND output by e-value (e-value < 1e-10) and coverage ( >= 30% percent ID)
awk '$3 >= 30 && $5 < 1e-10' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits.tsv > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_filtered.tsv


# get summary of filtered results to see what families hit to what OGs
awk -F'\t' '{split($2, a, "_"); print $1 "\t" a[2]}' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_filtered.tsv > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_filtered_summary.tsv
```

<br><br>
### Making Gene Family Calls

Even with the alignment thresholds applied, some orthogroups matched to multiple OR gene families, complicating the ability to assign a clean 1:1 match between orthogroup and gene family.

To address this, we implemented a classification scheme based on the following rules:

- **Single-hit orthogroups**: Orthogroups that matched only one OR gene family were classified as that family.
- **One-outlier orthogroups**: Orthogroups where all but one match belonged to the same gene family were classified as the majority (non-outlier) family.
- **Multi-family or ambiguous orthogroups**: Orthogroups with matches to more than one outlier or more than two gene families were labeled as **"Unclassified"**.

An awk script was used to execute this classification using if/else logic:
<br>

```bash
awk -F'\t' '{
    counts[$1][$2]++;  # count occurrences of each family (col 2) for each orthogroup (col 1)
    total[$1]++;       # track total hits for each orthogroup
} END {
    for (og in counts) {
        max_count = 0;
        majority_family = "";
        unique_families = 0;

        # determine the majority family and count unique families
        for (family in counts[og]) {
            unique_families++;
            if (counts[og][family] > max_count) {
                max_count = counts[og][family];
                majority_family = family;
            }
        }

        # decide classification based on the rules
        if (unique_families == 1) {
            print og "\t" majority_family;  # all hits same
        } else if (max_count == total[og] - 1) {
            print og "\t" majority_family;  # one outlier
        } else {
            print og "\tUnclassified";     # unclassified
        }
    }
}' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_filtered_summary.tsv > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_classified.tsv
```

<br><br>

## **Subsetting into final DF**

Finally, these orthogroup to OR gene family identities were merged with orthogroup gene copy number data from Orthofinder to create a final data frame with gene copy number for each OR gene family for each species.
<br>


```bash
 subet the OG output gene counts file to only the OGs that are in the OR gene family summary
awk 'NR==FNR {  # Process the classified file first
   split($1, a, "_");  # Extract the orthogroup part before the first underscore
   orthogroups[a[1]];  # Store the orthogroup in an array
   next;
}
FNR == 1 || $1 in orthogroups {  # Print the header (FNR == 1) or matching orthogroups
   print;
}' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_classified.tsv \
  /home/FCAM/eeb5300/usr3/GROUP_PROJECT/ORTHOFINDER_INPUT/OrthoFinder/Results_OrthoFinder_run1/Orthogroups/Orthogroups.GeneCount.tsv \
  > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/Orthogroups.GeneCount_subset.tsv

  # add in family names to this count summary
awk 'NR==FNR {  # Process the classified file first
   split($1, a, "_");  # Extract the orthogroup name only, before the first underscore
   families[a[1]] = $2;  # Store the family name in an array
   next;
}
FNR == 1 {  # For the gene count file, print the header with an added "Family" column
   print "Orthogroup\tFamily\t" substr($0, index($0, "\t") + 1);
   next;
}
$1 in families {  # For matching orthogroups, add the family name as the second column
   print $1 "\t" families[$1] "\t" substr($0, index($0, "\t") + 1);
}' /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/OG_hits_classified.tsv \
  /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/Orthogroups.GeneCount_subset.tsv \
  > /home/FCAM/eeb5300/usr3/GROUP_PROJECT/DIAMOND_OGs/Orthogroups.GeneCount_with_family.tsv
```
<br><br>
