imputehp
========

A protocol to impute haptoglobin haplotypes from surrounding genotypes computed from genotype array or whole genome sequence data using the reference panel from:

```
Boettger L., McCarroll S., et al. Recurring exon deletions in the HP (haptoglobin)
gene contribute to lower blood cholesterol levels. Nat. Genet. 48, 359–366 (2016)
```

This reference panel was generated as explained <a href="http://mccarrolllab.org/wp-content/uploads/2014/12/ng.3510.pdf">here</a>. For any feedback, send an email to giulio.genovese@gmail.com or mccarroll@genetics.med.harvard.edu

![](http://mccarrolllab.org/wp-content/uploads/2018/09/Slideshow-08-1170x500.jpg)

Installation
============

Install basic tools (Debian/Ubuntu specific):

```
sudo apt install wget gzip samtools bcftools plink1.9 openjdk-11-jre-headless
```

Preparation steps
```
mkdir -p $HOME/res
```

Download Beagle binary
```
wget -P $HOME/res/ https://faculty.washington.edu/browning/beagle/beagle.25Nov19.28d.jar
```

Download reference panels
```
wget -P $HOME/res/ https://personal.broadinstitute.org/giulio/panels/HP_Euroref_1kgOMNI_HM3_merged.GRCh3{7,8}.vcf.gz
```

Run haptoglobin imputation from an input VCF
============================================

Run imputation using Beagle
```
vcf="..."
out="..."
build=38 # build=37
declare -A reg=( ["37"]="16:71070878-73097663" ["38"]="chr16:71036975-73063764" )

bcftools view --no-version "$vcf" -r ${reg[$build]} | \
  java -Xmx8g -jar $HOME/res/beagle.25Nov19.28d.jar gt=/dev/stdin \
  ref=$HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh$build.vcf.gz out="$out" \
  map=<(bcftools query -f "%CHROM\t%POS\n" $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh$build.vcf.gz | \
  awk '{print $1"\t.\t"$2/1e7"\t"$2}')
```

Extract imputed haptoglobin alleles into a table
```
out="..."
build=38 # build=37
declare -A reg=( ["37"]="16:72092044-72092044" ["38"]="chr16:72058145-72058145" )

bcftools index -ft "$out.vcf.gz" && \
bcftools query -f "[%SAMPLE\t%ALT\t%GT\n]" "$out.vcf.gz" -r ${reg[$build]} | tr -d '[<>]' | \
  awk -F"\t" -v OFS="\t" '{split($2,a,","); a["0"]="NA"; split($3,b,"|"); \
  print $1,a[b[1]],a[b[2]]}' > "$out.tsv"
```

Build the reference panels yourself
===================================

This section is only in case you want to build the reference panels yourself, you can skip it otherwise

Download GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/res/human_g1k_v37.fasta
samtools faidx $HOME/res/human_g1k_v37.fasta
```

Download GRCh38 human genome reference
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download haptoglobin reference panel in Beagle 3 format with marker positions for the GRCh37 human genome reference
```
wget -P $HOME/res/ https://raw.githubusercontent.com/freeseek/imputehp/master/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased
```

Liftover marker positions for the GRCh38 human genome reference
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
tail -n+2 $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased | grep -v HP | \
    awk '{print "16\t"$2}' > $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh37
tail -n+2 $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased | grep -v HP | \
  awk '{print "chr16\t"$2-1"\t"$2}' | \
  ./liftOver /dev/stdin hg19ToHg38.over.chain.gz /dev/stdout /dev/stderr | \
  awk '{print $1"\t"$3}' > $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh38
```

Generate haptoglobin reference panels in VCF format for both the GRCh37 and GRCh38 human genome references
```
declare -A fasta=( ["37"]="$HOME/res/human_g1k_v37.fasta" ["38"]="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" )
declare -A chr=( ["37"]="16" ["38"]="chr16" )
declare -A pos=( ["37"]="72092044" ["38"]="72058145" )

grep HP $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased | \
  awk 'BEGIN {x["AAAAABAAA"]=1; x["AABBBAAAB"]=4; x["BAABBAABA"]=3; x["BBAAAABAA"]=2}
  {for (i=3; i<=NF; i++) y[i]=y[i]$i}
  END {printf "<HP1A>,<HP1B>,<HP2A>,<HP2B>\t.\t.\t.\tGT";
  for (i in y) printf "\t"x[y[i]]}' | \
  sed 's/\t\([1-9]\)\t\([1-9]\)/\t\1|\2/g' > tmp

for build in 37 38; do
  tail -n+2 $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased | grep -v HP | tr ' ' '\t' | cut -f3- | \
    paste $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh$build - | \
    sed 's/\t\([ACGT]\)\t\([ACGT]\)/\t\1\2/g' | \
    bcftools convert --no-version --tsv2vcf /dev/stdin -c CHROM,POS,AA -f ${fasta[$build]} --samples \
    $(head -n1 $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.bgl.phased | tr ' ' '\n' | \
    tail -n+3 | uniq | tr '\n' ',' | sed 's/,$//') | tr '/' '|' | \
    awk -v chr=${chr[$build]} -v pos=${pos[$build]} 'NR==FNR {line=$0}
    NR>FNR {if (line && $0!~"^#" && $2>pos) {print chr"\t"pos"\t.\tG\t"line; line=""} print}' tmp - | \
    bcftools view --no-version -Oz \
    -o $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh$build.vcf.gz
done
/bin/rm tmp
```

Check for consistency of the haptoglobin reference panel
========================================================

Convert the haptoglobin and 1000 Genomes project reference panels to plink and then merge to compute consistency
```
url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
bcftools view --no-version $url -r 16:71070878-73097663 | \
  awk 'NF==2 {print "##contig=<ID=16,length=90354753>"} {print}' | \
  bcftools view --no-version -v snps | \
  bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
  $HOME/bin/plink --vcf /dev/stdin --keep-allele-order --const-fid --make-bed \
  --out ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes

bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' \
  $HOME/res/HP_Euroref_1kgOMNI_HM3_merged.GRCh37.vcf.gz | \
  $HOME/bin/plink --vcf /dev/stdin --biallelic-only --keep-allele-order --const-fid --make-bed \
  --out HP_Euroref_1kgOMNI_HM3_merged.GRCh37

plink --bfile HP_Euroref_1kgOMNI_HM3_merged.GRCh37 \
  --bmerge ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes --merge-mode 6
```

You should get the following result
```
241696 overlapping calls, 241696 nonmissing in both filesets.
240611 concordant, for a concordance rate of 0.995511.
```
