# SARS-CoV-2 Nanopore sequencing data analysis:

# Activate conda environment for required packages:
conda activate nanopore_covid
wait

# Basecalling and Barcoding:
mkdir -p basecall_files/
guppy_basecaller -i fast5/ -s basecall_files/ -c dna_r9.4.1_450bps_hac.cfg -x "cuda:0 cuda:1"

wait

mkdir -p demultiplexed/
guppy_barcoder -i basecall_files/pass -s demultiplexed/ --barcode_kits SQK-RBK110-96 -x "cuda:0 cuda:1"

wait

# Merge fastq files:
mkdir -p fastq_files/
cd demultiplexed/
for i in barcode*; do cd $i; cat *fastq > ../../fastq_files/${i}.fastq; cd ..; done
wait
cd ..

# create list of files for input: 
ls fastq_files/*.fastq | sed 's/fastq_files\/\(.*\).fastq/\1/' > files.list

# Adapter trimming:
mkdir -p trimmed_fastq/
parallel --verbose -j 10 'porechop -i fastq_files/{}.fastq -o trimmed_fastq/{}_trimmed.fastq -t 16 --discard_middle' ::: $(cat files.list)
wait

# Fetching quality control files:
mkdir -p nanoplot_results/
parallel --verbose -j 10 'NanoPlot --fastq fastq_files/{}.fastq -t 16 -o nanoplot_results/{}_raw' ::: $(cat files.list)
wait
parallel --verbose -j 10 'NanoPlot --fastq trimmed_fastq/{}_trimmed.fastq -t 16 -o nanoplot_results/{}_trimmed' ::: $(cat files.list)
wait

# Mapping with SARS-CoV-2 reference genome:
mkdir -p alignment/
parallel --verbose -j 10 'minimap2 -a -x map-ont -t 16 NC_045512.fasta trimmed_fastq/{}_trimmed.fastq | samtools view -bS -F 4 - | samtools sort -o alignment/{}_sorted.bam' ::: $(cat files.list)
wait

# Indexing BAM files:
parallel --verbose -j 10 'samtools index alignment/{}_sorted.bam' ::: $(cat files.list)
wait

# Variant Calling:
mkdir -p results/bcftools_mpileup/
parallel --verbose -j 10 'bcftools mpileup --count-orphans --no-BAQ --max-depth 1000000 --max-idepth 100000 -m 20 --fasta-ref NC_045512.fasta --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR alignment/covid/{}_sorted.bam > results/bcftools_mpileup/{}.mpileup' ::: $(cat files.list)
wait

mkdir -p results/variant_analysis/filter_vcfs/
parallel --verbose -j 10 'bcftools call -c --ploidy 1 --keep-alts --keep-masked-ref --variants-only results/bcftools_mpileup/{}.mpileup -o results/variant_analysis/{}.vcf.gz --output-type z' ::: $(cat files.list)
wait

parallel --verbose -j 10 '/data2/Safal/COVID/tools/bcftools-1.13/bcftools index -t results/variant_analysis/{}.vcf.gz' ::: $(cat files.list)
wait

parallel --verbose -j 10 'bcftools view results/variant_analysis/{}.vcf.gz --output-file results/variant_analysis/filter_vcfs/{}_filter.vcf.gz --output-type z --include "INFO/DP>=5 && INFO/AF1 > 0.5"' ::: $(cat files.list)
wait

parallel --verbose -j 10 '/data2/Safal/COVID/tools/bcftools-1.13/bcftools index -t results/variant_analysis/filter_vcfs/{}_filter.vcf.gz' ::: $(cat files.list)
wait

# Generating consensus sequences:
mkdir -p results/new_consensus/
parallel --verbose -j 10 'bedtools genomecov -bga -split -ibam alignment/{}_sorted.bam -g NC_045512.fasta | awk "\$4 < 10" | bedtools merge > results/new_consensus/{}.mask.bed' ::: $(cat ../../files.list)
wait

parallel --verbose -j 10 'bcftools consensus results/variant_analysis/filter_vcfs/{}_filter.vcf.gz -f NC_045512.fasta -m results/new_consensus/{}.mask.bed -o results/new_consensus/{}.masked.consensus.fa' ::: $(cat ../../files.list)
wait

cd results/new_consensus/
for i in *.masked.consensus.fa; do ID_NAME=$(basename "$i" .masked.consensus.fa); sed -i 's/>\(.*\)/>'$ID_NAME'/' $i; done
wait
cd ../..


# Annotation:
mkdir -p results/variant_annotation/
bcftools merge --merge all results/variant_analysis/filter_vcfs/*.vcf.gz --missing-to-ref --output results/variant_annotation/Merged.vcf.gz --output-type z
wait

cd results/variant_annotation/

bcftools norm -m-any Merged.vcf.gz | bcftools view --include 'GT!="0"' -Oz -o Merged_norm.vcf.gz  #to split multiallelic variants
wait

snpEff NC_045512.2 Merged_norm.vcf.gz -canon -no-upstream -no-downstream -classic -csvStats variant_sites.snpEff.csv |bgzip -c > Merged.snpEff.vcf.gz
wait

# Extract and modify annotation files:
SnpSift extractFields Merged.snpEff.vcf.gz -s "," -e "." CHROM POS REF ALT ANN[*].GENE EFF[*].FEATUREID EFF[*].EFFECT EFF[*].IMPACT EFF[*].AA ANN[*].BIOTYPE > annotation.tsv
wait

vcftools --gzvcf Merged.snpEff.vcf.gz --extract-FORMAT-info GT --out genotypes
wait

paste annotation.tsv genotypes.GT.FORMAT | awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' > complete_annotation.tsv


