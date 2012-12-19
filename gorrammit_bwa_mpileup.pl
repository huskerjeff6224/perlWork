#!/usr/bin/perl
use warnings;
use strict;

usage() unless $#ARGV == 0;
my $in     = shift;
my $config = `cat $in`;
usage() if $config !~ /READS_DIR/;

######## Start Variables ########

my $bin       = ( $config =~ /BIN\s+(\S+)/ )       ? $1 : "/opt/var_calling";
my $snpEff    = ( $config =~ /SNPEFF\s+(\S+)/ )    ? $1 : "/opt/snpeff";
my $ref_dir   = ( $config =~ /REF_DIR\s+(\S+)/ )   ? $1 : "/data/genomes/Homo_sapiens/UCSC/hg18";
my $bt2_idx   = ( $config =~ /BT2\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/BowtieIndex/ucsc.hg19";
my $bwa_ref   = ( $config =~ /BWA\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19";
my $threads   = ( $config =~ /THREADS\s+(\S+)/ )   ? $1 : "24";
my $memory    = ( $config =~ /MEMORY\s+(\S+)/ )    ? $1 : "48";
my $reads_dir = ( $config =~ /READS_DIR\s+(\S+)/ ) ? $1 : ".";
my $step      = ( $config =~ /STEP\s+(\S+)/ )      ? $1 : "0";
my $gatk      = "$bin/GenomeAnalysisTK.jar";
my $dbsnp     = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $omni      = "$ref_dir/Annotation/Variation/omni.vcf";
my $hapmap    = "$ref_dir/Annotation/Variation/hapmap.vcf";
my $mills     = "$ref_dir/Annotation/Variation/indels.vcf";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
die "There are no read 1 fastq reads in $reads_dir. The read 1 reads must be formatted as follows: *_R1.fastq.\n" unless ( `ls $reads_dir/*_R1_001.fastq` );
die "There are no read 2 fastq reads in $reads_dir. The read 2 reads must be formatted as follows: *_R2.fastq.\n" unless ( `ls $reads_dir/*_R2_001.fastq` );
chomp ( my @reads  = `ls $reads_dir/*fastq` );
chomp ( my $time   = `date +%T` );

print "Options used     :\n",
      "\tBIN      : $bin\n",
      "\tSNPEFF   : $snpEff\n",
      "\tREF_DIR  : $ref_dir\n",
      "\tBT2_IDX  : $bt2_idx\n",
      "\tTHREADS  : $threads\n",
      "\tMEMORY   : $memory\n",
      "\tREADS_DIR: $reads_dir\n",
      "\tSTEP     : $step\n";

######## End Variables ########

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name) = $reads[$i] =~ /^.+\/(.+?)_/;
    my $R1     = $reads[$i];
    my $R2     = $reads[$i+1];
    my $JAVA_pre = "java -Xmx${memory}g -jar";
    my $GATK_pre   = "$JAVA_pre $gatk -T";
    my @steps      = (
"bwa aln -q 1 -t $threads $bwa_ref $R1 > R1.sai",
                      "bwa aln -q 1 -t $threads $bwa_ref $R2 > R2.sai",
                      "bwa sampe $bwa_ref R1.sai R2.sai $R1 $R2 > $name.sam",
"samtools view -bS $name.sam -o $name.bam",
                       "$JAVA_pre $bin/ReorderSam.jar INPUT=$name.bam OUTPUT=$name.sorted.bam REFERENCE=$ref VALIDATION_STRINGENCY=LENIENT",
                       "$JAVA_pre $bin/AddOrReplaceReadGroups.jar I=$name.sorted.bam O=$name.Not_Ref_Sorted_fixed_RG.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                       "$JAVA_pre $bin/ReorderSam.jar INPUT=$name.Not_Ref_Sorted_fixed_RG.bam OUTPUT=$name.fixed_RG.bam REFERENCE=$ref VALIDATION_STRINGENCY=LENIENT",
                       "samtools index $name.fixed_RG.bam",
                       "$GATK_pre RealignerTargetCreator -R $ref -I $name.fixed_RG.bam -known $dbsnp -o $name.indel_realigner.intervals",
                       "$GATK_pre IndelRealigner -R $ref -I $name.fixed_RG.bam -known $dbsnp -o $name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $name.indel_realigner.intervals",
                       "$JAVA_pre $bin/MarkDuplicates.jar INPUT=$name.indels_realigned.bam OUTPUT=$name.realigned_2.bam REMOVE_DUPLICATES=True METRICS_FILE=$name.txt VALIDATION_STRINGENCY=LENIENT",
                       "$JAVA_pre $bin/ReorderSam.jar INPUT=$name.realigned_2.bam OUTPUT=$name.recalibrated.bam REFERENCE=$ref VALIDATION_STRINGENCY=LENIENT",
                       "samtools index $name.recalibrated.bam",
                       "samtools mpileup  -uBf  $ref  $name.recalibrated.bam  > $name.mpileup",
                       "bcftools view -bvcg $name.mpileup > $name.bcf",
                       "bcftools view $name.bcf > $name.mpileup.raw.vcf",
                       "perl /safer/tools/annovar/convert2annovar.pl -format vcf4 $name.mpileup.raw.vcf > $name.mpileup.raw.annovar.input.txt",
                       "perl /safer/tools/annovar/summarize_annovar.pl $name.mpileup.raw.annovar.input.txt -buildver hg19 -outfile $name.mpileup.raw -verdbsnp 135 -ver1000g 1000g2012apr -veresp 6500 /safer/tools/annovar/humandb/"
                       );

    chomp ( $time = `date +%T` );
    print "[$time][- / -] Working on sample $name.\n";
    my $nom;
    for ( my $i = $step; $i < @steps; $i++ )
    {
        $nom = sprintf ( "%02d", $i );
        my $current_step = $steps[$i];
        chomp ( $time = `date +%T` );
        my ($clean_step) = $current_step;
        $clean_step =~ s/ -/\n                  -/g if length ($clean_step) > 256;
        print "[$time][$nom/$#steps] Running this step: \n\n", " "x18, "$clean_step\n\n";
        system ( $current_step );
    }
   #system ( "mkdir $name; mv $name.* $name;" ); 
}

sub usage
{
    die <<USAGE;

    Usage: perl $0 <configuration_file.txt>

    Your configuration file MUST have a READS_DIR specified.

    Configuration options available

      OPTION    Default                                  Description
      BIN       /opt/var_calling                         Absolute location of the Picard Tools and GATK jar files
      SNPEFF    /opt/snpeff                              Absolute location of snpEff and its requisite files
      REF_DIR   /data/genomes/Homo_sapiens/UCSC/hg19/    Absolute location of the reference directory
      BT2       REF_DIR/Sequence/BowtieIndex/ucsc.hg19   Absolute location of the Bowtie2 index
      THREADS   24                                       Number of threads to use in parallelizable modules
      MEMORY    48                                       Amount of memory, in gigabytes, to use
      READS_DIR N/A                                      Absolute location of the reads that are going to be used
      STEP      0                                        The step to start at in the pipeline (0-indexed).     

USAGE
}
