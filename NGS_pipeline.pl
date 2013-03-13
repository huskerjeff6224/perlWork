#!/usr/bin/perl
use warnings;
use strict;
use v5.14;
use lib "/miseqdata/tools/vcftools_0.1.10/perl/";
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
my $targetInterval = "$ref_dir/Annotation/Variation/";
my $dbsnp     = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $omni      = "$ref_dir/Annotation/Variation/omni.vcf";
my $hapmap    = "$ref_dir/Annotation/Variation/hapmap.vcf";
my $mills     = "$ref_dir/Annotation/Variation/indels.vcf";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
my $Gene_Panel = "/miseqdata/Annotation_files/52GenesWithExonCoordinatesPLUS10.bed";
my @recalibrated_files;
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
print join "\n",@reads;
for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name) = $reads[$i] =~ /^.+\/(.+?_.+?_.+?)_/;
    my $R1     = $reads[$i];
    my $R2     = $reads[$i+1];
    
    my $JAVA_pre = "java -Xmx${memory}g -jar";
    my $GATK_pre   = "$JAVA_pre $gatk -T";
    my @steps      = (
                       "bwa aln -q 1 -t $threads $bwa_ref $R1 > R1.sai",
                       "bwa aln -q 1 -t $threads $bwa_ref $R2 > R2.sai",
                       "bwa sampe -r \'\@RG\tID:Default\tLB:Library\tPL:Illumina\tSM:SAMPLE\' $bwa_ref R1.sai R2.sai $R1 $R2 > $name.sam",
                       "samtools view -bS $name.sam -o $name.bam",
                       "$JAVA_pre $bin/picard/SortSam.jar INPUT=$name.bam OUTPUT=$name.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT",
                       "samtools index $name.sorted.bam",
                       "$JAVA_pre $bin/picard/MarkDuplicates.jar INPUT=$name.sorted.bam OUTPUT=$name.sorted.nodups.bam REMOVE_DUPLICATES=True METRICS_FILE=$name.txt VALIDATION_STRINGENCY=LENIENT",
                       "samtools index $name.sorted.nodups.bam",
                       "$GATK_pre RealignerTargetCreator -R $ref -I $name.sorted.nodups.bam -known $dbsnp -o $name.indel_realigner.intervals",
                       "$GATK_pre IndelRealigner -R $ref -I $name.sorted.nodups.bam -known $dbsnp -o $name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $name.indel_realigner.intervals",
                       "samtools index $name.indels_realigned.bam",
                       "$GATK_pre -T BaseRecalibrator -R $ref -nct 3 -I $name.indels_realigned.bam -knownSites $dbsnp -o $name.recalibrated.grp",
                       "$GATK_pre -T PrintReads -R $ref -BQSR $name.recalibrated.grp -I $name.indels_realigned.bam  -o $name.recalibrated.bam",
                       "samtools index $name.recalibrated.bam",
                       "$GATK_pre UnifiedGenotyper -I $name.recalibrated.bam -nt 5 -R $ref -D $dbsnp -glm SNP -o $name.raw.snvs.GATK.vcf",
                       "$GATK_pre UnifiedGenotyper -I $name.recalibrated.bam -nt 5 -R $ref -D $mills -glm INDEL -o $name.raw.indel.GATK.vcf",
                       "$GATK_pre VariantFiltration -R $ref --variant $name.raw.snvs.GATK.vcf -o $name.filtered_variants.snvs.GATK.vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP < 15' --filterName 'GATK_BP_filters' --maskExtension 10 --mask $Gene_Panel --maskName PANEL",
                       "$GATK_pre VariantFiltration -R $ref  --variant $name.raw.indel.GATK.vcf -o $name.filtered_variants.indel.GATK.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8' --filterName 'GATK_BP_filters' --maskExtension 10 --mask $Gene_Panel --maskName PANEL",
                       "samtools mpileup  -uBf  $ref  $name.recalibrated.bam  > $name.mpileup",
                       "bcftools view -bvcg $name.mpileup > $name.bcf",
                       "bcftools view $name.bcf > $name.raw.MPILEUP.vcf",
                       "$GATK_pre VariantFiltration -R $ref --variant $name.raw.MPILEUP.vcf -o $name.filtered_variants.MPILEUP.vcf --mask $Gene_Panel --maskName PANEL",
                       "cat $name.filtered_variants.MPILEUP.vcf | grep -P \"#+\" > $name.snvs.PANEL.MPILEUP.vcf",
                       "cat $name.filtered_variants.MPILEUP.vcf | grep  -v  \"^#\" | grep -v \"INDEL\" | grep \"PANEL\" >> $name.snvs.PANEL.MPILEUP.vcf",
                       "cat $name.filtered_variants.MPILEUP.vcf | grep -P  \"#+\" > $name.indel.PANEL.MPILEUP.vcf",
                       "cat $name.filtered_variants.MPILEUP.vcf | grep  -v  \"^#\" | grep  \"INDEL\" | grep \"PANEL\" >> $name.indel.PANEL.MPILEUP.vcf",
                       "cat $name.filtered_variants.snvs.GATK.vcf | grep  -P \"#+\" > $name.snvs.PANEL.GATK.vcf",
                       "cat $name.filtered_variants.snvs.GATK.vcf | grep  -v  \"^#\" | grep -v \"GATK_BP_filters\" | grep \"PANEL\">> $name.snvs.PANEL.GATK.vcf",
                       "cat $name.filtered_variants.indel.GATK.vcf | grep  -P \"#+\" > $name.indel.PANEL.GATK.vcf",
                       "cat $name.filtered_variants.indel.GATK.vcf | grep  -v  \"^#\" | grep -v \"GATK_BP_filters\" | grep \"PANEL\" >> $name.indel.PANEL.GATK.vcf",
                       );
                        
                       
                     
                      
print join("\n",@steps);
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

 #Intersection of GATK and MPILEUP snvs files 
                      my $ExecStatement = VcfCommand("$name.snvs.PANEL.MPILEUP.vcf","$name.snvs.PANEL.GATK.vcf","$name.intersect.snvs.vcf","vcf-isec ");
                      LinuxExecute($ExecStatement);
                       #Intersection of GATK and MPILEUP indel files 
                      $ExecStatement = VcfCommand("$name.indel.PANEL.MPILEUP.vcf","$name.indel.PANEL.GATK.vcf","$name.intersect.indel.vcf","vcf-isec ");
                      LinuxExecute($ExecStatement);
                       #Merge all snvs together
                      $ExecStatement = VcfCommand("$name.snvs.PANEL.GATK.vcf.gz","$name.snvs.PANEL.MPILEUP.vcf.gz","$name.All.snvs.vcf", "vcf-merge ");
                       LinuxExecute($ExecStatement);
                       #Merge all indels together
                      LinuxExecute($ExecStatement);
                      $ExecStatement = VcfCommand("$name.indel.PANEL.GATK.vcf.gz","$name.indel.PANEL.MPILEUP.vcf.gz","$name.All.indel.vcf", "vcf-merge ");
                      LinuxExecute($ExecStatement); 
                       #"vcf-merge   $name.snvs.PANEL.MPILEUP.vcf.gz | bgzip -c > $name.All.snvs.vcf",
                       #"vcf-merge  $name.indel.PANEL.GATK.vcf.gz $name.indel.PANEL.MPILEUP.vcf.gz | bgzip -c > $name.All.indel.vcf.gz",
                       #"tabix -p vcf $name.All.snvs.vcf.gz  && tabix -p vcf $name.All.indel.vcf.gz",
                       
                       #Find snvs that are not in the intersection file
                      $ExecStatement = VcfCommand("$name.All.snvs.vcf.gz","$name.intersect.snvs.vcf.gz","$name.FromOnlyOneCaller.snvs.vcf","vcf-isec -c -f ");
                      LinuxExecute($ExecStatement);
                       #Find indels that are not in the intersection file
                      $ExecStatement = VcfCommand("$name.All.indel.vcf.gz","$name.intersect.indel.vcf.gz","$name.FromOnlyOneCaller.indel.vcf","vcf-isec -c -f ");
                      LinuxExecute($ExecStatement); 
                       #"vcf-isec -c $name.All.snvs.vcf.gz $name.merge.snvs.vcf.gz > $name.ORPHAN.snvs.vcf",
                       #"vcf-isec -c $name.All.indel.vcf.gz $name.merge.indel.vcf.gz > $name.ORPHAN.indel.vcf",
                       
                      $ExecStatement ="ls *.gz | xargs -I \{\} bgzip -d \{\}";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="rm *.tbi";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/convert2annovar.pl -format vcf4 $name.intersect.snvs.vcf > $name.intersect.snvs.annovar.input.txt";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/convert2annovar.pl -format vcf4 $name.intersect.indel.vcf > $name.intersect.indel.annovar.input.txt";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/convert2annovar.pl -format vcf4 $name.FromOnlyOneCaller.snvs.vcf > $name.FromOnlyOneCaller.snvs.annovar.input.txt";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/convert2annovar.pl -format vcf4 $name.FromOnlyOneCaller.indel.vcf > $name.FromOnlyOneCaller.indel.annovar.input.txt";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/summarize_annovar.pl  $name.intersect.snvs.annovar.input.txt -buildver hg19 -outfile $name.intersect.snvs -verdbsnp 137 -ver1000g 1000g2012apr -veresp 6500 -alltranscript -remove /miseqdata/tools/annovar/humandb/";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/summarize_annovar.pl  $name.intersect.indel.annovar.input.txt -buildver hg19 -outfile $name.intersect.indel -verdbsnp 137 -ver1000g 1000g2012apr -veresp 6500 -alltranscript -remove /miseqdata/tools/annovar/humandb/";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/summarize_annovar.pl  $name.FromOnlyOneCaller.snvs.annovar.input.txt -buildver hg19 -outfile $name.FromOnlyOneCaller.snvs -verdbsnp 137 -ver1000g 1000g2012apr -veresp 6500 -alltranscript -remove /miseqdata/tools/annovar/humandb/";
                      LinuxExecute($ExecStatement);
                      $ExecStatement ="perl /miseqdata/tools/annovar/summarize_annovar.pl  $name.FromOnlyOneCaller.indel.annovar.input.txt -buildver hg19 -outfile $name.FromOnlyOneCaller.indel -verdbsnp 137 -ver1000g 1000g2012apr -veresp 6500 -alltranscript -remove /miseqdata/tools/annovar/humandb/";
                      LinuxExecute($ExecStatement);

};
   #system ( "mkdir $name; mv $name.* $name;" ); 


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
sub LinuxExecute()
{
	my $ExecStatement = shift;
	say($ExecStatement);
    system ( $ExecStatement);
}

sub AnyMutationsInVCF_File
{
	my $file = shift;
	if ($file =~ /.gz$/)
		{
		#my ($file2) = ($file =~ /(.+)\.gz$/);
		return qx/bgzip -cd  $file | grep -vP "^#+" | wc -l /;
		}
	 else
	 {
		return qx/cat $file | grep -vP "^#+" | wc -l /;
	 }
}

sub VcfCommand
{
#Vcftools throws and error if you try to intersect a file with just a header and no records
#If only 1 file has records then just copy that file to the merge name
#If there are records in both file then retirn vcf-merge statement	
	my $VcfFile1 = shift;
	my $VcfFile2 = shift;
	my $OutputFileName = shift;
	my $VcfCmd = shift;
		
print "\n";
say "There ar emutation in file $VcfFile1 ".AnyMutationsInVCF_File($VcfFile1);
say "There ar emutation in file $VcfFile2 ".AnyMutationsInVCF_File($VcfFile2);
	
	if (((AnyMutationsInVCF_File($VcfFile1)+0) and (AnyMutationsInVCF_File($VcfFile2)+0)))
		{
			#both files have mutations
			my $returnString = CheckForCompressionStatus($VcfFile1,$VcfFile2);
            $VcfFile1 =~ s/.gz$// if $VcfFile1 =~ /.gz$/;
            $VcfFile2 =~ s/.gz$// if $VcfFile2 =~ /.gz$/;
            
            $returnString .= " /miseqdata/tools/vcftools_0.1.10/bin/$VcfCmd $VcfFile1.gz $VcfFile2.gz >  $OutputFileName; bgzip -f $OutputFileName; tabix -p vcf $OutputFileName.gz;"; 			 
			return $returnString
     		
		} elsif (AnyMutationsInVCF_File($VcfFile1))
		{
			#Only one file has mutations so copy it, rename it and zip up new file
			return IfCompressedChangeCommand($VcfFile2)." $OutputFileName;".
			"bgzip -f $OutputFileName; tabix -p vcf $OutputFileName.gz; ".
    		CheckForCompressionStatus($VcfFile1,$VcfFile2)
		} elsif (AnyMutationsInVCF_File($VcfFile2))              
		{
			return IfCompressedChangeCommand($VcfFile2)," $OutputFileName;".
			"bgzip -f $OutputFileName; tabix -p vcf $OutputFileName.gz; ".
    		CheckForCompressionStatus($VcfFile1,$VcfFile2)
		} else
		{
			#neither file has records
			IfCompressedChangeCommand($VcfFile2)." $OutputFileName;".
			"bgzip -f $OutputFileName; tabix -p vcf $OutputFileName.gz; ".
    		CheckForCompressionStatus($VcfFile1,$VcfFile2)
		}
}

sub CheckForCompressionStatus
{
	my $VcfFile1 = shift;
	my $VcfFile2 = shift;
	
	if ($VcfFile1 !~ /gz$/ and $VcfFile2 !~ /gz$/)
		{
			return "bgzip -f $VcfFile1; tabix -p vcf $VcfFile1.gz; bgzip -f $VcfFile2; tabix -p vcf $VcfFile2.gz; "
		}
}

sub IfCompressedChangeCommand
{
	my $file = shift;
	if ($file =~ /gz$/)
		{
		return "bgzip -cd $file > "
		}else
		{
		return "cp $file"
		}
	
}
