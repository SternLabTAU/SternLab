use strict;
#use lib '/sternadi/home/volume1/shared/tools/pipeline';
my $current_version_path = "/sternadi/home/volume1/daniellem1/SternLab/CirSeq/v5.2_power8"
use lib '/sternadi/home/volume1/shared/tools/pipeline/v5.2_power8';
use Create_cmd;

#########################################################################
# Pipeline running of raw read files (fastq.gz) to frequency files
# Begins with an input direcory containing fastq.gz files
# 1. Convert fastq.gz to fasta (mask low Q bases - V3), and
# 2. Split all fasta files into N equally sized smaller fasta files, with ~50K reads in each
# 3. Run formatdb on each of N files above, and blast against ref seq
# 4. run base calling script on each blast file above (output-freq files)
# 5. Join output of all freq files above
## V4: changes made to splitting of files (into more segments); no filtering or masking (these are performed at base_calling)
## V4.1: allow gap skipping. do not call mutations at ends of reads
# V5: changes made to 
#  1) reading .gz files from toFastaAndSplit_v5, (includes IO::Uncompress::Gunzip, and generating multiple quality files)
#  2) running megablast instead of blast
#  3) basecall quality_cutoff = 23; read quality taken from corresponding part's quality file; printing all contributing mutations
# V5.1.1
# added the option for % id for blast
# V5.2.1
# 
#
#
#########################################################################


die "usage pipeline_runner.pl  <input directory with fastq.gz files> <output dir> <reference genome seq (fasta)> <error output file (with path)> <start at step number> <end at step number> <type of input files, optional f if fastq and not zipped files> <refer to gaps? Y/N default Y> <% id for blast, default=85> <Q score (default Q=23)\n
1. Convert fastq.gz to fasta & split all fasta files into N equally sized smaller fasta files (50K reads per file)\n
2. Run formatdb on each of N files above, and blast against ref seq\n
3. run base calling script on each blast file above (output-freq files)\n
4. Join output of all freq files above

" unless (scalar(@ARGV)>=6);

my $in_dir = $ARGV[0];
$in_dir.="/" unless ($in_dir =~ /\/$/);

my $out_dir = $ARGV[1];
$out_dir.="/" unless ($out_dir =~ /\/$/);
unless (-e $out_dir) {system("mkdir $out_dir");}

my $file_to_write_in_dir = $out_dir."fastq_dir";
open FA, ">$file_to_write_in_dir" or die "cannot open file $file_to_write_in_dir\n";
print FA "input dir for fastq files is $in_dir\n";
close FA;


my $ref_genome = $ARGV[2];
unless (-e $ref_genome) {die "error, reference genome $ref_genome does not exist\n"};
my $err_file=$ARGV[3];

my $start_stage=$ARGV[4];
my $end_stage=$ARGV[5];

my $type_file="z"; # default - gzipped fastq files
if (defined $ARGV[6]) {
    $type_file=$ARGV[6];
}
my $do_gaps="Y"; # default - gzipped fastq files
if (defined $ARGV[7]) {
    $do_gaps=$ARGV[7];
}

my $pcID_blast = 85;

if (defined $ARGV[8]) {
    $pcID_blast=$ARGV[8];
}

my $qscore = 23;

if (defined $ARGV[9]) {
    $qscore=$ARGV[9];
}
die "unexpected error, start stage $start_stage is larger than end stage $end_stage\n" if ($start_stage>$end_stage);

my $scripts_dir = $current_version_path
my $ascii_file = "/sternadi/home/volume1/daniellem1/SternLab/pipelineNGS/ascii_table_processed.txt";
my $blast_dir ="/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin";

my $number_of_reads_per_fasta_file=50000;
my $list_files=$out_dir."list.fastq.txt";

my $num_fasta_files_in_this_run=-1;

# change by tal from increments of 60 it won't take that long
my $sleep_quantum=30;
my $sleep_max=1200000; 

open ERR, ">$err_file" or die "cannot open file $err_file\n";
&main;
close ERR;

sub main {
    for my $i($start_stage..$end_stage) {
	my $num_file=&process($i);
    }
}

sub process {
    my $num=shift;
    &splitFastq if ($num==1);
    &blast_fasta if ($num==2);
    &base_call if ($num==3);
     if ($num==4) {
	 &join_files;
	 &wrap_up;
    }
}

sub splitFastq {
    my $cmd="";
    if ($type_file eq "z") {
	$cmd="ls -l $in_dir"."*.gz \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| perl -lane \'s/\.fastq\.gz//\;print;' >$list_files";
    }
    else {
	$cmd="ls -l $in_dir"."*.fastq \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| perl -lane \'s/\.fastq//\;print;' >$list_files";
    }

	print "cmd is:\n\t$cmd";
	system($cmd);

    my $num_files=&get_number_of_fasta_files_there_should_be($list_files);

	print OUT $list_files;
    my $alias="toFastaAndSplit";
    my $cmd_file=$out_dir."splitFastq.cmd";
    my $cmd1="INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_files\)\n";
    my $cmd2="";
    if ($type_file eq "z") {
		$cmd2="perl $scripts_dir/toFastaAndSplit_v5.2.1.pl $in_dir\$INFILE.fastq.gz $out_dir\$INFILE $number_of_reads_per_fasta_file\n";	
    } else {
		$cmd2="perl $scripts_dir/toFastaAndSplit_v5.2.1.pl $in_dir\$INFILE.fastq $out_dir\$INFILE $number_of_reads_per_fasta_file f\n";
	}
    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,4,join(" ",$cmd1,$cmd2));
     $cmd="qsub $cmd_file";
   	 
#Your job-array 8416810.1-4:1 ("toFastaAndSplit") has been submitted
    my $stdin=`$cmd`;
    &sleep_manager(1,$stdin,$alias);

    sleep(60); # to ensure the glob test works... need to wait
    my @fasta_files=glob($out_dir."*fasta");
   	
    if ($num_files != scalar(@fasta_files)) {
	# TODO (bom) change to more informative error message
	print ERR "error in toFastaAndSplit, number of available fasta files in directory $out_dir = ".scalar(@fasta_files).", should be $num_files\n ";
    }
    if (scalar(@fasta_files)==0) {
	die "error in fastq2fasta, no fasta files created\n"; 
    } 

} ##splitFastq

sub blast_fasta {
    my $list_parts_fasta_files = $out_dir."list_parts_fasta_files.txt";
# all files with size >0
    my $glob_test=$out_dir."*part*fasta";
    my @test=glob($glob_test);
  	
    die "Error, no files produced from splitfasta, glob tests is $glob_test\n" if (scalar(@test)==0);
    my $cmd="ls -lrt $glob_test \| awk \'\$5>0\' \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' >$list_parts_fasta_files";
    
    system($cmd);
    $cmd="wc $list_parts_fasta_files \| awk \'{print \$1}\'";
    my $num_files=`$cmd`;

    my $alias="blast";
    my $cmd_file=$out_dir."blastfasta.cmd";
    my $cmd1="INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_parts_fasta_files\)\n";
    
    my $cmd2="$blast_dir/makeblastdb -in $out_dir\$INFILE -dbtype nucl\n";
	
    my $cmd3="$blast_dir/blastn -query $ref_genome -task megablast -db $out_dir\$INFILE -outfmt \"6 sseqid qstart qend qstrand sstart send sstrand length btop\" -num_alignments 1000000 -dust no -soft_masking F -perc_identity $pcID_blast -evalue 1e-7 -out $out_dir\$INFILE.blast";
	
    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,4,join(" ",$cmd1,$cmd2,$cmd3));
    
    $cmd="qsub $cmd_file";
   
    my $stdin=`$cmd`;
    &sleep_manager(3,$stdin,$alias);
    my @files=glob("$out_dir*.fasta");
  	
    my @files_blast=glob("$out_dir*.blast");
  	
    if (scalar(@files)!=scalar(@files_blast) ) {
	print ERR "number of blast output files".scalar(@files_blast) ." not compatible with number of input fasta files".scalar(@files)."\n";
    }
    if (scalar @files_blast == 0) {
	die "error in blastfasta: no blast output files produced\n";
    }

}


sub base_call {
    my $list_blast_results_files = $out_dir."list_blast_results_files.txt";
	my $cmd="ls -lc $out_dir"."*blast \| awk \'\$5>0\' \| awk \'{print \$NF}\' | sort >$list_blast_results_files";
    system($cmd);
	
    my $list_qual_files = $out_dir."list_qual_files.txt";
    my $qcmd="ls -lc $out_dir"."*qual \| awk \'\$5>0\' \| awk \'{print \$NF}\' | sort >$list_qual_files";
    system($qcmd);
    
    $cmd="wc $list_blast_results_files \| awk \'{print \$1}\'";
    my $num_files=`$cmd`;
	
    my $alias="basecall";
    my $cmd_file=$out_dir."basecall.cmd";
    my $cmd1="INFILE\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_blast_results_files\)\n";
    my $cmd2="FASTQ\=\$\(awk \"NR\=\=\$PBS_ARRAY_INDEX\" $list_qual_files\)\n";  
    my $cmd3="perl $scripts_dir/base_call_and_freqs.v5.2.1.pl \$INFILE \$FASTQ $ref_genome \$INFILE\.freqs $do_gaps 2 $qscore \n"; 
    my $mem_request = 8; #originally it was 8 
    Create_cmd::create_cmd_file($cmd_file,$alias,$num_files,$mem_request,join(" ",$cmd1,$cmd2,$cmd3));
    $cmd="qsub $cmd_file";
  
    my $stdin=`$cmd`;
    &sleep_manager(4,$stdin,$alias);

    my @files_blast=glob("$out_dir*.blast");
  
    my @files_freqs=glob("$out_dir*.freqs");
  	
    if (scalar(@files_freqs)!=scalar(@files_blast) ) {
	print ERR "number of blast output files ".scalar(@files_blast) ." not compatible with number of freqs files created: ".scalar(@files_freqs)."\n";
    }
    if (scalar @files_freqs == 0) {
	die "error in blastfasta: no freqs output files produced\n";
    }
}


sub join_files {
    my $list_freqs_files = $out_dir."*part*.freqs";
# list all output blast files with file size > 0
    #my $cmd="ls -l $list_freqs_files \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| head -1 \| awk -F \"\\"."\.\" \'{print \$1}\'";
    my $cmd="ls -l $list_freqs_files \| awk \'{print \$NF}\' \| awk -F \"/\" \'{print \$NF}\' \| head -1 \| awk -F \"\.\" \'{print \$1}\'";
    my $prefix=`$cmd`;

    chomp $prefix;
    $prefix =~ s/\.fastq//;
    $prefix =~ s/\.gz.*//;
    $prefix =~ s/^\.\///;
    $prefix =~ s/\_.+//;

    my $alias="join";
    my $cmd_file=$out_dir."join.cmd";
    my $output_file = $out_dir.$prefix.".freqs";

    my $cmd1="perl $scripts_dir/join_freqs_files.v5.1.pl $out_dir $output_file $do_gaps";
    Create_cmd::create_cmd_file($cmd_file,$alias,1,2,$cmd1);
    $cmd="qsub $cmd_file";
   
    my $stdin=`$cmd`;
    &sleep_manager(5,$stdin,$alias);
    unless (-e $output_file) {
	die "ERROR! Final output file was not created\n";
    }
}

sub get_number_of_fasta_files_there_should_be {
    my $list_file=shift;
    
    if ($num_fasta_files_in_this_run<0) {
	$num_fasta_files_in_this_run=&get_num_fasta_files($list_file);
    }
    return $num_fasta_files_in_this_run;
}


sub get_num_fasta_files {
    my $list_file=shift;
    my $cmd="wc $list_file \| awk \'{print \$1}\'";
    my $num_files=`$cmd`;
    die "error ,get_num_fasta_files returned zero\n" if ($num_files==0);
    return $num_files;
}

sub sleep_manager {
    my $stage_number=shift;
    my $stdin=shift;
    my $job_name=shift;
#4797045.power2.tau.ac.il
    die "unxpected error, cannot find job number, qsub returned $stdin\n" unless ($stdin =~ m/^(\d+).*/);
    my $job_num=$1;
    my $i=0;
    sleep(30);
    while ($i<$sleep_max) {
	if (&test_qstat($job_num)>0)  {
	    print "sleeping $i...";
	    sleep ($i);
	}
	else {
	    print "$job_name done, no more jobs in q\n";
	    return;
	}
	$i+=$sleep_quantum;
    }
    sleep(10); # sleep another ten seconds to ensure that files appear for ls (mounting bug!! on qb3 cluster)
    print "\nexceeded max sleeping time $sleep_max\n";
}

sub test_qstat {
    my $job_num=shift;
# 2016-08-27 tal adapt to new pbs, print all subjobs instead of just the array's name added switch -t
    my $cmd = "qstat -t \| awk \'\$1 ~ \/^$job_num\/' \| wc \| awk \'{print \$1}\'";
    my $num_q_jobs=`$cmd`; chomp $num_q_jobs;
    return $num_q_jobs;

}


sub wrap_up {
    #my $cmd="grep -i killed *.o*";
    #my $res=`$cmd`;
    #if ($res ne "") {
	#print ERR "WARNING!! The following jobs were killed\n";
	#print ERR $res."\n";
    #}
    my $tmp_dir = $out_dir."tmp/";
    system("mkdir $tmp_dir") unless (-e $tmp_dir);
    my $cmd = "mv $out_dir"."*part* "."$tmp_dir";
    system($cmd);
    $cmd = "mv *.o* "."$tmp_dir";
    system($cmd);

    my $mutations_files = $tmp_dir.'*mutations.list';
    my $mutations_outfile = $out_dir.'mutations_all.txt'; 
    system("cat $mutations_files > $mutations_outfile");
}
