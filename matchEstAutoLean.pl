#!/usr/bin/perl -w
#
# COPYRIGHT NOTICE
#
# (c)2000-2003 All rights are reserved by the author (Daniel Nilsson) of the A GUI program.
#
# The program is provided as-is, without any claims of accuracy or support. The author will
# not be held liable for any damage, mental, reputational or otherwise, caused by using the program.
# You are free to use the program in your non-profit research organisation, provided that you cite the article 
# describing the program in any publications where the program has been used to produce or curate the data.
#
# Scan nt-db for hits with a cDNA seq.
#
# Usage: matchEsts.pl querySeqencesFasta subjectDbName
#
# 991203 DN
# 000117 DN Removed quality features to gain speed. 
# 000128 Modified to take subject db filename
#

#my $genomicDbFileName="/home/daniel/test/genomic/DB3n";
#my $tempHitParseFileName="/home/daniel/test/genomic/temporaryHitParseFile";

(@ARGV==5)||die "Usage: matchEsts.pl querySeqencesFasta subjectDbName pathToBlastparse blastBinaryPath tmpPath\n";

my $est_seq_file_name = $ARGV[0];
my $genomic_db_file_name = $ARGV[1];
my $path= $ARGV[2]; # "/home/daniel/work/gui";
my $blast2path= $ARGV[3]; #"/usr/remote/blast2";
my $tmppath = $ARGV[4];

my $temp_est_query_file_name="$tmppath/singleEstQuery.$$".int(rand 100000).".fasta";

# Get all query sequences and correspondnig quality values
my @est_seq;
my @est_name;
getFastaSeq($est_seq_file_name, \@est_seq, \@est_name);
my $nSeq=@est_seq;

# Pick a single EST..
#my $i=$ARGV[2];

for($i=0;$i<@est_seq;$i++) {

  open(ESTQUERY, ">$temp_est_query_file_name");
  print ESTQUERY ">$est_name[$i]\n";
  print ESTQUERY "$est_seq[$i]\n";
  close(ESTQUERY);

  #  system("/usr/remote/blast2/blastn /home/daniel/test/genomic/DB3n $temp_est_query_file_name -filter dust -echofilter");
  print ">$est_name[$i]\n";
  
  #$queryLength=length($est_seq[$i]);
  #  print "Query length: $queryLength\n";
  # Run blast, roughly parse the output and give it to me! ;) 
  open(BLASTOUT, "$blast2path/blastn $genomic_db_file_name $temp_est_query_file_name B=50000 V=50000 | $path/blastparse.pl |");
  
  while(<BLASTOUT>) {
    chop;
    # Parse blast-parser output. 
  ($subjectName,$score,$p,$identities,$hitLength,$queryHitBegin,$queryHitEnd,$queryHitStrand,$subjectHitBegin,$subjectHitEnd,$subjectHitStrand)=m/Name\s+\=\s+(.+)\s*,\s+Score\s+\=\s+(\d+),\s+P\s+\=\s+([0-9\.e\-\+]+),\s+Ids\s+\=\s+(\d+)\/(\d+),\s+Query hit\s+(\d+)\-(\d+)\s+\((\w+)\),\s+Subject hit\s+(\d+)\-(\d+)\s+\((\w+)\)/;
    # blastparse.pl: print PARSE "Name = $subjectName[$i]\t, Score = $score[$i], P = $p[$i], Ids = ",$identities[$i],"/",$alignLength[$i],", Query hit $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]), Subject hit $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i])\n";
    
    # Average quality over hit sequence.
#    @quals=split(/,+/,$estQual[$i]);
#    if($queryHitStrand eq "Plus") {
#      $averageHitQuality=(sum(@quals[($queryHitBegin-1)..($queryHitEnd-1)]))/($queryHitEnd-$queryHitBegin);
#    } elsif ($queryHitStrand eq "Minus") {
#      $averageHitQuality=(sum(@quals[($queryHitEnd-1)..($queryHitBegin-1)]))/($queryHitBegin-$queryHitEnd);
#    }
    
    $freqIds = $identities/$hitLength;
#    $perr = pErr($averageHitQuality);
    
    # Report results.
    print "Subject=$subjectName\t,Score=$score,P=$p,Ids=",$identities,"/",$hitLength,",($freqIds) Query hit $queryHitBegin-$queryHitEnd ($queryHitStrand), Subject hit $subjectHitBegin-$subjectHitEnd ($subjectHitStrand)\n";
  }
}

unlink($temp_est_query_file_name);

############################################# END MAIN #################################################

sub getFastaSeq {
  my $fastaFileName=shift;
  open(FASTAFILE,"<$fastaFileName") || die "Sequence fasta input file $fastaFileName open failed.\n";

  my @fasta;
  my @name;

  # First, get the sequences
 # print "Reading sequences...\n";
  my $nFASTA=0;
  while(<FASTAFILE>) {
    chop;
    if(/^\>/) {
      # On fasta description line
      $nFASTA++;
      ($name[$nFASTA-1])=/^\>(.+)/;    
      $fasta[$nFASTA-1]="";
    } else {
      # Well, either the input is broken, or this is sequence data. Let us assume a sentient user.. :)
      # Get all genomic sequence chars into that $fasta string..
      s/[^atgcnxATGCNX]//g;
      $fasta[$nFASTA-1].=$_;
    }
  }

  # Done processing fasta sequence file
  close(FASTAFILE);

  return ($nFASTA, @fasta,@name); 
}

sub getFastaQual {
  my $nQual=0;

  my $fastaQualFileName=shift;
  open(QUALFILE, "<$fastaQualFileName") || die "Quality fasta input file $fastaQualFileName open failed\n";

#  print "Reading quality values...\n";
  while(<QUALFILE>) {
    chop;
    if(/^\>/) {
      # On description line.. Name field should be pretty much equal to the fasta file, so ignore it.
      $nQual++;
      $qual[$nQual-1]="";
    } else {
      $qual[$nQual-1].=join(',',split(/\s+/));
      $qual[$nQual-1].=",";
    }  
  }

  # Done processing quality file
  close(QUALFILE);
  
  return($nQual,@qual);
}

sub sum {
  $sum=0;
  # The pop-version looked nice, but choked on qual==0...
  for($n=0;$n<@_;$n++) {
    $sum+=$_[$n];
  }
  return $sum;
}

sub pErr {
  $Q=shift;
  $pErr=10**(-$Q/10);
  return($pErr);
}

