#!/usr/bin/perl -w
#
# 000111 DN
# 
# Usage: glim.pl icmFile seqFastaFile minExonLength
#
# For Glimmer v2.10 
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

# Global debug flag
$DEBUG = 1;

(@ARGV>2)|| die "Usage: glim.pl icmFile seqFastaFile minExonLength glimmerPath\n";

my $icmFileName = $ARGV[0];
my $seqFastaFileName =$ARGV[1];
my $minExonLen = $ARGV[2];
my $glimmer_path = $ARGV[3];  #= "/home/daniel/install/Glimmer2.0";

my $nSeq;
my @genSeq;
my @genName;

@_=getFastaSeq($seqFastaFileName);
$nSeq=shift(@_);
@genSeq=@_[0..($nSeq-1)];
@genName=@_[$nSeq..(2*$nSeq-1)];

for($i=0;$i<$nSeq;$i++) {
    
  # No spaces in filename!      
  # Well, if we have quotes around the filename we can (and must, to get the name correctly parsed back into the gui...)
  #$genName[$i]=~s/ //g;
  
  # Write genomic (cosmid) sequence to file
  open(GENQUERY, ">$genName[$i].fasta.$$.tmp");
  print GENQUERY ">$genName[$i]\n";
  # Enhance readablity - besides, glimmer might have a row-length limit?
  for($j=0;$j<(length($genSeq[$i]));$j+=60) {
    print GENQUERY substr($genSeq[$i],$j,60),"\n";
  }
  close(GENQUERY);

  # Run glimmer on it..  
  $DEBUG && print STDERR "DEBUG: Trying $glimmer_path/glimmer2 \"$genName[$i].fasta.$$.tmp\" \"$icmFileName\" | tee DEBUG.glimmerout |...\n";
  open(GLIMMER,"$glimmer_path/glimmer2 \"$genName[$i].fasta.$$.tmp\" \"$icmFileName\" | tee DEBUG.glimmerout |");
  
  # Right now, we're most interested in the actual gene predictions made:
  
  my $inPutative=0;
  my $nPutative=0;
  my @putativeBegin;
  my @putativeEnd;
  my @frame;
  my @comment;

  while(<GLIMMER>) {
    chop;
    if($inPutative) {
      ($glimDigit,$putativeBegin[$nPutative],$putativeEnd[$nPutative],$frame[$nPutative],$len,$r,$comment[$nPutative])=m/^\s+(\d+)\s+(\d+)\s+(\d+)\s+\[([+-]+\d+)\s+L\=\s*(\d+)\s+r\=([+-.\d]+)\]\s*(.*)/;

      # Add glimmer exon number as a comment (to make the other comments intelligible)
      $comment[$nPutative].="[#$glimDigit] r=$r";
      if($len>$minExonLen) {
	# Translate Glimmer style frame number into something that is easier to handle
	if($frame[$nPutative]<0) {
	  $frame[$nPutative]=3+abs($frame[$nPutative]);
	}	
	$frame[$nPutative]--; # Based on 0
	
	# Assert that reading frame and start/stop positioning are consistent. 
	# Errors here typically pertain to Glimmer percieving DNA as circular. 
	if($frame[$nPutative]<=2) {
	  if ($putativeEnd[$nPutative]<$putativeBegin[$nPutative]) {
	    $DEBUG && print STDERR "DEBUG: reading frame orientation incorrect! Linear genome. frame=$frame[$nPutative] begin=$putativeBegin[$nPutative] end=$putativeEnd[$nPutative]\n";
	    $nPutative--;
	  }
	} elsif($frame[$nPutative]>2) {
	  if ($putativeEnd[$nPutative]>$putativeBegin[$nPutative]) { 
	    $DEBUG && print STDERR "DEBUG: reading frame orientation incorrect! Linear genome. frame=$frame[$nPutative] begin=$putativeBegin[$nPutative] end=$putativeEnd[$nPutative]\n";
	    $nPutative--;
	  }
	  # We might wanna flip begin and end? Actually, flipping those in pmest.pl makes more sense.
	  # Currently done in compare.pl.
	}
	
	$nPutative++;
      } else {
	$DEBUG && print STDERR "DEBUG: Short putative (#$glimDigit; $len bp) discarded.\n";
      }
    } elsif(m/Putative Genes:/) {
      $inPutative=1;
    }
  }
        
  # Output glimmer results...
  print "Putatives in $genName[$i] [start end frame comments]\n";
  
  for($k=0;$k<$nPutative;$k++) {
    print "$putativeBegin[$k]\t$putativeEnd[$k]\t$frame[$k]\t$comment[$k]\n";    
  }

  unlink("$genName[$i].fasta.$$.tmp");
}




################################# END MAIN ###################################

sub getFastaSeq {
  my $fastaFileName=shift;
  open(FASTAFILE,"<$fastaFileName") || die "Sequence fasta input file $fastaFileName open failed.\n";

  my @fasta;
  my @name;

  # First, get the sequences
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

#echo "run Glimmer2"
#clear
#echo "Genome is " $1
#echo "Find non-overlapping orfs in  tmp.coord"
#rm -f tmp.coord
#long-orfs $1 -l | get-putative >tmp.coord
#echo "Extract training sequences to  tmp.train"
#rm -f tmp.train
#extract $1 tmp.coord >tmp.train
#wc tmp.train
#echo "Build interpolated context model in  tmp.model"
#rm -f tmp.model
#build-icm <tmp.train >tmp.model
#echo "Predict genes with Glimmer2 and output putatives to g2.coord"
#rm -f g2.coord
#glimmer2 $1 tmp.model -l | tee glimmer2.out | get-putative >g2.coord3
