#!/usr/bin/perl -w
#
# DN 991213,991215 Heavy file useage..
# DN 991216 Incorporated whatOrf(). 
# DN 000117 Disabled some quality features and save-variables to gain speed.
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
# On a personal note, I enjoy working with open software, and if you have interest in developing the 
# program further I will certainly consider releasing the copyright.
#

(@ARGV==2) || die "Usage: pmestLean.pl genomicSeqFastaFileName estSeqFastaFileName matchScriptPath blastBinaryPath\n";

#my $vectorFileName="";             # vector sequences to remove in fasta format

my $genomic_seq_file_name=$ARGV[0];
my $est_seq_file_name=$ARGV[1];       # EST sequences in fasta format
my $binpath = $ARGV[2]; #"/home/daniel/work/gui";
my $ncbiblastpath = $ARGV[3]; # "/usr/remote/ncbiBlast21";

my $shortestOrf=150;

my $shortHit=70;

my @match;

my $nEST=-1;
my @hit;
my @orf;

my @nOrf;

my @orfStart;
my @orfStop;

my @subjectName;
my @score;
my @p;
my @identities;
my @hitLength;
my @freqIds;
my @queryHitBegin;
my @queryHitEnd;
my @queryHitStrand;
my @subjectHitBegin;
my @subjectHitEnd;
my @subjectHitStrand;

# Read in that genomic sequence file (for use with whatorf()) to relieve the disks a bit.. =)
my @genSeq;
my @genName;
getFastaSeq($genomic_seq_file_name, \@genSeq, \@genName);
my $nSeq = @genSeq;

$genomic_db_file_name = $genomic_seq_file_name; # assume seqname eq dbfilebasename)

print "DEBUG pmestMegablast.pl; testing $ncbiblastpath/megablast -d $genomic_db_file_name -i $est_seq_file_name -b 50000 -F N -D 2 | $binpath/ncbiblastparse.pl |\n";

my $status = 1;
open(BLASTOUT, "$ncbiblastpath/megablast -d $genomic_db_file_name -i $est_seq_file_name -b 50000 -F D -D 2 |$binpath/ncbiblastparse.pl |") or $status = 0;
if ($status == 0) {
    print "ERROR in pmestMegablast.pl: $ncbiblastpath/megablast -d $genomic_db_file_name -i $est_seq_file_name -b 50000 -F N -D 2 | $binpath/ncbiblastparse.pl | failed with $!.\nReturning to caller.\n";
    exit 1;
}

while(<BLASTOUT>) {
    if(/^>/) {
	$nEST++;
	$hit[$nEST]=-1;
	($estName[$nEST])=/^>(.+)/;
    } else {
	if($nEST==-1) {
	    die "Error parsing ncbiblastparse.pl output. Fix me!\n";
	}
	$hit[$nEST]++;
	# Read info on this match..
	($subjectName[$nEST][$hit[$nEST]],$score[$nEST][$hit[$nEST]],$p[$nEST][$hit[$nEST]],$identities[$nEST][$hit[$nEST]],$hitLength[$nEST][$hit[$nEST]],$freqIds[$nEST][$hit[$nEST]],$queryHitBegin[$nEST][$hit[$nEST]],$queryHitEnd[$nEST][$hit[$nEST]],$queryHitStrand[$nEST][$hit[$nEST]],$subjectHitBegin[$nEST][$hit[$nEST]],$subjectHitEnd[$nEST][$hit[$nEST]],$subjectHitStrand[$nEST][$hit[$nEST]])=m/Name\s+\=\s*(.+?)\s*,\s+Score\s+\=\s*(\d+),\s+Expect\s+\=\s*([0-9\.e\-\+]+),\s+Ids\s+\=\s+(\d+)\/(\d+)\s*\(([\d\.]+)\),\s+Query hit\s+(\d+)\-(\d+)\s+\((\w+)\),\s+Subject hit\s+(\d+)\-(\d+)\s+\((\w+)\)/;

	print "Hit: $subjectName[$nEST][$hit[$nEST]] $subjectHitBegin[$nEST][$hit[$nEST]]-$subjectHitEnd[$nEST][$hit[$nEST]]($subjectHitStrand[$nEST][$hit[$nEST]]) P= $p[$nEST][$hit[$nEST]],ids= $identities[$nEST][$hit[$nEST]]/$hitLength[$nEST][$hit[$nEST]]($freqIds[$nEST][$hit[$nEST]]%) Query: $queryHitBegin[$nEST][$hit[$nEST]]-$queryHitEnd[$nEST][$hit[$nEST]]($queryHitStrand[$nEST][$hit[$nEST]]) Subject: $subjectHitBegin[$nEST][$hit[$nEST]]-$subjectHitEnd[$nEST][$hit[$nEST]]($subjectHitStrand[$nEST][$hit[$nEST]])\n";
	
	# Hits are classified as "expressed DNA","EST-homologies" or simply false hits. 
	if( ($p[$nEST][$hit[$nEST]]>0.01) or ($freqIds[$nEST][$hit[$nEST]]<0.5) ) {
	    print "$estName[$nEST] on $subjectName[$nEST][$hit[$nEST]] (Shaky hit.)\n";
	    $hitCategory[$nEST][$hit[$nEST]]=3; # A bit shaky...
	} elsif ( ($p[$nEST][$hit[$nEST]]<1e-50) or ($freqIds[$nEST][$hit[$nEST]]>0.97 and (abs($subjectHitEnd[$nEST][$hit[$nEST]]-$subjectHitBegin[$nEST][$hit[$nEST]]) > $shortHit) ) ) {
	    print "$estName[$nEST] on $subjectName[$nEST][$hit[$nEST]] (Real good hit!)\n";
	    $hitCategory[$nEST][$hit[$nEST]]=1; # Parent candidate
	} else {
	    print "$estName[$nEST] on $subjectName[$nEST][$hit[$nEST]] (Pretty similar hit.)\n";
	    $hitCategory[$nEST][$hit[$nEST]]=2; # Homolog
	}

	if(($hitCategory[$nEST][$hit[$nEST]]==1) or ($hitCategory[$nEST][$hit[$nEST]]==2)) {
	    # Usage: whatORF genome.fasta.filename cosmidName hitBegin hitEnd shortestOrf
#      open(WHATORF,"$matchHome/whatORF.pl $genomic_seq_file_name $subjectName[$nEST][$hit[$nEST]] $subjectHitBegin[$nEST][$hit[$nEST]] $subjectHitEnd[$nEST][$hit[$nEST]] $shortestOrf|");
	    
	    # Set to zero initially; there is some (perl?) bug here: $nOrfs->[$j] = 0 in whatOrf, but until
	    # that particular $nOrf[$j] has been set a second time it refuses to return properly. Odd..
	    @nOrf=(0,0,0,0,0,0);
	    
	    whatORF($subjectName[$nEST][$hit[$nEST]],$subjectHitBegin[$nEST][$hit[$nEST]],$subjectHitEnd[$nEST][$hit[$nEST]],$shortestOrf,\@nOrf,\@orfStart,\@orfStop,$nSeq,@genSeq,@genName);
	    
	    # A call-by-reference _might_ have been better.. But, hey, there it is..
	    # (This design really sucks, but its also really late, and way of reusing a piece code made by accident...)

	    $orf[$nEST][$hit[$nEST]]=0;
	    
	    my $orfCount=0;
	    for($j=0;$j<6;$j++) {
		if($nOrf[$j]) {
		    $orf[$nEST][$hit[$nEST]]+=$nOrf[$j];
		    for($k=0;$k<$nOrf[$j];$k++) {
			print "$subjectName[$nEST][$hit[$nEST]] touches ORF $orfStart[$j][$k] - $orfStop[$j][$k] $j\n";
			$orfCount++;
		    }
		}
	    }
	    
	    if($orfCount!=$orf[$nEST][$hit[$nEST]]) {
		die "Weird error. Number of hit orf begins ($orfCount) returned and caught in pmest do not agree with the returned number of orfs ($orf[$nEST][$hit[$nEST]]).\n";
	    }
	}
    }
}

######################################## END MAIN ##########################################


sub whatORF {
    # Note: the reason for not separating the ORF finding and match-finding is
    # really the idea of using a mother-of-all-regexps for this... 
  
  use POSIX qw(floor);
  
  my $cosmidName=shift;
  my $hitBegin=shift(@_)-1;		# That would be blastn-coordinates or some other position
  my $hitEnd=shift(@_)-1;		# starting at 1. We will want to index perl-strings that start at 0.
  my $shortestOrf=shift;

  my $nOrfs=shift(@_);
  my $orfStart=shift(@_);
  my $orfStop=shift(@_);

  my $nSeq =shift;
  my @genSeq=@_[0..($nSeq-1)];
  my @genName=@_[$nSeq..(2*$nSeq-1)];
  
  my @trouble; # indicate start/stop codons obscured by hit start/stop symbols
  
  # What number was that named sequence?
  my $cosmidNr=0;
  foreach(@genName) {
    if($_=~m/$cosmidName/) {
      last;
    }
    $cosmidNr++;
  }
  
  # Rename that sequence for ease of use..
  $cosmidSeq=$genSeq[$cosmidNr];
  
  # Assume: stops are absolute. 
  # => Check what inter-stop regions this sequence hit.
  
  # Translate trinucleuotides -> aa in all six reading frames.
  my @phase;
  my @phaseOffset=(0,1,2,2,1,0);
  my @phaseEndOffset;
  
  $phase[0]=nt2aa($cosmidSeq);
  $phase[1]=nt2aa(substr($cosmidSeq,1));
  $phase[2]=nt2aa(substr($cosmidSeq,2));
  # Same story with the complementary strand.
  my $complementSeq=$cosmidSeq;
  # The complement has complementary bases, and in the 
  # translation perspective, reverse orientation.
 #print "Seq: $complementSeq\n";
  $complementSeq=~tr/atgcATGC/tacgTACG/;
  $phase[3]=nt2aa(substr(join('',reverse(split(/ */,$complementSeq))),2));
  #print "Phase 3: $phase[3]\n";
  $phase[4]=nt2aa(substr(join('',reverse(split(/ */,$complementSeq))),1));
  $phase[5]=nt2aa(join('',reverse(split(/ */,$complementSeq))));
  
  # If its not a pattern matching problem, make it so. =)
  
  #First, check if hitBegin and hitEnd are start or stop themselves?
  for($j=0;$j<6;$j++) {
    
    # How many bases was left over (in the end-of-sequence) when translating?
    # $phaseEndOffset[$j]=length($cosmidSeq)-$phaseOffset[$j]-3*length($phase[0]);
    
    if($j<=2) {
      $hitBeginAaPos[$j]=POSIX::floor( ($hitBegin-1-$phaseOffset[$j])/3 ); # Hit begin is a 1-based 5' coordinate
      $hitEndAaPos[$j]=POSIX::floor( ($hitEnd-1-$phaseOffset[$j])/3 );
    } else {
      # If we read the complement, the begining is in the end, sortof.
      $hitBeginAaPos[$j]=POSIX::floor((length($cosmidSeq)-$hitEnd-1-$phaseOffset[$j])/3); # Hit begin is a 1-based 5' coordinate; AaPos are phase-internal (was -offset compensated.. shouldnt be - a hit position is there no matter what frame it resides in.), 0-based.
      $hitEndAaPos[$j]=POSIX::floor((length($cosmidSeq)-$hitBegin-1-$phaseOffset[$j])/3);
    }
    # print "DEBUG: Phase $j: $hitBeginAaPos[$j] - $hitEndAaPos[$j]\n";    

    $hitBeginAa[$j]=substr($phase[$j],$hitBeginAaPos[$j],1);
    # print "DEBUG: Phase $j hitBeginAa is a $hitBeginAa[$j].\n"; 

    if($hitBeginAa[$j] eq "M" or $hitBeginAa[$j] eq "+") {
      # print "WARNING: Phase $j hitBeginAa is a $hitBeginAa[$j]!\n";
      $trouble[$j]=1;
      # print "DEBUG: Phase $j hitBeginAa is $hitBeginAa[$j]  ";
    }    
    $hitEndAa[$j]=substr($phase[$j],$hitEndAaPos[$j],1);
    # print "DEBUG: Phase $j hitEndAa is a $hitEndAa[$j].\n"; 
    if($hitEndAa[$j] eq "M" or $hitEndAa[$j] eq "+" ) {  
      # print "WARNING: Phase $j hitEndAa is a $hitEndAa[$j]!\n";
      $trouble[$j]=1;
    }
    # Then set the aa that has hitBegin in it to "*"..
    substr($phase[$j],$hitBeginAaPos[$j],1)="*";
    # ..and correspondingly for hitEnd "#".
    substr($phase[$j],$hitEndAaPos[$j],1)="#"; 
  }

  my @orfStartAa;
  my @orfStopAa;
  
  for($j=0;$j<6;$j++) {  
    
    my @putativeOrf;
    
    if( !$trouble[$j] ) {
      # ORF-finding explained.
      # The third (?:) handles hit internal orfs, the second (?:) cares for orfs over the hit start,
      # the fourth (?:) deals with orfs over hit end and the first (?:) meets the criteria of having an orf completely containing the hit.
      #  while($phase[$j]=~m/\+{1}\w*(?:(M\w*\*\w*\++)|\*)[\w\+]*(?:(M\w*\++)|[\w\+]*)[\w\+]*(?:(M\w*e\w*\++)|[\w*\+]*e\w*)\+{1}/g) {

      # Note: for speed, you should probably change \* and \+ to ordinary letters...
      # 000128 Changed rule 1, 2 & 4 early [\w\+]* into \w* - first there is a stop, then any aa BUT a stop or the e.
      # 000131 Introduced non-greedy matching on some non-orf stretches to catch longest ORFs..
      #      while($phase[$j]=~m/(?:\+{1}\w*?(M\w*\*\w*e\w*\+))|(?:\+{1}\w*\*[\w\+]*?(M\w*\+)[\w\+]*?e\w*\+{1})|(?:\+{1}\w*?(M\w*\*\w*\+)[\w\+]*?e\w*\+)|(?:\+\w*\*[\w\+]*?(M\w*e\w*\+))/g) {
      # 000207 Tentative introduction of "hit over two consecutive ORFS in same phase" by modification of rule 4 to be less fuzzy about the hit structure. Changed hit end letter to # to avoid [\w] reading through hit ends. Changed hit-internal rule to avoid reading hit-end.
      # Dumbo - either match once using smart-longest pattern, or match repeatedly using a while(//g). Using latter now.

      while($phase[$j]=~m/(M\w*\*\w*\#\w*\+)|(M\w*\*\w*\+)|(?:\*[\w\+]*?(M\w*\+))|(M\w*\#\w*\+)/g) {
	# Naive approach for getting all possible, not just the longest: back up yourself once every hit *disabled*     
	if($1) {
	  push @putativeOrf,$1;
	  # print "DEBUG: $j spanning: $1\n";
	  #pos($phase[$j])-=(length($1)-1);
	} elsif($2) {
	  push @putativeOrf,$2;
	  # print "DEBUG: $j init: $2\n";
	  # pos($phase[$j])-=(length($2)-1);
	} elsif($3) {
	  push @putativeOrf,$3;
	  # print "DEBUG: $j inter: $3\n";
	  #pos($phase[$j])-=(length($3)-1);
	} elsif($4) {
	  push @putativeOrf,$4;
	  # print "DEBUG: $j end: $4\n";
	  #pos($phase[$j])-=(length($4)-1);
	}      
      }
    } else {
      my $hitRegion;

      # First, locate the open region around the hit.
      # hitBegin/hitEnd eq stopcodon need special treatement.
      if ($hitBeginAa[$j] eq '+') {
	($hitRegion)= $phase[$j]=~m/(\*[\w\+]*?\#\w*\+)/;
      } elsif ($hitEndAa[$j] eq '+') {
	($hitRegion)= $phase[$j]=~m/(\+\w*\*\w*\#)/;
      } else {
	($hitRegion)= $phase[$j]=~m/(\+\w*\*[\w\+]*?\#\w*\+)/;
      }
      # If no open region is found, dump it.
      # I realize there may be more special cases of "open regions" like start of sequence etc.
      # but they are uncommon, and thus ignored. In the long run, cosmids are hopefully assembled into
      # longer stretches of sequence, and thus even fewer of these cases will arise.
      if(!$hitRegion) {
	next;
      } 

      # print "DEBUG: Hit region phase $j: $hitRegion\n";
      
      # Put the initial Aa's back
      $phase[$j]=~s/\*/$hitBeginAa[$j]/;
      $hitRegion=~s/\*/$hitBeginAa[$j]/;
      $phase[$j]=~s/\#/$hitEndAa[$j]/;
      $hitRegion=~s/\#/$hitEndAa[$j]/;

      # print "DEBUG: Resubst it phase $j: $hitRegion\n";
      
      #   $hitRegionStart=index($phase[$j],$hitRegion);
      #   $hitRegionStop=$hitRegionStart+length(hitRegion);
      
      # Then simply get the longest ORF in the hit-region.
      while($hitRegion=~m/(M\w*\+)/g) {
	push @putativeOrf,$1;
	# print "DEBUG: longest from hit region phase $j: $1\n";
      }
    }

    # The length is of importance.
    $nOrfs->[$j]=0;

    foreach(@putativeOrf) {
      if ((3*length($_)) > $shortestOrf) { 
	#print "DEBUG: Putative: $_\n";
	$orfStartAaPos[$j][$nOrfs->[$j]]=index($phase[$j],$_); # index gives 0-based result
	$orfStopAaPos[$j][$nOrfs->[$j]]=$orfStartAaPos[$j][$nOrfs->[$j]]+length($_)-1; # so, stop should also be 0-based.
	# Retransform start and stop coordinates.
	if($j<=2) {
	  $orfStart->[$j][$nOrfs->[$j]]=3*$orfStartAaPos[$j][$nOrfs->[$j]]+$phaseOffset[$j];
	  # Rational for treating start/stop differently: startAA is on correct 1st-in-codon pos   
	  # wheras stopAA should include two extra nt's. 000131 DN
	  $orfStop->[$j][$nOrfs->[$j]]=3*$orfStopAaPos[$j][$nOrfs->[$j]]+2+$phaseOffset[$j];
	} else {
	  $orfStop->[$j][$nOrfs->[$j]]=-3*$orfStartAaPos[$j][$nOrfs->[$j]]+length($cosmidSeq)-1-$phaseOffset[$j];
	  $orfStart->[$j][$nOrfs->[$j]]=-3*($orfStopAaPos[$j][$nOrfs->[$j]]+1)+length($cosmidSeq)-$phaseOffset[$j]; # return 0-based 5'-coordinate
	}
	$nOrfs->[$j]++;
	
	# and just to be repeat-insured:
	
	#if( index($phase[$j],$_,$orfStartAaPos[$j][$nOrfs->[$j]-1]+1) != -1 ) {
	#  die "ERROR: Multiple identical ORFs found in phase $j searching with ORF-pattern $_.\n";
	#}
	
	# Wow, it actually happened! So, let's deal with it. First, get all ORFs that looked like 
	# our hit one. They are as good gene-candidates as the one we came up with in the first place,
	# right? (Otherwise we could just filter them out depending on their actual position..

	while(index($phase[$j],$_,$orfStartAaPos[$j][$nOrfs->[$j]-1]+1) > -1) {
	  $orfStartAaPos[$j][$nOrfs->[$j]]=index($phase[$j],$_,$orfStartAaPos[$j][$nOrfs->[$j]-1]+1);
	  $orfStopAaPos[$j][$nOrfs->[$j]]=$orfStartAaPos[$j][$nOrfs->[$j]]+length($_)-1; # 0 based
	  # print "DEBUG: Another duplicate was found ($orfStartAaPos[$j][$nOrfs->[$j]]-$orfStopAaPos[$j][$nOrfs->[$j]])";
	  if($j<=2) {
	    $orfStart->[$j][$nOrfs->[$j]]=3*$orfStartAaPos[$j][$nOrfs->[$j]]+$phaseOffset[$j];
	    # Rational for treating start/stop differently: startAA is on correct 1st-in-codon pos   
	    # wheras stopAA should include two extra nt's. 000131 DN
	    $orfStop->[$j][$nOrfs->[$j]]=3*$orfStopAaPos[$j][$nOrfs->[$j]]+2+$phaseOffset[$j];
	  } else {
	    # Rational for treating start/stop differently: startAA is on correct 1st-in-codon pos   
	    # wheras stopAA should include two extra nt's. 000131 DN
	    # Nope? 000204 DN
	    $orfStop->[$j][$nOrfs->[$j]]=-3*$orfStartAaPos[$j][$nOrfs->[$j]]-1+length($cosmidSeq)-$phaseOffset[$j]; # return 0-based 5'-coordinate
	    $orfStart->[$j][$nOrfs->[$j]]=-3*($orfStopAaPos[$j][$nOrfs->[$j]]+1)+length($cosmidSeq)-$phaseOffset[$j]; # return 0-based 5'-coordinate
	  }
	  $nOrfs->[$j]++;	  
	}
	# print "DEBUG: Found ORFS. nOrfs[$j] is now $nOrfs->[$j].\n";
	
#	print "$cosmidName $orfStart->[$j][$nOrfs->[$j]-1] - $orfStop->[$j][$nOrfs->[$j]-1] (phase $j)\n"; #Expected output. Change with caution!
	#      print "Found at $orfStart[$j][$nOrfs[$j]-1] - $orfStop[$j][$nOrfs[$j]-1] (phase $j).\nNt-sequence is ",substr($cosmidSeq,$orfStart[$j][$nOrfs[$j]-1],$orfStop[$j][$nOrfs[$j]-1]-$orfStart[$j][$nOrfs[$j]-1]),"\n"; #DEBUG
      } else {
#       print "VERBOSE: Discarded ph $j hit on ",length($_)," bp ORF.\n"; 
      }
    }   
  }
}

sub getFastaSeq {
    my $fastaFileName=shift;
#  print "DEBUG: getFastaSeq tries to open $fastaFileName.\n";
    open(FASTAFILE,"<$fastaFileName") || die "Sequence fasta input file $fastaFileName open failed.\n";
    
    my $seq = shift; # ref
    my $name = shift; # ref
    
    # First, get the sequences
    my $nFASTA=0;
    while(<FASTAFILE>) {
	chop;
	if(/^\>/) {
	    # On fasta description line
	    $nFASTA++;
	    ($$name[$nFASTA-1])=/^\>(.+?)\s*$/; # skip trailing whitespaces
	    $$seq[$nFASTA-1]="";
	} else {
	    # Well, either the input is broken, or this is sequence data. Let us assume a sentient user.. :)
	    # Get all genomic sequence chars into that $seq string..
	    s/[^atgcnxATGCNX]//g;
	    $$seq[$nFASTA-1].=$_;
	}
    }
    
    # Done processing fasta sequence file
    close(FASTAFILE);
}

sub nt2aa {
  $tri=$_[0];
  $aa="";
  # Make nucleotides lowercase, and have all AAs uppercase.
  $tri=~tr/ATCGNX/atcgnx/;

  # DEBUG: Strict??
  for($k=0;$k<(length($tri)/3);$k++) {
    # Make trinucleotide codons... 
    $_=substr($tri,3*$k,3);

    # Genetic standard code according to NCBI
    #   AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    # Starts = ---M---------------M---------------M----------------------------
    # Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    # Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    # Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    # According to Björn Andersson, this should be used instead of the "Protozoan, mold, etc.",
    # that is ONLY used in the Kinetoplast (miniDNAs).

    s/tt[tc]{1}/F/g;
    s/tt[ag]{1}/L/g;
    s/tc[atgcnx]{1}/S/g;
    s/ta[tc]{1}/Y/g;
    s/ta[ag]{1}/+/g;
    s/tg[tc]{1}/C/g;
    s/tga/+/g;
    s/tgg/W/g;
    s/ct[tcagnx]{1}/L/g;
    s/cc[tcagnx]{1}/P/g;
    s/ca[tc]{1}/H/g;
    s/ca[ag]{1}/Q/g;
    s/cg[tcagnx]{1}/R/g;
    s/at[tca]{1}/I/g;
    s/atg/M/g;
    s/ac[tcagnx]{1}/T/g;
    s/aa[tc]{1}/N/g;
    s/aa[ag]{1}/K/g;
    s/ag[tc]{1}/S/g;
    s/ag[ag]{1}/R/g;
    s/gt[tcagnx]{1}/V/g;
    s/gc[tcagnx]{1}/A/g;
    s/ga[tc]{1}/D/g;
    s/ga[ag]{1}/E/g;
    s/gg[atgcnx]{1}/G/g;
    s/[atgc]*[nx]+[atgc]*[nx]*/U/g;

    # Append the aa the codon encoded to the virtual protein.
    $aa.=$_;
  }
  # Remove the extra nucleotides...

  $aa=~s/[atgcnx]//g; 
  return $aa;
}

sub logo {
  print "-oO0O---O0Oo-\n  | |   | |  \n  | |\"\"\"| |  \n  | |O O| |  \n  | | ^ | |  \n   \\ \\0/ /   \n (c)1999-DN  \n";
}
