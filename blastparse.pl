#!/usr/bin/perl -w
#
# Daniel Nilsson, 991122
#
# Parse wu-blastn output 
#
# 
# Rewritten 991206 to check EST-hits.
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

#open(PARSE, ">$ARGV[0]") || die "Output file $ARGV[0] open failed\n";

my $nHits=0;
my $inHit=0;
my $inDetail=0;
my $detailField=0;

# Process input..
while(<STDIN>) {
  chop;

  # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO DO CHANGE THEM ALL!
  if($inHit) {
    if($inDetail) { 
      if(/^Sbjct\:/) {
	if($detailField==0) {
	  # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
	  ($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[ATGCNX-]+\s+(\d+)/;
	} else { 
	  ($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[ATGCNX-]+\s+(\d+)/;
	}
      } elsif(/^Query\:/) {
	# A new block of hit alignment rows was found.
	$detailField++;
	($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[ATGCNX-]+\s+(\d+)/;		
      } elsif(/Score\s{1}\=/) {
	# Parse of alignment rows found a new hit on the same subject sequence.
	$subjectName[$nHits]=$subjectName[$nHits-1];
	$inDetail=0;
	$nHits++;
	($score[$nHits-1])=/Score\s+=\s+(\d+)/;
	($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
      } elsif(/^\>/) {
	# Parse of alignment rows found hit on new subject sequence.
	$inDetail=0;
	# Hits supposedly begin with a row "> FASTA_NAME"
	($subjectName[$nHits])=/^\>(.+)/;
	$nHits++;	     
      } 
      # in detail ends
      #Parse a hit header..     
    } elsif(/Score\s{1}\=/) {
      ($score[$nHits-1])=/Score\s+=\s+(\d+)/;
      # Syntax of the P value varies btw runs in blastn.. *sigh*
      ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
    } elsif (/Identities/) {
      ($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
      ($queryHitStrand[$nHits-1],$subjectHitStrand[$nHits-1])=/Strand\s+\=\s+(\w+)\s+\/\s+(\w+)/;
    } elsif(/^Query\:/) {
      $inDetail=1;
      $detailField=0;
      # If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
      ($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[ATGCNX-]+\s+(\d+)/;
      # Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
    }   
  } elsif(/^\>/) {
    # Hits supposedly begin with a row "> FASTA_NAME"
    ($subjectName[$nHits])=/^\>(.+)/;
    $nHits++;
    $inHit=1;
  } elsif(/^Query\=/) {
    # Actually just evaluated once to get query name..
    #($queryName)=/^Query\=\s*(.+)/; 
  }
}

for($i=0;$i<$nHits;$i++) {

  print "Name = $subjectName[$i]\t, Score = $score[$i], P = $p[$i], Ids = ",$identities[$i],"/",$alignLength[$i],", Query hit $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]), Subject hit $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i])\n";

  # NOTE: if reverse mapping to original ESTs is required, the trashed ESTs 
  # must be accounted for!
}

#close(PARSE);

#################################### END MAIN #######################################

