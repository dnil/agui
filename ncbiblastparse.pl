#!/usr/bin/perl -w
#
# Parse ncbi-blast results, in particular MegaBlast hits
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

$DEBUG = 0;

my $n_hits = 0;
my %hit;

# There is probably som better way to do this, but, hey, it works! =O)
my ($ncbi_reply, $program);

while(<STDIN>) {
    $ncbi_reply .= $_;
    if(/^MEGABLAST/) {
	$program = 'megablast';
    } elsif(/^BLASTN/) {
	$program = 'blastn';
    } elsif(/^BLASTP/) {
	$program = 'blastp';
    } elsif(/^BLASTX/) {
	$program = 'blastx';
    }    
}

# program could be got from first row..
if($program eq 'blastn') {
    parse_blastn_reply(\%hit, $ncbi_reply);
} elsif ($program eq 'blastx') {
    parse_blastx_reply(\%hit, $ncbi_reply);
} elsif ($program eq 'blastp') {
    parse_blastp_reply(\%hit, $ncbi_reply);
} elsif ($program eq 'megablast') {
    parse_megablast_reply(\%hit, $ncbi_reply);
} else {
    print "ERROR: Could not find a parsable NCBI blast format input on stdin.\n";
    exit 1;
}

$n_hits = @{$hit{subjectName}};
# print "Found $n_hits hits.\n";

# init
my $currentQueryName = $hit{queryName}->[0];
print ">$currentQueryName\n";

for($i=0;$i<$n_hits;$i++) {    
    if ($currentQueryName ne $hit{queryName}->[$i]) {
	$currentQueryName = $hit{queryName}->[$i];
	print ">$currentQueryName\n";
    }

    if($hit{expect}->[$i] =~ /^e/) {
	# silly NCBI format output error
	($exp) = $hit{expect}->[$i] =~ /^(e.+)\s*$/;
	$hit{expect}->[$i] = "1$exp";
    }
    print "Name = $hit{subjectName}->[$i]\t, Score = $hit{score}->[$i], Expect = $hit{expect}->[$i], Ids = ",$hit{identities}->[$i],"/",$hit{alignLength}->[$i]," (",$hit{identities}->[$i]/$hit{alignLength}->[$i]*100,"), Query hit $hit{queryHitBegin}->[$i]-$hit{queryHitEnd}->[$i] ($hit{queryHitStrand}->[$i]), Subject hit $hit{subjectHitBegin}->[$i]-$hit{subjectHitEnd}->[$i] ($hit{subjectHitStrand}->[$i])\n";
    
    # NOTE: if reverse mapping to original ESTs is required, the trashed ESTs
    # must be accounted for!   
}

sub parse_megablast_reply {
    my $hit=shift;
    my $ncbi_reply=shift;

    # If HTML-formatted we translate the name tag into a fasta-header, and remove remaining tags
    $ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    $ncbi_reply=~s/\<.+?\>//g;
    $DEBUG && print "DEBUG: $ncbi_reply\n";
    
    # This version ok for multiple queries in one file.
    # Ok for multi line subject-names.

    my $nHits=0;
    my $inHit=0;
    my $inHitName=0;
    my $inQueryName=0;
    my $inDetail=0;
    my $detailField=0;
    
    my @subjectHitBegin;
    my @subjectHitEnd;
    my @subjectHitStrand;
    my @subjectName;
    my @queryHitBegin;
    my @queryHitEnd;
    my @queryHitStrand;
    my $currentQuery;
    my @queryName;
    my @identities;
    my @alignLength;
    my @score;
    my @expect;
    
    my $next_row_is_alignment = 0;

    foreach $_ (split(/\n{1}/,$ncbi_reply)) {
	
	# TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
	
	# Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
	# The rest of the output is kindly enough different.. =)
	if($inQueryName) {
	    if( /\(\s*\d+\s+letters\)/ ) {
		$inQueryName = 0;
	    } else {
		($queryNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		$currentQuery .= " $queryNamePart";
	    }
	}
	if($inHit) {
	    if($inHitName) {
		if( /Length\s+\=\s+\d+/ ) {
		    $inHitName = 0;
		} else {
		    ($subjectNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		    $subjectName[$nHits-1]  .= " $subjectNamePart";
		}
	    }
	    if($inQueryName) {
		if( /\(\s*\d+\s+letters\)/ ) {
		    $inQueryName = 0;
		} else {
		    ($queryNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		    $currentQuery .= " $queryNamePart";
		}
	    }
	    if($inDetail) {
		if($next_row_is_alignment) {
		    $alignment[$nHits-1] .= "$_\n";
		    $next_row_is_alignment = 0;
		}
		if(/^Sbjct\:/) {
		    if($detailField==0) {
			# Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
			($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    } else { 
			($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    }
		    $alignment[$nHits-1] .= "$_\n";
		} elsif(/^Query\:/) {
		    # A new block of hit alignment rows was found.
		    $detailField++;
		    ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    $alignment[$nHits-1] .= "\n$_\n"; #the extra newline only inside alignments, not before!
		    $next_row_is_alignment = 1;
		} elsif(/Score\s{1}\=/) {
		    # Parse of alignment rows found a new hit on the same subject sequence.
		    $subjectName[$nHits]=$subjectName[$nHits-1];
		    $queryName[$nHits]=$currentQuery;
		    $inDetail=0;
		    # Apparently, the ORF-finding thing needs smallest-first type of coordinates?!
		    if($subjectHitBegin[$nHits-1] > $subjectHitEnd[$nHits-1]) {
			$tmp = $subjectHitBegin[$nHits-1];
			$subjectHitBegin[$nHits-1] = $subjectHitEnd[$nHits-1];
			$subjectHitEnd[$nHits-1] = $tmp ;
		    }
		    if($queryHitBegin[$nHits-1] > $queryHitEnd[$nHits-1]) {
			$tmp = $queryHitBegin[$nHits-1];
			$queryHitBegin[$nHits-1] = $queryHitEnd[$nHits-1];
			$queryHitEnd[$nHits-1] = $tmp ;
		    }
		    $nHits++;
		    ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
		    ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
		} elsif(/^\>/) {
		    # Parse of alignment rows found hit on naew subject sequence.
		    $inDetail = 0;
		    $inHitName = 1;
		    # Hits supposedly begin with a row "> FASTA_NAME"
		    ($subjectName[$nHits]) = /^\>(.+?)\s*$/;
		    $queryName[$nHits] = $currentQuery;
		    # Apparently, the ORF-finding thing needs smallest-first type of coordinates?!
		    if($subjectHitBegin[$nHits-1] > $subjectHitEnd[$nHits-1]) {
			$tmp = $subjectHitBegin[$nHits-1];
			$subjectHitBegin[$nHits-1] = $subjectHitEnd[$nHits-1];
			$subjectHitEnd[$nHits-1] = $tmp ;
		    }
		    if($queryHitBegin[$nHits-1] > $queryHitEnd[$nHits-1]) {
			$tmp = $queryHitBegin[$nHits-1];
			$queryHitBegin[$nHits-1] = $queryHitEnd[$nHits-1];
			$queryHitEnd[$nHits-1] = $tmp ;
		    }
		    $nHits++;
		} elsif(/^Query\=/) {
		    # End of this query sequence. End of detail, end of hit.
		    $inDetail = 0;
		    $inHit = 0;
		    ($currentQuery) = /^Query\=\s*(.+)\s*/;
		    $inQueryName = 1;
		} 
		# in detail ends
		#Parse a hit header..
	    } elsif(/Score\s{1}\=/) {
		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		# Syntax of the P value varies btw runs in blastn.. *sigh*
		($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
	    } elsif (/Identities/) {
		($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
	    } elsif (/Strand\s+\=/) {
		($queryHitStrand[$nHits-1],$subjectHitStrand[$nHits-1])=/Strand\s+\=\s+(\w+)\s+\/\s+(\w+)/;
	    } elsif(/^Query\:/) {
		$inDetail=1;
		$detailField=0;
		# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
		($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
		# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
		$alignment[$nHits-1] .= "$_\n";
		$next_row_is_alignment = 1;
	    }   
	} elsif(/^\>/) {
	    # Hits supposedly begin with a row "> FASTA_NAME"
	    ($subjectName[$nHits])=/^\>(.+?)\s*$/;
	    $inHitName = 1;
	    $queryName[$nHits]=$currentQuery;

	    # Apparently, the ORF-finding thing needs smallest-first type of coordinates?!
	    if( $nHits != 0 ) {
		if($subjectHitBegin[$nHits-1] > $subjectHitEnd[$nHits-1]) {
		    $tmp = $subjectHitBegin[$nHits-1];
		    $subjectHitBegin[$nHits-1] = $subjectHitEnd[$nHits-1];
		    $subjectHitEnd[$nHits-1] = $tmp ;
		}
		if($queryHitBegin[$nHits-1] > $queryHitEnd[$nHits-1]) {
		    $tmp = $queryHitBegin[$nHits-1];
		    $queryHitBegin[$nHits-1] = $queryHitEnd[$nHits-1];
		    $queryHitEnd[$nHits-1] = $tmp ;
		}
	    }
	    $nHits++;
	    $inHit=1;
	    $alignment[$nHits-1]="";
	} elsif(/^Query\=/) { # notice: megablast/blastall *WILL* contain multiple queries..
	    ($currentQuery) = /^Query\=\s*(.+)\s*/;	    
	    $inQueryName = 1;
	    # print "DEBUG: Query $currentQuery.\n";
	}
    }

    for($i=0; $i<$nHits;$i++) {
	$hit{subjectName}->[$i] = $subjectName[$i];
	$hit{subjectHitBegin}->[$i] = $subjectHitBegin[$i];
	$hit{subjectHitEnd}->[$i] = $subjectHitEnd[$i];
	$hit{subjectHitStrand}->[$i] = $subjectHitStrand[$i];
	$hit{queryName}->[$i] = $queryName[$i];
	$hit{queryHitBegin}->[$i] = $queryHitBegin[$i];
	$hit{queryHitEnd}->[$i] = $queryHitEnd[$i];
	$hit{queryHitStrand}->[$i] = $queryHitStrand[$i];
	$hit{identities}->[$i] = $identities[$i];
	$hit{alignLength}->[$i] = $alignLength[$i];
	$hit{score}->[$i] = $score[$i];
	$hit{expect}->[$i] = $expect[$i];
	
	$DEBUG && print "DEBUG: $queryName[$i] $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
    }

    if($nHits == 0) {
	print STDERR "WARNING: no hits found when parsing megablast results.\n";
    }
}

sub parse_blastn_reply {
    my $hit=shift;
    my $ncbi_reply=shift;

    # head
    # ref, about query, 
    # areamap
    # Sequences producing significant alignments
    # <PRE>
    #  <a name = > </a><a href>id</a>name
    # Length
    # score, ids, strand  
    # </PRE>
    #<PRE> 
    # Database:
    
    # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
    # but on the other hand it is possible to use previously written WU-blast parsing code..

    $ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    $ncbi_reply=~s/\<.+?\>//g;
    $DEBUG && print "DEBUG: $ncbi_reply\n";

    # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

    my $nHits=0;
    my $inHit=0;
    my $inDetail=0;
    my $detailField=0;

    my @subjectHitBegin;
    my @subjectHitEnd;
    my @subjectHitStrand;
    my @subjectName;
    my @queryHitBegin;
    my @queryHitEnd;
    my @queryHitStrand;
    my $queryName;
    my @identities;
    my @alignLength;
    my @score;
    my @expect;
    
    my $next_row_is_alignment = 0;

    foreach $_ (split(/\n{1}/,$ncbi_reply)) {    

	# TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
	
	# Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
	# The rest of the output is kindly enough different.. =)
	if($inHit) {
	    if($inDetail) { 
		if($next_row_is_alignment) {
		    $alignment[$nHits-1] .= "$_\n";
		    $next_row_is_alignment = 0;
		}
		if(/^Sbjct\:/) {
		    if($detailField==0) {
			# Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
			($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    } else { 
			($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    }
		    $alignment[$nHits-1] .= "$_\n";
		} elsif(/^Query\:/) {
		    # A new block of hit alignment rows was found.
		    $detailField++;
		    ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
		    $alignment[$nHits-1] .= "\n$_\n"; #the extra newline only inside alignments, not before!
		    $next_row_is_alignment = 1;
		} elsif(/Score\s{1}\=/) {
		    # Parse of alignment rows found a new hit on the same subject sequence.
		    $subjectName[$nHits]=$subjectName[$nHits-1];
		    $inDetail=0;
		    $nHits++;
		    ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
		    ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
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
		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		# Syntax of the P value varies btw runs in blastn.. *sigh*
		($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
	    } elsif (/Identities/) {
		($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
	    } elsif (/Strand\s+\=/) {
		($queryHitStrand[$nHits-1],$subjectHitStrand[$nHits-1])=/Strand\s+\=\s+(\w+)\s+\/\s+(\w+)/;
	    } elsif(/^Query\:/) {
		$inDetail=1;
		$detailField=0;
		# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
		($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
		# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
		$alignment[$nHits-1] .= "$_\n";
		$next_row_is_alignment = 1;
	    }   
	} elsif(/^\>/) {
	    # Hits supposedly begin with a row "> FASTA_NAME"
	    ($subjectName[$nHits])=/^\>(.+)/;
	    $nHits++;
	    $inHit=1;
	    $alignment[$nHits-1]="";
	} elsif(/^Query\=/) {
	    # Actually just evaluated once to get query name..
	    ($queryName)=/^Query\=\s*(.+)/; 
#      print "Query $queryName.\n";
	}
    }

    for($i=0; $i<$nHits;$i++) {
	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
    }

    # Then pop up a user interaction window for choosing hits for use in annotation.
    my $qblast_win=$main->Toplevel;
    $qblast_win->title("qblastn at NCBI results");
    $qblast_win->geometry('+300+300');
    $qblast_win->configure(-background=>'linen',-width=>'600');
    
    my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
    my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
    my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
    for($i=0; $i<$nHits;$i++) {
	$hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
	$DEBUG && print "DEBUG: $hitentry\n";
	$qblast_hits_list->insert('end',$hitentry);
    }  

    my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
    my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
    $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
    # scrollbars... connect commands...
    
    my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
    my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
	#$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
	$annotation_nr=annotation_what_nr($annotation,$annotation_id);
	$note=$$annotation{note}->[$annotation_nr];
	foreach $selected ($qblast_hits_list->curselection) {
	    $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitStrand[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ($subjectHitStrand[$selected]) ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
	}
	main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
    },-text=>"Annotate query")->pack(-side=>'left');
    my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
	foreach $selected ($qblast_hits_list->curselection) {
	    $annotation_nr=annotation_what_nr($annotation,$annotation_id);
	    add_blasthit($canvas,$annotation,$sheet,$seq,'blastn',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$queryHitStrand[$selected],$subjectName[$selected], $subjectHitBegin[$selected],$subjectHitEnd[$selected],$subjectHitStrand[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
	}
	  main::level_layout($annotation, 'blast'); 
	  main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..

    })->pack(-side=>'left');

    my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
	# fancy view?
	$annotation_nr=annotation_what_nr($annotation,$annotation_id);

	print "Alignments (qblastn) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

	foreach $selected ($qblast_hits_list->curselection) {
	    print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitStrand[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ($subjectHitStrand[$selected]) ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
	    print "\n$alignment[$selected]\n";
	}
    })->pack(-side=>'left');


    my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

    # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
    # Suggest name accordning to Gen. nomenclature "whitepaper"?

}

sub parse_blastx_reply {
  my $main=shift;
  my $canvas=shift;
  my $sheet=shift;
  my $seq=shift;
  my $annotation=shift;
  my $annotation_id=shift;
  my $ncbi_reply=shift;
  
  # head
  # ref, about query, 
  # areamap
  # Sequences producing significant alignments
  # <PRE>
  #  <a name = > </a><a href>id</a>name
  # Length
  # score, ids, strand  
  # </PRE>
  #<PRE> 
  # Database:
  
  # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
  # but on the other hand it is possible to use previously written WU-blast parsing code..

  $ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
  $ncbi_reply=~s/\<.+?\>//g;
  $DEBUG && print "DEBUG: $ncbi_reply\n";

  # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

  my $nHits=0;
  my $inHit=0;
  my $inDetail=0;
  my $detailField=0;

  my @subjectHitBegin;
  my @subjectHitEnd;
  my @subjectName;
  my @queryHitBegin;
  my @queryHitEnd;
  my @queryHitFrame;
  my $queryName;
  my @identities;
  my @alignLength;
  my @score;
  my @expect;

  my $next_row_is_alignment = 0;
    
  foreach $_ (split(/\n{1}/,$ncbi_reply)) {    

    # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
    
    # Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
    # The rest of the output is kindly enough different.. =)
    if($inHit) {
      if($inDetail) { 
	  if($next_row_is_alignment) {
	      $alignment[$nHits-1] .= "$_\n";
	      $next_row_is_alignment = 0;
	  }
	if(/^Sbjct\:/) {
	  if($detailField==0) {
	    # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
	    ($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;	    
	  } else { 
	    ($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[\w*-]+\s+(\d+)/;
	  }
	  $alignment[$nHits-1] .= "$_\n";
	} elsif(/^Query\:/) {
	  # A new block of hit alignment rows was found.
	  $detailField++;
	  ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[\w*-]+\s+(\d+)/;
	  $alignment[$nHits-1] .= "\n$_\n";
	  $next_row_is_alignment = 1;
	} elsif(/Score\s{1}\=/) {
	  # Parse of alignment rows found a new hit on the same subject sequence.
	  $subjectName[$nHits]=$subjectName[$nHits-1];
	  $inDetail=0;
	  $nHits++;
	  ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
	  ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
	} elsif(/^\>/) {
	  # Parse of alignment rows found hit on new subject sequence.
	  $inDetail=0;
	  # Hits supposedly begin with a row "> FASTA_NAME"
	  ($subjectName[$nHits])=/^\>(.+)/;
	  $nHits++;
	  $alignment[$nHits-1]="";
	} 
	# in detail ends
	#Parse a hit header..     
      } elsif(/Score\s{1}\=/) {
	($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
	# Syntax of the P value varies btw runs in blastn.. *sigh*
	($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
      } elsif (/Identities/) {
	($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
      } elsif (/Frame\s+\=/) {
	($queryHitFrame[$nHits-1])=/Frame\s+\=\s+([\d+-]+)/;
      } elsif(/^Query\:/) {
	$inDetail=1;
	$detailField=0;
	# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
	($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
	# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
	$alignment[$nHits-1] .= "$_\n";
	$next_row_is_alignment = 1;
      }   
    } elsif(/^\>/) {
      # Hits supposedly begin with a row "> FASTA_NAME"
      ($subjectName[$nHits])=/^\>(.+)/;
      $nHits++;
      $inHit=1;
    } elsif(/^Query\=/) {
      # Actually just evaluated once to get query name..
      ($queryName)=/^Query\=\s*(.+)/; 
#      print "Query $queryName.\n";
    }
  }

  for($i=0; $i<$nHits;$i++) {
    $DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
  }

  # Then pop up a user interaction window for choosing hits for use in annotation.
  my $qblast_win=$main->Toplevel;
  $qblast_win->title("qblastx at NCBI results");
  $qblast_win->geometry('+300+300');
  $qblast_win->configure(-background=>'linen',-width=>'600');
  
  my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
  my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
  my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
  for($i=0; $i<$nHits;$i++) {
    $hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
    $DEBUG && print "DEBUG: $hitentry\n";
    $qblast_hits_list->insert('end',$hitentry);
  }  

  my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
  my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
  $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
  # scrollbars... connect commands...
  
  my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
  my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
						     #$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
						     $annotation_nr=annotation_what_nr($annotation,$annotation_id);
						     $note=$$annotation{note}->[$annotation_nr];
						     foreach $selected ($qblast_hits_list->curselection) {
						       $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitFrame[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected]ex Expect=$expect[$selected]\n";
						     }
						     main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
						   },-text=>"Annotate query")->pack(-side=>'left');
  my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
      foreach $selected ($qblast_hits_list->curselection) {
	  $annotation_nr=annotation_what_nr($annotation,$annotation_id);
	  add_blasthit($canvas,$annotation,$sheet,$seq,'blastx',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$queryHitFrame[$selected],$subjectName[$selected], $subjectHitBegin[$selected],$subjectHitEnd[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
      }
	  main::level_layout($annotation, 'blast'); 
	  main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..
  })->pack(-side=>'left');

  my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
      # fancy view?
      $annotation_nr=annotation_what_nr($annotation,$annotation_id);

      print "Alignments (qblastx) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

      foreach $selected ($qblast_hits_list->curselection) {
	  print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitFrame[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
	  print "\n$alignment[$selected]\n";
      }
  })->pack(-side=>'left');

  my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

  # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
  # Suggest name accordning to Gen. nomenclature "whitepaper"?

}

sub parse_blastp_reply {
  my $main=shift;
  my $canvas=shift;
  my $sheet=shift;
  my $seq=shift;
  my $annotation=shift;
  my $annotation_id=shift;
  my $ncbi_reply=shift;
  
  # head
  # ref, about query, 
  # areamap
  # Sequences producing significant alignments
  # <PRE>
  #  <a name = > </a><a href>id</a>name
  # Length
  # score, ids, strand  
  # </PRE>
  #<PRE> 
  # Database:
  
  # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
  # but on the other hand it is possible to use previously written WU-blast parsing code..

  $ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
  $ncbi_reply=~s/\<.+?\>//g;
  $DEBUG && print "DEBUG: $ncbi_reply\n";

  # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

  my $nHits=0;
  my $inHit=0;
  my $inDetail=0;
  my $detailField=0;

  my @subjectHitBegin;
  my @subjectHitEnd;
  my @subjectName;
  my @queryHitBegin;
  my @queryHitEnd;
  my $queryName;
  my @identities;
  my @alignLength;
  my @score;
  my @expect;

  my $next_row_is_alignment = 0;
    
  foreach $_ (split(/\n{1}/,$ncbi_reply)) {    

    # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
    
    # Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
    # The rest of the output is kindly enough different.. =)
    if($inHit) {
      if($inDetail) { 
	  if($next_row_is_alignment) {
	      $alignment[$nHits-1] .= "$_\n";
	      $next_row_is_alignment = 0;
	  }
	if(/^Sbjct\:/) {
	  if($detailField==0) {
	    # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
	    ($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
	  } else { 
	    ($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[\w*-]+\s+(\d+)/;
	  }
	  $alignment[$nHits-1] .= "$_\n";
	} elsif(/^Query\:/) {
	  # A new block of hit alignment rows was found.
	  $detailField++;
	  ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[\w*-]+\s+(\d+)/;
	  $alignment[$nHits-1] .= "\n$_\n";
	  $next_row_is_alignment = 1;
	} elsif(/Score\s{1}\=/) {
	  # Parse of alignment rows found a new hit on the same subject sequence.
	  $subjectName[$nHits]=$subjectName[$nHits-1];
	  $inDetail=0;
	  $nHits++;
	  ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
	  ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
	} elsif(/^\>/) {
	  # Parse of alignment rows found hit on new subject sequence.
	  $inDetail=0;
	  # Hits supposedly begin with a row "> FASTA_NAME"
	  ($subjectName[$nHits])=/^\>(.+)/;
	  $nHits++;
	  $alignment[$nHits-1]="";
	} 
	# in detail ends
	#Parse a hit header..     
      } elsif(/Score\s{1}\=/) {
	($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
	# Syntax of the P value varies btw runs in blastn.. *sigh*
	($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
      } elsif (/Identities/) {
	($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
      } elsif(/^Query\:/) {
	$inDetail=1;
	$detailField=0;
	# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
	($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
	# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
	$alignment[$nHits-1] .= "$_\n";
	$next_row_is_alignment = 1;
      }   
    } elsif(/^\>/) {
      # Hits supposedly begin with a row "> FASTA_NAME"
      ($subjectName[$nHits])=/^\>(.+)/;
      $nHits++;
      $inHit=1;
    } elsif(/^Query\=/) {
      # Actually just evaluated once to get query name..
      ($queryName)=/^Query\=\s*(.+)/; 
#      print "Query $queryName.\n";
    }
  }

  for($i=0; $i<$nHits;$i++) {
    $DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
  }

  # Then pop up a user interaction window for choosing hits for use in annotation.
  my $qblast_win=$main->Toplevel;
  $qblast_win->title("qblastp at NCBI results");
  $qblast_win->geometry('+300+300');
  $qblast_win->configure(-background=>'linen',-width=>'600');
  
  my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
  my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
  my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
  for($i=0; $i<$nHits;$i++) {
    $hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i]  hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
    $DEBUG && print "DEBUG: $hitentry\n";
    $qblast_hits_list->insert('end',$hitentry);
  }  

  my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
  my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
  $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
  # scrollbars... connect commands...
  
  my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
  my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
						     #$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
						     $annotation_nr=annotation_what_nr($annotation,$annotation_id);
						     $note=$$annotation{note}->[$annotation_nr];
						     foreach $selected ($qblast_hits_list->curselection) {
						       $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
						     }
						     main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
						   },-text=>"Annotate query")->pack(-side=>'left');

  my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
      foreach $selected ($qblast_hits_list->curselection) {
	  $annotation_nr=annotation_what_nr($annotation,$annotation_id);
	  add_blasthit($canvas,$annotation,$sheet,$seq,'blastp',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$$annotation{frame}->[$annotation_nr],$subjectName[$selected],$subjectHitBegin[$selected],$subjectHitEnd[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
	  $qblast_win->update();
      }
    main::level_layout($annotation, 'blast'); 
    main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..
      
  })->pack(-side=>'left');

  my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
      # fancy view?
      $annotation_nr=annotation_what_nr($annotation,$annotation_id);

      print "Alignments (qblastp) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

      foreach $selected ($qblast_hits_list->curselection) {
	  print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
	  print "\n$alignment[$selected]\n";
      }
  })->pack(-side=>'left');

  my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

  # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
  # Suggest name accordning to Gen. nomenclature "whitepaper"?

}
