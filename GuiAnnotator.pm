#
# GuiAnnotator module for the A GUI program.
#
# The annotator module contains some silly global variables that don't really fit very well elsewhere...
# Also, annotation viewing and merger is dealt with here.
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

package GuiAnnotator;

my ($DEBUG, $DEVEL, $WARNING);

# BEGIN {

use Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(%annotatorName 
	     %annotatorBioframed 
	     %annotatorBrev 
	     %annotatorDirected 
	     %annotatorDisplayLevel 
	     %annotatorHeight 
	     %annotatorOrder 
	     %annotatorDefaultColor 
	     %annotatorPrefs 
	     %annotatorProtein 
	     $AnnotationBoxFrame
	     add_gff_annotator
	     bioframe
	     merge
	     view_annotation view_merged 
	     annotation_what_nr annotation_what_nr_for_name annotation_what_nr_for_uid);
# }

$DEBUG = $main::DEBUG;
$DEVEL = $main::DEVEL;
$WARNING = $main::WARNING;

#use warnings;
#use strict;

# annotation helper subs
sub annotation_what_nr;
sub annotation_what_nr_for_name;
sub annotation_what_nr_for_uid;

# view helpers
sub bioframe;

sub view_merged;
sub view_annotation;

# merge
sub merge;

sub exists_annotator;
sub add_gff_annotator;

our %annotatorDisplayLevel = ();

our %annotatorOrder = (glimmer2 => 2,  # ordering of annotation levels
		       testcode => 3,
		       EST => 5,
		       merge => 10,
		       manual => 9,
		       ESTORF => 4,
		       ORF => 1,
		       ST => 1,
		       polypy => 6,
		       blast => 7,
		       genbankcds => 12,
		       gff =>11, 
		       regexp => 8);

our %annotatorHeight = (glimmer2 => 3,  # level 0 and sum+1 (29) are used for tickmarks and some spacing
			testcode => 3,
			EST => 5,
			merge => 3,
			manual => 3,
			ESTORF => 3,
			ORF => 3,
			ST => 3,
			polypy => 1,
			blast => 4,
			genbankcds => 3,
			gff => 4,
			regexp => 2);

our $AnnotationBoxFrame = 6;

my @annotators = keys %annotatorOrder;
@annotators = sort { $annotatorOrder{$a} <=> $annotatorOrder{$b} } @annotators;

my $currentLevelSum = 1;
foreach my $annotator (@annotators) {
    if($annotator eq 'ST') {
	# these occupy same levels as ORF
    } else {
	$annotatorDisplayLevel{$annotator} = $currentLevelSum;
	# $DEBUG && print "DEBUG: level for $annotator is $annotatorDisplayLevel{$annotator}.\n";
	$currentLevelSum += $annotatorHeight{$annotator};
    }
}
$annotatorDisplayLevel{ST} = $annotatorDisplayLevel{ORF};
#$DEBUG && print "DEBUG: ST level is $annotatorDisplayLevel{ST} since ORF level is $annotatorDisplayLevel{ORF}.\n";

our %annotatorName = (glimmer2 => 'GLIMMER',
		      testcode => 'Testcode',
		      EST => 'EST blastn',
		      merge=> 'Merged from others',
		      manual => 'Manually added annotation',
		      ESTORF => 'ORF hit by EST',
		      ORF => 'Longest ORFs',
		      ST => 'Start and stop codons',
		      polypy => 'Poly pyrimidine stretch',
		      blast => 'Blast hit',
		      genbankcds => 'GenBank CDS',
		      regexp => 'Regular expression',
		      gff => 'General Feature Format');

our %annotatorBrev = (glimmer2 => ['GLIMMER', 0], 
		      testcode => ['Testcode', 0],
		      EST => ['EST',0],
		      merge=> ['Merged',2],
		      manual => ['Manual',0],
		      ESTORF => ['EST ORF',1],
		      ORF => ['ORF',0],
		      ST => ['Start/stop',1],
		      polypy => ['PolyPy',0],
		      blast => ['Blast',0],
		      genbankcds => ['GBK CDS',0],

		      regexp => ['Regexp',0],
		      gff => ['GFF',0]);

# Bioframed 0: level layout (stacking on available frames), 1: bioframed, i e according to reading frame on 6 frames, 2: neither, depending on reading frame field 0-5 framed, 6 undirected - e g manual)
our %annotatorBioframed = (glimmer2 => 1,
		      testcode => 1,
		      EST => 0,
		      merge=> 1,
		      manual => 1,
		      ESTORF => 1,
		      ORF => 1,
		      ST => 1,
		      polypy => 0,
		      blast => 0,
		      genbankcds => 1,
		      regexp => 0,
		      gff => 0);

# 0: no 1: yes 2: conditional on frame
# 2: blast, manual, EST, merge(!?)
our %annotatorDirected = (glimmer2 => 1,
		      testcode => 1,
		      EST => 2,
		      merge=> 2,
		      manual => 2,
		      ESTORF => 1,
		      ORF => 1,
		      ST => 0,
		      polypy => 1,
		      blast => 2,
		      genbankcds => 1,
		      regexp => 1,
		      gff => 2);
# 
our %annotatorProtein = (glimmer2 => 1,
		      testcode => 1,
		      EST => 0,
		      merge=> 2,
		      manual => 2,
		      ESTORF => 1,
		      ORF => 1,
		      ST => 3,
		      polypy => 0,
		      blast => 1,
		      genbankcds => 1,
		      regexp => 0,
		      gff => 2);

our %annotatorDefaultColor = (glimmer2 => '',
			      testcode => '',
			      EST => '',
			      merge=> '',
			      manual => '',
			      ESTORF => '',
			      ORF => '',
			      ST => '',
			      polypy => '',
			      blast => '',
			      genbankcds => '',
			      regexp => '',
			      gff => 'aquamarine3');

our %annotatorPrefs = (glimmer2 => [""], 
		      testcode => [""],
		      EST => [""],
		      merge=> [""],
		      manual => [""],
		      ESTORF => [""],
		      ORF => [""],
		      ST => ["nosizeremove"],
		      polypy => [""],
		      blast => [""],
		      genbankcds => [""],
		      regexp => [""],
		      gff => [""]);

my @gff_colors = ("#d89248c448c4", "#3ce1769bc146", "DarkSlateGray4", "MediumPurple1", "#6da7f1a90362","chocolate2","PaleVioletRed3", "OliveDrab", "DarkGreen");
my $current_color = 0;

# sub annotator_is_protein { return %annotatorProtein{$_} }; etc. etc. would be nice...

sub exists_annotator {
    my $name = shift;

    if(defined ($annotatorName{$name})) {
	return 1; # exists
    } else {
	return 0;
    }
    
}

sub add_gff_annotator { # cheating here - should really generalise this...
    my $sheet = shift;
    my $name = shift;
    my $height = shift;
    my $bioframed = shift; # default stacked, can do bioframed...

    if(defined ($annotatorName{$name})) {
	$DEBUG && $WARNING && print "WARNING: an annotator by the name of $name (".$annotatorName{$name}.") already exists. NOT adding..\n";
	return 0;
    }

    $annotatorPrefs{$name} = [""];
    $annotatorHeight{$name} = $height;
    $annotatorName{$name} = $name;
    $annotatorBrev{$name} = [$name, 0];
    $annotatorDirected{$name} = 2;
    $annotatorProtein{$name} = 2;
    $annotatorBioframed{$name} = $bioframed;

    if($name eq 'gff_CDS' ) {
	$annotatorDefaultColor{$name} = 'DarkSlateGray4';
    } elsif ($name eq 'gff_uORF') {
	$annotatorDefaultColor{$name} = 'MediumPurple1';
    } elsif ($name eq 'gff_pPyAg') {
	$annotatorDefaultColor{$name} = '#6da7f1a90362';
    } else {
	$annotatorDefaultColor{$name} = $gff_colors[$current_color];

	$current_color++;
	if($current_color == @gff_colors) {
	    $current_color = 0;
	}
    }

    my @annotators = keys %annotatorOrder;
    @annotators = sort { $annotatorOrder{$a} <=> $annotatorOrder{$b} } @annotators;
    foreach my $annotator (@annotators) {
	if($annotatorOrder{$annotator} >= $annotatorOrder{gff} && $annotator ne 'gff') {
	    $annotatorOrder{$annotator}++;
	}
    }
    $annotatorOrder{$name} = $annotatorOrder{gff};

    $$sheet{"display_${name}_level"} = 1;

    $WARNING && print "Adding a new gff annotator: name ".$annotatorName{$name}."\n";

    main::set_main_view_height($sheet);
    ${main::view_m}->invoke("Refresh");

    return 1;
}

sub bioframe {
    my $frame=shift;
    my $type=shift;

    if ($type eq 'EST' or $type eq 'polypy') {
	if ($frame < 3) {
	    return "+";
	} elsif ($frame < 6) {
	    return "-";
	} elsif ($frame == $AnnotationBoxFrame) {
	    return "unoriented";
	}
    } else {
	if ($frame < 3) {
	    return ("+".($frame+1));
	} elsif ($frame < 6) {
	    return ("-".(6-$frame));
	} elsif ($frame == $AnnotationBoxFrame) {
	    return "unoriented";
	}
    }
}

sub view_merged {
  my $main=shift;
  my $annotation=shift; # ref
  my $annotation_id=shift;
  my $status_text=shift; # ref

  my $annotation_nr=annotation_what_nr($annotation,$annotation_id);
  if ($annotation_nr!=-1) {
      my $bioframe = bioframe($$annotation{frame}->[$annotation_nr],$$annotation{type}->[$annotation_nr]);
  
      print "View (merged) $$annotation{id}->[$annotation_nr]\n$$annotation{seqName}->[$annotation_nr]\n";
      print "Made with $annotatorName{$$annotation{type}->[$annotation_nr]}.\n";
      print "$$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] ($bioframe) ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.\n";
      $DEVEL && print "DEVEL: annotation uid $$annotation{uid}->[$annotation_nr].\n";
      $DEVEL && defined($$annotation{frame}->[$annotation_nr]) && print "DEVEL: Internal frame $$annotation{frame}->[$annotation_nr]\n";
      if(defined($$annotation{name}->[$annotation_nr])) {
	  print "Name: $$annotation{name}->[$annotation_nr]\n";
      }
      if(defined($$annotation{comment}->[$annotation_nr])) {
	  print "Comment: $$annotation{comment}->[$annotation_nr]\n";
      }
      if(defined($$annotation{note}->[$annotation_nr])) {
	  print "Note: $$annotation{note}->[$annotation_nr]\n";
      }
      $$status_text="$annotatorName{$$annotation{type}->[$annotation_nr]}: $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] ($bioframe) ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.";
      if(defined($$annotation{name}->[$annotation_nr])) {
	  $$status_text .= " $$annotation{name}->[$annotation_nr]";
      }
      if (defined($$annotation{comment}->[$annotation_nr])) {
	  $$status_text .= " $$annotation{comment}->[$annotation_nr]";
      }
  } else {
      print "Error! No such annotation id in annotations.\n";
  }
}

sub view_annotation {
  my $main=shift;
  my $annotation=shift;
  my $annotation_id=shift;
  my $status_text=shift;

  # Find annotation of that id..
  my $annotation_nr=annotation_what_nr($annotation,$annotation_id);

  if ($annotation_nr ==-1) {
      $DEBUG && print "DEBUG: No such annotation $annotation_id.\n";
      return;
  }

  print "View ($$annotation{type}->[$annotation_nr]) $$annotation{id}->[$annotation_nr]\n";
  defined($$annotation{seqName}->[$annotation_nr]) && print "On sequence $$annotation{seqName}->[$annotation_nr]\n";
  print "Made with $annotatorName{$$annotation{type}->[$annotation_nr]}.\n";

  my $bioframed = $annotatorBioframed{$$annotation{type}->[$annotation_nr]};

  if ($$annotation{type}->[$annotation_nr] eq "ST") {
      $$status_text = "$annotatorName{$$annotation{type}->[$annotation_nr]}.";
      return; 
  }
  
  my $bioframe;
  if($bioframed) {
     $bioframe = bioframe($$annotation{frame}->[$annotation_nr],$$annotation{type}->[$annotation_nr]);
      print "$$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] ($bioframe) ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.\n";
  } else {
      print "$$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr]  ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.\n";
  }
  $DEVEL && print "DEVEL: annotation uid $$annotation{uid}->[$annotation_nr].\n";

  $DEVEL && defined($$annotation{frame}->[$annotation_nr]) && print "DEVEL: Internal frame $$annotation{frame}->[$annotation_nr]\n";  
  if($bioframed) {
      $$status_text="$annotatorName{$$annotation{type}->[$annotation_nr]}: $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] ($bioframe) ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.";
  } else {
      $$status_text="$annotatorName{$$annotation{type}->[$annotation_nr]}: $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] ".($$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1)."bp.";
  } 

  if(defined($$annotation{name}->[$annotation_nr])) { 
      print "Name: $$annotation{name}->[$annotation_nr]\n";
      $$status_text .= " $$annotation{name}->[$annotation_nr]";
  }
  if(defined($$annotation{comment}->[$annotation_nr])) { 
      print "Comment: $$annotation{comment}->[$annotation_nr]\n";
      $$status_text .= " $$annotation{comment}->[$annotation_nr]";
  }
  if(defined($$annotation{note}->[$annotation_nr])) {
      if($$annotation{type}->[$annotation_nr] eq "blast") {
	  print "Note (alignment):\n$$annotation{note}->[$annotation_nr]\n";
      } else {
	  print "Note:\n$$annotation{note}->[$annotation_nr]\n";
      }
  }
}

# annotation helpers

sub annotation_what_nr {
    my $annotation=shift;
    my $annotation_id=shift;

    if(!defined($annotation_id) or $annotation_id eq '') {
	$DEBUG && print "DEBUG: called annotation_what_nr with blank or undef annotation_id ($annotation_id).\n";
	return -1;
    }

    # O(1) - small enough, and well used enough to justify storage and rewrite
    if(defined($main::annotation_nr_cache{$annotation_id})) {	 
	return $main::annotation_nr_cache{$annotation_id};
    } else {	
	$WARNING && print "WARNING: annotation id $annotation_id not in annotation_nr_cache!\n";
	# Error - not found!
	# Find annotation of that id.. worstcase O(n)..
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{id}->[$i] eq $annotation_id) {
		# Found it, so there is no need to look any further..
		$WARNING && print "WARNING: annotation id $annotation_id not in annotation_nr_cache, but was found in the annotation hash. Please check the cache handling!\n";
		return $i;
	    }
	}
	$WARNING && print "WARNING: annotation id $annotation_id not in annotation_nr_cache, nor in the annotation hash. Please check where the calling routine got that id. Incorrect deletion or entry?\n";
	return -1;
    }
}

sub annotation_what_nr_for_name {
    my $annotation=shift;
    my $annotation_name=shift;
    
    # Find annotation of that id..
    for(my $i=0;$i<$$annotation{nr};$i++)  {
	if($$annotation{name}->[$i] eq $annotation_name) {
	    # Found it, so there is no need to look any further..
	    return $i;
	}
    }
    # Error - not found!
    return -1;
}

sub annotation_what_nr_for_uid {
    my $annotation = shift;
    my $annotation_uid = shift;
    
    # Find annotation of that id..
    for(my $i=0;$i<$$annotation{nr};$i++) {
	if($$annotation{uid}->[$i] == $annotation_uid) {
	    # Found it, so there is no need to look any further..
	    return $i;
	}
    }
    # Error - not found!
    return -1;
}

sub merge {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $merge_status=shift;
    my $sheet=shift;
    my $seqName=shift;   
    
    my @mergeable_types = ("glimmer2", "testcode", "ESTORF", "ORF", "manual", "blast", "genbankcds");
    # first pass - mark annotations for merger, and don't list doubles.

    my %e_list=(); # Export list
    my $ex_nr=0;

    my $nr=$$annotation{nr};

    # Check for a job to do before we invest any cycles in this... =)
    if(!($$merge_status{glimmer2}||$$merge_status{testcode}||$$merge_status{ESTORF}||$$merge_status{ORF}||$$merge_status{manual}||$$merge_status{blast}||$$merge_status{genbankcds}) ) {
	print "No annotations were selected for export.\n";
	return;
    } else {
	$DEBUG && print "DEBUG: Merging... ($nr annotations)\n";
    }
    
    foreach my $atype (@mergeable_types) {
	if($$merge_status{$atype}) {
	    for(my $i=0;$i<$nr;$i++) {
		if ($$annotation{type}->[$i] eq $atype) {
		    $DEVEL && print "DEVEL: Checking $atype hit $$annotation{id}->[$i]...\n";
		    
		    my $present=0;
		    my $dupe_nr;
		    
		    if($ex_nr==0) {
			# $DEBUG && print "DEBUG: First annotation!\n";
			# premium merger
			foreach my $field (keys %$annotation) {
			    if($field eq "note") {
				if($$annotation{note}->[$i]) {
				    # experimenting to see if we could get something like what Esteban wants..
				    # this is an irregular operation, so we can affort to spend quite a few cycle.
				    $e_list{note}->[$ex_nr] = "*** NOTE FROM $atype ***\n".$$annotation{note}->[$i];
				}
			    } else {
				$e_list{$field}->[$ex_nr] = $$annotation{$field}->[$i];
			    }
			}
			$e_list{$atype}->[$ex_nr]=1;
			$ex_nr++;
		    } else {
			for(my $j=0;$j<$ex_nr;$j++) {
			    # Is there any other annotation already in the export list with these coordinates?
			    if( ($$annotation{start}->[$i] == $e_list{start}->[$j]) && ($$annotation{stop}->[$i] == $e_list{stop}->[$j]) && ($$annotation{frame}->[$i] == $e_list{frame}->[$j])) {
				$present=1;
				$dupe_nr=$j;
			    }
			}

			if($present) {
			    $DEVEL && print "DEVEL: ...present in export list already.\n";
			    # Ok, this was already there - lets note that we had a  hit!
			    $e_list{$atype}->[$dupe_nr]=1;
			    if($$annotation{note}->[$i]) {
				$e_list{note}->[$dupe_nr] .= "*** NOTE FROM $atype ***\n".$$annotation{note}->[$i];
				$DEVEL && print "DEVEL: note for nr $dupe_nr from nr $i (type $atype) is ".$$annotation{note}->[$i]."\n";
			    }
			    if($$annotation{comment}->[$i]) {
				$e_list{note}->[$dupe_nr] .= "*** COMMENT FROM $atype ***\n".$$annotation{comment}->[$i]."\n";
			    }
			} else {
			    # No previous annotation with the same coordinates!
			    $DEVEL && print "DEVEL: Novel coordinates - adding as tentative export candidate.\n";
			    # copy an annotation. M-% object-orient-buffer<ret> =]
			    foreach my $field (keys %$annotation) {
				if($field eq "comment") {
				    if($$annotation{comment}->[$i]) {
					$e_list{note}->[$ex_nr] .= "*** COMMENT FROM $atype ***\n".$$annotation{comment}->[$i]."\n";
				    }
				} elsif($field eq "note") {
				    if($$annotation{note}->[$i]) {
					$e_list{note}->[$ex_nr] .= "*** NOTE FROM $atype FEATURE ***\n".$$annotation{note}->[$i];
				    }
				} else {
				    $e_list{$field}->[$ex_nr]=$$annotation{$field}->[$i];
				}
			    }
			    $e_list{$atype}->[$ex_nr]=1;
			    $ex_nr++;
			}
		    }
		}
	    }
	}
    }

    # second pass - actually export

    if($$merge_status{and}) {
	for(my $j=0;$j<$ex_nr;$j++) {
	    # We know that at least one of the annotations has been selected, so we can check for 
	    # falsifications of the and criteria (recall the shortcircuit property of perl logic-ops)
	    if($$merge_status{glimmer2} && !$e_list{glimmer2}->[$j]) {
	    } elsif($$merge_status{testcode} && !$e_list{testcode}->[$j]) {
	    } elsif($$merge_status{ESTORF} && !$e_list{ESTORF}->[$j]) {
	    } elsif($$merge_status{ORF} && !$e_list{ORF}->[$j]) {
	    } elsif($$merge_status{blast} && !$e_list{blast}->[$j]) {
	    } elsif($$merge_status{manual} && !$e_list{manual}->[$j]) {
	    } elsif($$merge_status{genbankcds} && !$e_list{genbankcds}->[$j]) {
	    } else {
		$nr=$$annotation{nr};
		# Alright! This annotation goes on export.
		$DEBUG && print "DEBUG: Merging elist nr $j as new annotation nr $nr.\n";
		# Account!!
		foreach my $field (keys %$annotation) {
		    $$annotation{$field}->[$nr]=$e_list{$field}->[$j];
		}

		$$annotation{uid}->[$nr] = $$annotation{unr};
		$$annotation{unr}++;

		$$annotation{id}->[$nr]="mb_$nr";
		$main::annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;
		$$annotation{type}->[$nr]="merge";
		$$annotation{color}->[$nr]="purple";

		if($e_list{frame}->[$j]<3) {
		    $$annotation{level}->[$nr]=$$annotation{frame}->[$nr];
		} elsif ($e_list{frame}->[$j]>2) {
		    $$annotation{level}->[$nr]=$$annotation{frame}->[$nr]-5;
		}

		$$annotation{comment}->[$nr] = "Merged with and from ";
		foreach my $allowed_type (@mergeable_types) {
		    if($$merge_status{$allowed_type} && $e_list{$allowed_type}->[$j]) {
			$$annotation{comment}->[$nr] .= "$allowed_type ";
		    }
		}
		$$annotation{comment}->[$nr] .= ".";

		# This is an and-exactly type of annotation. The notes and comments of the individual merged items are concatenated and added to the new annotation.
		# Filtering of this accordning to what a human annotator would do should be attempted!

		$$annotation{nr}++;
	    }
	}
    } elsif ($$merge_status{or}) {
	for(my $j=0;$j<$ex_nr;$j++) {
	    if(($$merge_status{glimmer2} && $e_list{glimmer2}->[$j])|| ($$merge_status{testcode} && $e_list{testcode}->[$j]) || ($$merge_status{ESTORF} && $e_list{ESTORF}->[$j]) || ($$merge_status{ORF} && $e_list{ORF}->[$j]) || ($$merge_status{manual} && $e_list{manual}->[$j]) || ($$merge_status{blast} && $e_list{blast}->[$j]) || ($$merge_status{genbankcds} && $e_list{genbankcds}->[$j])) {
		$nr=$$annotation{nr};
		# Alright! This annotation goes on export.
		$DEBUG && print "DEBUG: Merging elist nr $j as new annotation nr $nr.\n";
		foreach my $field (keys %$annotation) {
		    $$annotation{$field}->[$nr]=$e_list{$field}->[$j];
		}

		$$annotation{uid}->[$nr] = $$annotation{unr};
		$$annotation{unr}++;

		$$annotation{id}->[$nr]="mb_$nr";
		$main::annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;
		$$annotation{type}->[$nr]="merge";
		$$annotation{color}->[$nr]="purple";

		$$annotation{comment}->[$nr] = "Merged with or from ";
		foreach my $allowed_type (@mergeable_types) {
		    if($$merge_status{$allowed_type} && $e_list{$allowed_type}->[$j]) {
			$$annotation{comment}->[$nr] .= "$allowed_type ";
		    }
		}
		$$annotation{comment}->[$nr] .= ".";

		if($e_list{frame}->[$j]<3) {
		    $$annotation{level}->[$nr]=$$annotation{frame}->[$nr];
		} elsif ($e_list{frame}->[$j]>2) {
		    $$annotation{level}->[$nr]=$$annotation{frame}->[$nr]-5;
		}
		$$annotation{nr}++;
	    }
	}
    } elsif ($$merge_status{comp}) {
	# Run in entirely separate track, or run within the "second pass"?
	# Pro 2nd: only "singlet" annotations of interesting types in e_list.
	# Con 2nd: potentially ineffective. 
        # Alternatives: * mark olaps in a string/array: $hit{glimmer2}->[sequencepos] and extract regions
	# fulfilling the selected criteria; then export annotations overlapping these. M+N+Ñ^2 
	# * ORFbased; note overlaps from other types on each e_listed annotation. O(N^2) expensive.
	for(my $j=0;$j<$ex_nr;$j++) {
	    # importing compare. Ignoring hit-classes; that is already implemented in removeAll.
	    my $olap=0;
	    #if($$merge_status{glimmer2} && !$e_list{glimmer2}->[$j]) {
	    # Not a valid exact and.
	    #} elsif($$merge_status{ESTORF} && !$e_list{ESTORF}->[$j]) {
	    # Not a valid exact and.
	    #} elsif($$merge_status{ORF} && !$e_list{ORF}->[$j]) {	
	    # Not a valid exact and.
	    #} else {
	    #$DEBUG && print "DEBUG: Exact match. \n";
	    #$olap=10;
	    #}
	    for(my $k=$j;$k<$ex_nr;$k++) {
		# Check against remaining annotations on list for olaps, starting with this in case it 
		# has had an exact match.

		# overlaps:
		# Case #0
		# ------------
		# ++++++++++++
		# Case #1                   #2               #3                  #4
		# -----------------    ------------   ------------            ---------
		#     ++++++++        ++++++++++            +++++++++++    +++++++++++++++++
		# Case #5                   #6               #7                  #8
		# -----------------   ----------------  ---------              -----------
		# +++++++++++               ++++++++++  +++++++++++++++   ++++++++++++++++

		# Obviously enough, exact matches are overlaps.
		# Check match-length for "shorter-than-interesting". Pretty uninteresting, though.
		my @notMyType;

	      TYPE:      foreach my $type (@mergeable_types) {

		  # The merge_status check is actually superfluous right now. (Done by the export routine.)
		  #  
		  if(($e_list{start}->[$k]==$e_list{start}->[$j])&&($e_list{stop}->[$k]==$e_list{stop}->[$j])) {
		      if($j == $k) {
			  # Ok, this is the hit itself. Any dupe of the current type?
			  if($e_list{$type}->[$k] && ($e_list{type}->[$j] ne $type)) {
			      # Case 0: exact match; yes, dupe that is not of this
			      $DEBUG && print "DEBUG: Exact match.\n";
			      $olap=10;
			  } else {
			      # Hey, I was compared to myself!
			  }
		      } else {
			  # Case 0: exact match
			  $DEBUG && print "DEBUG: Exact match.\n";
			  $olap=10;
		      }
		  }

		  # Check if the other annotation ($k) is of the correct type (or has dupe of correct type) 
		  # and that this is not (only) a hit between equal types. (non-overlap!)
		  # Dupes of other kinds appear only with _exact_ hits, so checking for "annotation not of this kind"
		  # should be sufficent if we do an exact match check first.
		  if($$merge_status{$type} && ($e_list{$type}->[$k] or ($e_list{type}->[$k] eq $type)) && ($e_list{type}->[$j] ne $type) ) {
		      if (($e_list{start}->[$k] > $e_list{start}->[$j] ) and ( $e_list{stop}->[$k] < $e_list{stop}->[$j])) {
			  # Case 1: hit covers hit
			  # Mark orf as overlapping.
			  $olap=1;
		      } elsif ( ($e_list{start}->[$k] < $e_list{start}->[$j]) && ($e_list{stop}->[$k] < $e_list{stop}->[$j]) && ($e_list{stop}->[$k] > $e_list{start}->[$j] ) ) {
			  if($e_list{stop}->[$k]-$e_list{start}->[$j]>$$merge_status{shortest_olap}) {
			      $DEBUG && print "DEBUG: Type 2 olap detected.\n";
			      $olap=2;
			  } else {
			      $DEBUG && print "DEBUG: Type 2 olap detected but discarded since olap only ",$e_list{stop}->[$k]-$e_list{start}->[$j], " bp.\n";
			  }
		      } elsif ( ($e_list{start}->[$k] > $e_list{start}->[$j] ) && ($e_list{start}->[$k]<$e_list{stop}->[$j] ) && ( $e_list{stop}->[$k]>$e_list{stop}->[$j] ) ) {		
			  if($e_list{stop}->[$j]-$e_list{start}->[$k] > $$merge_status{shortest_olap}) {
			      $DEBUG && print "DEBUG: Type 3 olap detected.\n";
			      $olap=3;
			  } else {
			      $DEBUG && print "DEBUG: Type 3 olap detected but discarded since olap only ",$e_list{stop}->[$j]-$e_list{start}->[$k], " bp.\n";
			  } 
		      } elsif (($e_list{start}->[$k] < $e_list{start}->[$j]) && ( $e_list{stop}->[$k] > $e_list{stop}->[$j] ) ) {
			  $DEBUG && print "DEBUG: Type 4 olap detected.\n";
			  $olap=4;
		      } elsif (($e_list{start}->[$k] ==$e_list{start}->[$j] ) && ( $e_list{stop}->[$k] < $e_list{stop}->[$j] ) ) {
			  $DEBUG && print "DEBUG: Type 5 olap detected.\n";
			  $olap=5;
		      } elsif (($e_list{start}->[$k] > $e_list{start}->[$j] ) && ( $e_list{stop}->[$k]==$e_list{stop}->[$j] ) ) {
			  $DEBUG && print "DEBUG: Type 6 olap detected.\n";
			  $olap=6;
		      } elsif (($e_list{start}->[$k] ==$e_list{start}->[$j] ) && ( $e_list{stop}->[$k] > $e_list{stop}->[$j] ) ) {
			  $DEBUG && print "DEBUG: Type 7 olap detected.\n";
			  $olap=7;
		      } elsif (($e_list{start}->[$k] < $e_list{start}->[$j] ) && ( $e_list{stop}->[$k]==$e_list{stop}->[$j] ) ) {
			  $DEBUG && print "DEBUG: Type 8 olap detected.\n";
			  $olap=8;
		      }
		      if($olap) {
			  if($$merge_status{frame}) {
			      if(!(($e_list{frame}->[$k] < 3 and $e_list{frame}->[$j] < 3) or ($e_list{frame}->[$k] >=3 and $e_list{frame}->[$j]>=3 ) )) {
				  # Mismatched frame orientation
				  $olap=0; # lower flag. Would have been done after "cross-matches" were added anyway.
				  next TYPE;
			      }
			  } elsif($$merge_status{framex}) {
			      if(!($e_list{frame}->[$k] == $e_list{frame}->[$j])) {
				  # Mismatched frame
				  $olap=0; # lower flag. Would have been done after "cross-matches" were added anyway.
				  next TYPE;
			      }
			  }
			  $e_list{"olap_$type"}->[$j]=$olap;
			  foreach my $thisType (@mergeable_types) {
			      if($e_list{$thisType}->[$j]) {
				  $e_list{"olap_$thisType"}->[$k]=10+$olap; # Note other-type
			      }
			  }
			  # lower flag so as not to "add" the remaining annotations in the inner loop..
			  $olap=0;
		      }
		  }
	      }
	    }

	    foreach my $thisType (@mergeable_types) {
		if($e_list{"olap_$thisType"}->[$j]) {
		    # It should now be safe to add the $j:th annotation. It has been checked against all others.	    
		    $nr=$$annotation{nr};
		    
		    if($DEBUG) {
			$DEBUG && print "DEBUG: Merging since ";
			$olap && print "olap = $olap ";
			foreach my $olapmergeable (@mergeable_types) {
			    $e_list{"olap_$olapmergeable"}->[$j] && print "olap_$olapmergeable= ".$e_list{olap_$olapmergeable}->[$j]." ";
			}
			print "\n";
		    }
		    $DEBUG && print "DEBUG: Merging e_list nr $j ($e_list{id}->[$j]) as new annotation nr $nr.\n";
		    foreach my $field (keys %$annotation) {
			$$annotation{$field}->[$nr]=$e_list{$field}->[$j];
		    }
		    
		    $$annotation{uid}->[$nr] = $$annotation{unr};
		    $$annotation{unr}++;
			
		    $$annotation{id}->[$nr]="mb_$nr";
		    $$annotation{type}->[$nr]="merge";
		    $$annotation{color}->[$nr]="purple";
		    if($e_list{frame}->[$j]<3) {
			$$annotation{level}->[$nr]=$$annotation{frame}->[$nr];
		    } elsif ($e_list{frame}->[$j]>2) {
			$$annotation{level}->[$nr]=$$annotation{frame}->[$nr]-5;		   
		    }
		    
		    $$annotation{comment}->[$nr] = "Merged with olap from ";
		    foreach my $allowed_type (@mergeable_types) {
			if($e_list{"olap_$allowed_type"}->[$j]) {
			    $$annotation{comment}->[$nr] .= "$allowed_type ";
			}
		    }
		    $$annotation{comment}->[$nr] .= ".";
		    
		    $$annotation{nr}++;
		    if($$merge_status{clone_name}) {
			annotation_name($annotation,$nr,$seqName,"clone");
		    }
		}
	    }
	}
    }
}

return 1;
