#!/usr/bin/perl -w
# 
# IO routines for semi-automatic genomic feature annotation GUI
# Daniel Nilsson, 010619
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

package GuiIO;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read_prefs 
	     save_or_discard 
	     open_sequence select_sequence open_sequence_fasta_file 
	     last_opened_files add_last_opened_file
	     getFastaSeq getFastaQual 
	     open_genbank import_genbank
	     annotation_est_serialize annotation_estorf_serialize annotation_glimmer_serialize annotation_manual_serialize annotation_orf_serialize annotation_polypy_serialize annotation_serialize annotation_merge_serialize annotation_testcode_serialize annotation_genbankcds_serialize annotation_st_serialize sheet_serialize plotstate_serialize 
	     file_close file_load file_save 
	     export_sequin export_gff export_genbank 
	     qblast qblastcd check_waiting_time parse_blastn_reply parse_blastx_reply parse_blastp_reply add_blasthit 
	     run_external 
	     max min maxref minref log2 sum 
	     testcode 
	     dinucleotide_descriptor trinucleotide_descriptor 
	     nt2aa 
	     logo);

use Tk;
#use Tk::Dialog;
#use Tk::FileSelect;
use Bio::Tools::GFF;

use GuiAnnotator;

# use strict;

my $DEBUG = $main::DEBUG;
my $DEVEL = $main::DEVEL;
my $WARNING = $main::WARNING;

# file helper subs
sub read_prefs;

sub last_opened_files;
sub add_last_opened_file;

sub save_or_discard;
sub open_sequence;

use Bio::Seq;
use Bio::SeqIO;

sub open_genbank;
sub import_genbank;

sub file_load;
sub file_save;
sub file_close;

sub annotation_est_serialize;
sub annotation_estorf_serialize;
sub annotation_glimmer_serialize;
sub annotation_manual_serialize;
sub annotation_orf_serialize;

sub annotation_polypy_serialize;
sub annotation_st_serialize;
sub annotation_testcode_serialize;
sub annotation_serialize;

sub sheet_serialize;
sub plotstate_serialize;

# file i/o
sub getFastaSeq;
sub getFastaQual;
sub select_sequence;
sub open_sequence_fasta_file;

sub export_sequin;
sub export_gff;
sub export_genbank;

# qblast at ncbi
use Socket;

sub qblast;
sub check_waiting_time;
sub parse_blastn_reply;
sub parse_blastx_reply;
sub parse_blastp_reply;
sub blastresult_interact;

sub add_blasthit;

# unexported blasthit helpers
sub permanent_temporary;
sub remove_temporary;
my $qblast_norm_win; # package-wide var to ensure only one open "blasthits_normalise" window.

# run external program

sub run_external;

# math helpers
sub max;
sub min;

sub maxref;
sub minref;

sub log2;
sub sum;

# misc
sub nt2aa;
sub logo;

# file helper subs

sub read_prefs {
    my $prefs_file = shift;
    my $sheet = shift;

    if(! -e $prefs_file) {
	$DEBUG && print "DEBUG: no preferences file $prefs_file found. Continuing anyway.\n";
	return;
    }
    
    open PREFS, $prefs_file;
    while(<PREFS>) {
	chop;
	if (/^\s*\#/) {
	    # comment!
	    $DEBUG && print "DEBUG: prefs comment: $_\n";
	} elsif (/^\s*$/) {
	    # pretty blank line
	} else {
	    my ($key, $arg) = /\s*([\w_.\d]+?)\s+(.+)/;
	    $DEVEL && $DEBUG && print STDERR "DEVEL: key $key arg $arg found in prefs file $prefs_file.\n";
	    if($key ne "" && defined($arg) ) {
		$DEBUG && print STDERR "DEBUG: setting preferences, key $key arg $arg from prefs file $prefs_file.\n";
		$$sheet{$key} = $arg;
	    }
	}
    }
    close PREFS;
}

sub last_opened_files {
    my $sheet = shift;

    my @lastfiles;

    if (!defined ($$sheet{lastfilesfile}) or $$sheet{lastfilesfile} eq "") {
	$WARNING && print "WARNING: no lastfilesfile defined. Please specify one in your user aguirc-file, or disable the uselastfiles feature altogheter in the same file.\n";
	return;
    }
    if(! -e $$sheet{lastfilesfile}) {
	return;
    }

    my $error = 0;
    open LASTFILES, "<$$sheet{lastfilesfile}" or $error = 1;
    while(<LASTFILES>) {
	chop;
	push @lastfiles, $_;
    }
    close LASTFILES;
    return @lastfiles;
}

sub add_last_opened_file {
    my $sheet = shift;
    my $newfile = shift;
    
    my @last_files = last_opened_files($sheet);   
    
    if (!defined ($$sheet{lastfilesfile}) or $$sheet{lastfilesfile} eq "") {
	$WARNING && print "WARNING: no lastfilesfile defined. Please specify one in your user aguirc-file, or disable the uselastfiles feature altogheter in the same file.\n";
	return;
    }
    
    # checking for duplicate entries on lastfiles menu..
    my $dupe = 0;
    my $dupenr;
    my $i = 0;
    foreach(@last_files) {
	if($newfile eq $_) {
	    $dupe = 1;
	    $dupenr = $i;
	}
	$i++;
    }

    if($dupe) {
	$DEBUG && print "DEBUG: file already on lastfiles menu - moving  $newfile to to top of list.\n";
	splice @last_files, $dupenr, 1;
	unshift @last_files, $newfile;
    } else { 
	$DEBUG && print "DEBUG: adding $newfile to list of last opened/written files.\n";
	unshift @last_files, $newfile;
    }

    if(@last_files > $$sheet{nolastfiles}) {
	my $dumpme = pop @last_files;
	$DEBUG && print "DEBUG: removing $dumpme from list of last opened/written files.\n";
    }

    my $error = 0;
    open LASTFILES, ">$$sheet{lastfilesfile}" or $error = 1;
    if ($error == 0) {
	foreach(@last_files) {
	    print LASTFILES  "$_\n";
	}
	close LASTFILES;
    }
}

sub save_or_discard {
    my $main=shift;
    my $sheet=shift;
    my $annotation=shift;
    my $plotstate=shift;
    
    my($save, $discard, $cancel) = ('Save', 'Discard', 'Cancel');
    if (not Exists($save_or_discard_dialog)) {
	$save_or_discard_dialog = $$main->Dialog(
						-title          => 'Worksheet changed?',
						-text           => 'Worksheet has been changed.  Save changes to disk?',
						-bitmap         => 'info',
						-default_button => $save,
						-buttons        => [$save, $discard, $cancel],
						);
    }
    
    # local-lock & wait
    my $ok=$save_or_discard_dialog->Show;
    
    # Return value dependent on user choice.
    if($ok eq $save) {
	file_save($main,$sheet,$annotation,$plotstate);
    } elsif ($ok eq $discard) {
	return 1;
    } elsif ($ok eq $cancel) {
	return 0;
    }
}

sub open_sequence {
    my $main=shift;
    my $seqName=shift;		# catch references to seq data
    my $seq=shift;
    my $sheet=shift;

    my @filetypes=(["Fasta files",'*']);

    my $selectedFile= $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => "",-title=>"Open sequence file");

    if(!$selectedFile) {
	print "Open sequence cancelled.\n";
	return 0;
    } else {    
	print "Open sequence file $selectedFile.\n";      
	
	$WARNING && print "WARNING: Assuming fasta format...\n";
	
	# Open sequence fasta file... Returns true only when file selected.
	if(open_sequence_fasta_file($main,$selectedFile,$seqName,$seq)) {
	    $$sheet{status}="dirty"; # Sheet is now dirty.
	    $DEBUG && print"DEBUG: Open sequence fasta file returns true (sequence selected).\n";
	    return 1;
	} else {
	    $DEBUG && print"DEBUG: Open sequence fasta file returns false (sequence selection error/cancel).\n";
	    return 0;
	}
    }   
}

sub open_genbank {
    my $main=shift;

    my $annotation = shift;

    my $seqName=shift;		# catch references to seq data
    my $seq=shift;
    my $sheet=shift;

    my $selected_file = shift;

    my @filetypes=(["Genbank files",'*']);
    !defined($selected_file) && ($selected_file = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => "",-title=>"Open genbank file"));

    if(!$selected_file) {
	print "Open sequence cancelled.\n";
	return 0;
    } else {
	print "Open genbank file $selected_file.\n";   

	my $gbkin = Bio::SeqIO->new('-format' => 'GenBank', -file => $selected_file);

	if($gbkin == 0) {
	    $WARNING && print "WARNING: genbank seq open returns non-true.\n";
	    return 1;
	}

      main::add_last_opened_file($sheet, $selected_file);

	$DEBUG && print "DEBUG: Read genbank sequence..\n";
	my $gbkseq = $gbkin->next_seq();

	if($gbkseq) {
	    $DEBUG && print "DEBUG: Read genbank sequence ok ..\n"; 
	    $$seq = $gbkseq->seq;
	    $$seqName = $gbkseq->display_id();

	    import_genbank($main, $gbkseq, $annotation, $seqName, $seq, $sheet);
	    # unhide genbank cds level -- note that this counters explicit user instructions from prefs file...
	    $$sheet{display_genbankcds_level}=1;
	    # main::set_main_view_height($sheet);       

	    $$sheet{file} = $selected_file;
	    $$sheet{status}="dirty"; # Sheet is now dirty.

	    $DEBUG && print"DEBUG: Open genbank file returns true.\n";
	    return 1;
	} else {
	    $DEBUG && print"DEBUG: Open genbank file returns false.\n";
	    return 0;
	}
    }
}

sub import_genbank {
    my $main=shift;
    my $gbkseq = shift;
    my $annotation=shift;
    my $seqName = shift;
    my $seq = shift;
    my $sheet = shift;

    $DEBUG && print "DEBUG: Importing genbank.. \n";
    foreach my $feat ( $gbkseq->all_SeqFeatures() ) {
	if( ! defined($feat) ) {
	    # Skip feature-less entries, such as contig scaffolds..
	    next;
	}

	if($feat->primary_tag eq "CDS") {
	    # assign id & store annotation
	    my $nr;
	    
	    if($$annotation{nr}) {
		$nr=$$annotation{nr};
	    } else {
		$nr=0;
		$$annotation{nr}=0;
	    }
	    
	    $$annotation{uid}->[$nr] = $$annotation{unr};
	    $$annotation{unr}++;

      	    $$annotation{id}->[$nr]="ab_$nr";
	    $main::annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;
	    
	    $$annotation{seqName}->[$nr]=$$seqName; # Pretty redundant.. Remove?
	    $$annotation{type}->[$nr]='genbankcds';
	    $$annotation{start}->[$nr]=$feat->start;
	    $$annotation{stop}->[$nr]=$feat->end;

	    if($feat->strand==0 || $feat->strand==1) {
		$$annotation{frame}->[$nr] = (($feat->start % 3)?($feat->start % 3):3) - 1;
		$DEBUG && print "DEBUG: Feat start ",$feat->start,", mod 3 is ", $feat->start % 3, " and frame $$annotation{frame}->[$nr].\n";
	    } elsif ($feat->strand==-1) {
	        $$annotation{frame}->[$nr] = 5 - ((length($$seq) - ($feat->end)) % 3);
		$DEBUG && print "DEBUG: Feat end ",$feat->end,", seq len ", length($$seq)," diff ",length($$seq) - $feat->end ,"  mod 3 is ", $feat->end % 3," and frame $$annotation{frame}->[$nr].\n";
	    } else {
		$WARNING && print "WARNING: No strand indication for cds feature (uid $$annotation{uid}->[$nr], id $$annotation{id}->[$nr])\n";
	    }

	    if($$annotation{frame}->[$nr] < 3) {
		$$annotation{level}->[$nr] = $$annotation{frame}->[$nr];    
	    } else {
		$$annotation{level}->[$nr] = $$annotation{frame}->[$nr] - 5;
	    }

	    if($feat->has_tag('gene')) {
		my @gene = $feat->each_tag_value('gene');
		$$annotation{name}->[$nr] = $gene[0];
	    }
	    if($feat->has_tag('product')) {
		my @product = $feat->each_tag_value('product');
		$$annotation{comment}->[$nr] .= " ".$product[0];
		# print " ", $product[0];
	    }
	    if($feat->has_tag('protein_id')) {
		my @protein_id = $feat->each_tag_value('protein_id');
		$$annotation{comment}->[$nr] .= " ".$protein_id[0];
	    }
	    if($feat->has_tag('note')) {
		my @note = $feat->each_tag_value('note');
		$$annotation{note}->[$nr] .= join("\n", @note);
	    }	

	    $$annotation{color}->[$nr]='dark turquoise';
	    
	    $$annotation{nr}++;	# Adding an entry, so increase counter..

	    $$sheet{status}="dirty"; # Sheet is now dirty.
	} else {
	    # what gbk features to keep? all?
	    $DEVEL && print "DEVEL: Unknown feature with primary tag ", $feat->primary_tag, " at ",$feat->start,"-",$feat->end," (",$feat->strand,") encountered while parsing genbank entry.\n";
	}
    }
}

sub file_load {
    my $main=shift;
    my $canvas=shift;  
    my $annotation=shift;
    my $sheet=shift;
    my $plotstatus=shift;
    my $seqName=shift;
    my $seq=shift;
    my $selectedFile = shift;

    my $ok = file_close($main,$canvas,$sheet,$annotation,$seqName,$seq,$plotstatus);
    if($ok) {
	$WARNING && print "File saved ok, or discarded from file_load.\n";
    } else {
	$WARNING && print "Cancel new.\n";
	return 0;			# Ok, so cancel new.
    }

#  my $loadDialog = $$main->FileSelect(-directory=>'.',-filelabel=>'Load file',-filelistlabel=>'Files');
#  $loadDialog->title("A GUI - Load worksheet");
#  my $selectedFile = $loadDialog->Show;

    if(!defined($selectedFile)) {
	my @filetypes=(['GUI WorkSheet','.gws'],["All files",'*']);

	# $DEBUG && print "DEBUG: Trying FBox getOpenFile!\n";
	$selectedFile = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '.gws',-title=>"Load worksheet");
	
	if(!$selectedFile) {
	    $WARNING && print "Cancel load.\n";
	    return 0;
	} 
    }

    $WARNING && print "Load $selectedFile.\n";

    $$sheet{file}=$selectedFile;

    # delete old canvas items..
    $$canvas->delete('all');
    
    # Empty annotation hash, taking care to save the highest unique identifier..
    my $unr = $$annotation{unr};
    %$annotation=();
    $$annotation{unr}=$unr;
    
    $$sheet{selected} = '';

    $$seqName="";
    $$seq="";

    open LOAD, "<$selectedFile" || die "DEBUG ERROR: Could not open file $selectedFile for reading.\n";

    # Statevars  
    my $fileOk=0;
    my $inVersion=0;
    my $inSequence=0;
    my $inAnnotation=0;

    my $inNote=0;
    my $firstNoteLine=0;
    my $inComment=0;
    my ($noteline,$commentline);

    my $version="010611.1314";

    my $fileVersion="";
    my $nr=0;

    while(<LOAD>) {
	chop;

	if(!$fileOk) {
	    # Check file "magic".. 
	    if(m/^A\s{1}GUI/) {
		$fileOk=1;
	    } else {
		$WARNING && print "Oopps... Incorrect file format in $selectedFile.\nNot a worksheet file.\n";
		close LOAD;
		return 0;
	    }
	}

	if($inVersion) {
	    # version number of file..
	    if (m/^\};/) {
		$inVersion=0;
		if(!$fileVersion) {
		    $WARNING && print "File format error reading file version number in $selectedFile.\n";
		    close LOAD;
		    return 0;
		}
	    } else {
		($fileVersion)=m/\s{1}([\d\.]+)/;
		if($fileVersion ne $version) {
		    $WARNING && print "File version error: $selectedFile is version $fileVersion, and GUI I/O is version $version.\n";
		    if($fileVersion < $version) {
			$WARNING && print "It seems that $fileVersion is older than the current. Attempting to import anyway.. \n";
		    } else {
			$WARNING && print "It seems that $fileVersion is newer than the current. Please upgrade this software to be able to open that file!\n";
			close LOAD;
			return 0;
		    }
		}
	    }
	} elsif ($inSequence) {
	    if(m/^\};/) {
		$inSequence=0;
	    } else {
		if(!$$seqName) {
		    ($$seqName)=m/^\s{1}(.+)/;
		    $DEVEL && print "DEVEL: Load seq $$seqName.\n";
		} elsif (!$$seq) {
		    ($$seq)=m/^\s{1}(.+)/;
		}
	    }
	} elsif ($inAnnotation) {
	    if($inNote) {
		if(m/^\s{1}\}/) {
		    $inNote=0;
		} else {	    
		    if($firstNoteLine) {	
			# strip any initial ws on first note line
			($noteline)=m/\s*(.+)/;
			$firstNoteLine = 0;
		    } else {
			$noteline=$_;
		    }

		    $$annotation{note}->[$nr].="$noteline\n";
		}
	    } elsif ($inComment) {
		if(m/^\s{1}\}/) {
		    $inComment=0;
		} else {   
		    # strip any initial ws.. 
		    ($commentline)=m/\s*(.+)/;
		    if($$annotation{comment}->[$nr]) {
			$$annotation{comment}->[$nr].="\n$commentline";
		    } else {
			# Last line of comment needs no CR, and an initial 
			# CR on first line is no good either.
			$$annotation{comment}->[$nr].="$commentline";
		    }
		}
	    } elsif (m/^\};/) {
		# End of annotation. Set some non-saved data..

		if(! defined ($$annotation{note}->[$nr]) ) {
		    $$annotation{note}->[$nr] = "";
		}	 
		if(! defined ($$annotation{name}->[$nr]) ) {
		    $$annotation{name}->[$nr] = "";
		}   
		if(! defined ($$annotation{comment}->[$nr]) ) {
		    $$annotation{comment}->[$nr] = "";
		}

		$$annotation{uid}->[$nr] = $$annotation{unr};
		$$annotation{unr}++;

		$$annotation{id}->[$nr]="ab_$nr";
		$main::annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;

		$$annotation{seqName}->[$nr]=$$seqName;

		if($$annotation{type}->[$nr] eq 'EST') { 
		    my ($estName) = $$annotation{comment}->[$nr]=~m/(.+)/;
		    $DEVEL && print "DEVEL: Got estname $estName from comment $$annotation{comment}->[$nr].\n";
		    
		    if($estName =~ /5{1}\'{1}/ ) {
			$$annotation{frame}->[$nr] = 0;	 
		    } elsif ($estName =~ /3{1}\'{1}/ ) {
			$$annotation{frame}->[$nr] = 3;
		    } else {
			$$annotation{frame}->[$nr] = $AnnotationBoxFrame; # Box
			$DEVEL && print "DEVEL: Could not determine orientation for $estName. Boxing it!\n";
		    }
		}

		# check frames, and modify to modulobased framing if need be
		if(defined($$annotation{frame}->[$nr])) {
		    my $frame;
		    if($$annotation{frame}->[$nr] < 3) {
			# update frame -- forward pass
			$frame = ($$annotation{start}->[$nr]%3)?($$annotation{start}->[$nr]%3)-1:2;
		    } else {
			$frame = 5-((length($$seq)-$$annotation{stop}->[$nr])%3);
		    }
		    
		    if($frame != $$annotation{frame}->[$nr]) {
			$WARNING && print "WARNING: incorrect frame ($$annotation{frame}->[$nr]) found when loading annotation nr $nr type $$annotation{type}->[$nr]. Setting frame $frame.\n";
			$$annotation{frame}->[$nr] = $frame;
		    }
		} # no frame is dealt with in the following...
		
		# re-locate the annotation levels upon loading?	
		if(defined($$annotation{frame}->[$nr]) && ($$annotation{frame}->[$nr] < 3)) {
		    if ($annotatorHeight{$$annotation{type}->[$nr]} == 3) {
			$$annotation{level}->[$nr] = $$annotation{frame}->[$nr];
		    } elsif ($annotatorHeight{$$annotation{type}->[$nr]} == 1) {
			$$annotation{level}->[$nr] = 0;
		    } else {
			$$annotation{level}->[$nr] = 0; # awaiting level_layout
		    }
		} elsif (defined($$annotation{frame}->[$nr]) && ($$annotation{frame}->[$nr] < 6)) {
		    if ($annotatorHeight{$$annotation{type}->[$nr]} == 3) {
			$$annotation{level}->[$nr] = $$annotation{frame}->[$nr] - 5;
		    } elsif ($annotatorHeight{$$annotation{type}->[$nr]} == 1) {
			$$annotation{level}->[$nr] = 0;
		    } else {
			$$annotation{level}->[$nr] = 0; # awaiting level_layout
		    }
		} elsif (defined($$annotation{frame}->[$nr]) && ($$annotation{frame}->[$nr] == $AnnotationBoxFrame)) {
		    $$annotation{level}->[$nr] = 0;
		} else {
		    $WARNING && print "WARNING: no frame information found when loading annotation nr $nr of type ".($$annotation{type}->[$nr])." from file. Using unoriented box on display relative level 0.\n";
		    $$annotation{frame}->[$nr] = $AnnotationBoxFrame;
		    $$annotation{level}->[$nr] = 0;
		}

		$DEVEL && print "DEVEL: Load to level ".($$annotation{level}->[$nr]).", since frame ".($$annotation{frame}->[$nr]).".\n";
		
		# Ready for next annotation..
		$nr++;
		$inAnnotation=0;
	    } elsif (m/^\s{1}type\s{1}\{$/) {
		$_=<LOAD>; 
		chop;
		($$annotation{type}->[$nr])=m/^\s{2}(.+)/;
		if($$annotation{type}->[$nr] =~ m/^gff_/) {
		    my $height = 3;
		    my $bioframed = 0;
		    add_gff_annotator($sheet, $$annotation{type}->[$nr], $height, $bioframed); #stacked or framed, bro.. Need to put that in the ff i guess - we could use the obsolete level field for that info to avoid incompatiblity..
		}
	    } elsif (m/^\s{1}start\s{1}\{$/) {
		$_=<LOAD>; 
		chop; 
		($$annotation{start}->[$nr])=m/^\s{2}(\d+)/;
	    } elsif (m/^\s{1}stop\s{1}\{$/) {
		$_=<LOAD>; 
		chop; 
		($$annotation{stop}->[$nr])=m/^\s{2}(\d+)/;
	    } elsif (m/^\s{1}frame\s{1}\{$/) {
		$_=<LOAD>; 
		chop; 
		($$annotation{frame}->[$nr])=m/^\s{2}(\d+)/;
	    } elsif (m/^\s{1}level\s{1}\{$/) {
		$_=<LOAD>; 
		chop;
		# Level is re-calculated upon load.
	    } elsif (m/^\s{1}color\s{1}\{$/) {
		$_=<LOAD>; 
		chop; 
		($$annotation{color}->[$nr])=m/^\s{2}(.+)/;
	    } elsif (m/^\s{1}name\s{1}\{$/) {
		$_=<LOAD>;
		chop; 
		($$annotation{name}->[$nr])=m/^\s{2}(.+)/;
	    } elsif (m/^\s{1}comment\s{1}\{$/) {
		$$annotation{comment}->[$nr]="";
		# May be multi-line. Push a state..
		$inComment=1;
	    } elsif (m/^\s{1}note\s{1}\{$/) {
		$$annotation{note}->[$nr]="";
		# May be multi-line. Push a state..
		$firstNoteLine=1;
		$inNote=1;	
	    } elsif (m/^\s{1}\}/) {
		# End of type/start/stop/frame/level/color
	    }
	} elsif (m/^version\s{1}\{$/) {
	    $DEBUG && print "DEBUG: Loading version nr..\n";
	    $inVersion=1;
	} elsif (m/^sequence\s{1}\{$/) {
	    $DEBUG &&  print "DEBUG: Loading sequence...\n";
	    $inSequence=1;
	} elsif (m/^annotation\s{1}\{$/) {
	    $DEVEL && print "DEVEL: Loading annotation $nr..\n";
	    $inAnnotation=1;
	} elsif (m/^annotations\s{1}\{/) {
	    $_=<LOAD>; 
	    chop;
	    my ($annotations)=m/^\s{1}(\d+)/;
	    if($nr != $annotations) {
		$WARNING && print "Annotation count check failed with nr in mem = $nr and nr in file $annotations\n";
		close LOAD;
		return 0;
	    } else {
		$DEBUG && print "DEBUG: Annotation count check pass with $nr annotations read.\n"
		}
	    $_=<LOAD>; # "\};\n";
	}
    }
    $$annotation{nr}=$nr;

    close LOAD;

  main::add_last_opened_file($sheet, $$sheet{file});

  main::level_layout($annotation,'EST'); # Would it be nice to re-locate the annotation levels upon loading?
  main::level_layout($annotation,'blast');
  main::level_layout($annotation,'regexp');
  main::level_layout_gff($annotation);

    $$sheet{status} = "clean";
    
    return 1;
}

sub file_save {
    my $main=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $plotstate=shift;
    my $name_change=shift;

#  my $seqName=shift; # include 'em in sheet?
#  my $seq=shift;

    my $version="010611.1314";

    my $selectedFile = undef;

    if($name_change or !defined($$sheet{file}) or $$sheet{file} eq '') {
#    my $saveDialog = $$main->FileSelect(-directory=>'.',-filelabel=>'Save file');
#    $saveDialog->title("A GUI - Save worksheet");
#    my $selectedFile=$saveDialog->Show;
	my @filetypes=(['GUI WorkSheet','.gws'],["All files",'*']);

	# $DEBUG && print "DEBUG: Trying FBox getSaveFile!\n";
	$selectedFile= $$main->getSaveFile(-filetypes=>\@filetypes,-initialfile => $main::seqName.".gws",-defaultextension => '.gws'); # using GLOBAL seqName

	if(!$selectedFile) {
	    $WARNING && print "WARNING: Cancel save.\n";
	    return;
	} else {
	    $WARNING && print "WARNING: Saving to file $selectedFile.\n";
	    $$sheet{file}=$selectedFile;
	}
    }

    open SAVE, ">$$sheet{file}" || die "DEBUG ERROR Could not open save-file for writing: $!\n";
    
    print SAVE "A GUI\n";
    print SAVE "version\t{\n\t$version\n};\n"; # Version
    print SAVE "sequence\t{\n\t".$main::seqName."\n\t".$main::seq."\n};\n"; # Sequence # using GLOBAL seq, seqName
    
    my $nr_of_st=0; # Keep until we decide to save Start/stops, if ever.
    my $nr;
    for($nr=0;$nr<$$annotation{nr};$nr++) {
	# annotation->serialize() for the different annotation types?
	# Would have to bless them then... Would be like real fun, but supposedly something for
	# the Next Version... (If that is in perl at all, and not in c++/Qt ie.)
	# I'll go for a poor-mans version of object->serialize(): object_serialize().. =)
	if($$annotation{type}->[$nr] eq "glimmer2") {
	    print SAVE annotation_glimmer_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "EST") {
	    print SAVE annotation_est_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "ESTORF") {
	    print SAVE annotation_estorf_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "ORF") {
	    print SAVE annotation_orf_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "ST") {
	    print SAVE annotation_st_serialize($annotation,$nr);
	    $nr_of_st++; # Keep until we decide to save Start/stops, if ever.       
	} elsif($$annotation{type}->[$nr] eq "polypy") {
	    print SAVE annotation_polypy_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "manual") {
	    print SAVE annotation_manual_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "merge") {
	    print SAVE annotation_merge_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "testcode") {
	    print SAVE annotation_testcode_serialize($annotation,$nr);
	} elsif($$annotation{type}->[$nr] eq "genbankcds") {
	    print SAVE annotation_genbankcds_serialize($annotation,$nr);
	} else {
	    print SAVE annotation_serialize($annotation,$nr); # generic	(gff*, regexp, blast)
	}
    }

    print SAVE "annotations\t{\n\t",$nr-$nr_of_st,"\n};"; # Keep until we decide to save Start/stops, if ever.
    print SAVE sheet_serialize($sheet);
    print SAVE plotstate_serialize($plotstate);
    
    close SAVE;

    print "File $$sheet{file} should be saved by now...\n";

  main::add_last_opened_file($sheet, $$sheet{file});

    $$sheet{status}='clean';
}

sub file_close {
    my $main=shift;
    my $canvas=shift;
    my $sheet=shift;
    my $annotation = shift;
    my $seqName=shift;
    my $seq=shift;
    my $plotstate = shift;

    # save-if-dirty-dialog, then clear...
    if($$sheet{status} eq "dirty") {
	my $ok = save_or_discard($main,$annotation,$sheet,$plotstate);
	if($ok) {
	    print "File saved ok, or discarded.\n";
	} else {
	    print "Cancel close.\n";
	    return 0;		# Ok, so cancel close.
	}
    }

    # Ok, deleting everything...  
    $$canvas->configure(-height=>$$sheet{canvas_seq_height});
    $$canvas->delete('all');
    
    $$sheet{selected}='';

    my $unr = $$annotation{unr};
    %$annotation=();
    $$annotation{unr}=$unr;

    $$seqName='';
    $$seq='';
    
    $$sheet{file}='';
    $$sheet{status}='closed';
}

sub annotation_est_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";
    
    # ids are not saved - new id's are issued upon load instead.
    # seqName is considered redundant so far. Issued upon load.
    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load??? heavy, but maby nice) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_estorf_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";
    
    # ids are not saved - new id's are issued upon load instead.
    # seqName is considered redundant so far. Issued upon load.
    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?)
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?)

    $serialize.="};\n";

    return $serialize;
}

sub annotation_polypy_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";
    
    # ids are not saved - new id's are issued upon load instead.
    # seqName is considered redundant so far. Issued upon load.
    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type  
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_orf_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";
    
    # ids are not saved - new id's are issued upon load instead.
    # seqName is considered redundant so far. Issued upon load.
    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_genbankcds_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";
    
    # ids are not saved - new id's are issued upon load instead.
    # seqName is considered redundant so far. Issued upon load.
    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_glimmer_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";

    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_testcode_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";

    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";

    return $serialize;
}

sub annotation_st_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";
    print "Start and stop annotations are not saved.\n";
    return $serialize;
}

sub annotation_manual_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";

    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";
    return $serialize;
}

sub annotation_merge_serialize {
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";

    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";
    return $serialize;
}

sub annotation_serialize {
    # generic serialize -- blast, regexp
    my $annotation=shift;
    my $nr=shift;
    
    my $serialize="";

    $serialize.="annotation\t{\n";

    $serialize.="\ttype\t{\n\t\t$$annotation{type}->[$nr]\n\t}\n"; # type
    $serialize.="\tstart\t{\n\t\t$$annotation{start}->[$nr]\n\t}\n"; # start
    $serialize.="\tstop\t{\n\t\t$$annotation{stop}->[$nr]\n\t}\n"; # start
    $serialize.="\tframe\t{\n\t\t$$annotation{frame}->[$nr]\n\t}\n"; # frame
    if($$annotation{comment}->[$nr]) {
	$serialize.="\tcomment\t{\n\t\t$$annotation{comment}->[$nr]\n\t}\n"; # comment
    }
    if($$annotation{note}->[$nr]) {
	$serialize.="\tnote\t{\n\t\t$$annotation{note}->[$nr]\n\t}\n"; # note
    }
    if($$annotation{name}->[$nr]) {
	$serialize.="\tname\t{\n\t\t$$annotation{name}->[$nr]\n\t}\n"; # name
    }
    $serialize.="\tlevel\t{\n\t\t$$annotation{level}->[$nr]\n\t}\n"; # level (issue on load?) 
    $serialize.="\tcolor\t{\n\t\t$$annotation{color}->[$nr]\n\t}\n"; # color (issue on load?) 

    $serialize.="};\n";
    return $serialize;
}

sub sheet_serialize {
    my $sheet=shift;
    my $serialize="";

    print "Sheet status not saved. No need yet..\n";  
    return $serialize;
}

sub plotstate_serialize {
    my $plotstate=shift;
    my $serialize="";

    print "Plot status not saved. No need yet..\n";  
    return $serialize;
}

# File I/O

sub export_genbank {

    my $main=shift;
    my $seqName=shift; # ref 
    my $seq=shift;
    my $annotation=shift; # ref
    my $annotator_selection = shift; # ref

    # Let the user choose a file for export, and open it
    my @filetypes=(['General Feature Format 2 file','.gbk'],["All files",'*']);
    
    my $selectedFile = $$main->getSaveFile(-filetypes=>\@filetypes,-initialfile => "$$seqName.gbk",-defaultextension => '.gbk');
    
    if(!$selectedFile) {
	print "Cancel export.\n";
	return;
    } else {
	print "Export $selectedFile.\n";
    }

    my $gbkout = Bio::SeqIO->new('-file' => ">$selectedFile", '-format' => 'GenBank');

    my $gbkseq = Bio::Seq->new(-seq => $$seq, -display_id => $$seqName, -id => $$seqName, -moltype => 'dna', -accession_number => $$seqName, -primary_id  => $$seqName, -desc => $$seqName);
#      my Bio::Seq, write to file

    ANNOT: for(my $i=0;$i<$$annotation{nr};$i++) {

	if(!$$annotator_selection{$$annotation{type}->[$i]}) { # Check type of annotation against export selection
	    # and, seriously, as cdregion only protein types of annotation should go... we leave that to the user at the moment
	    next ANNOT;
	}

	# create a temporary feature to write
	my $strand = ($$annotation{frame}->[$i] < 3 || $$annotation{frame}->[$i] > 5) ? 1 : -1;
	my $primary_tag = 'CDS';
	my $source_tag = $$annotation{type}->[$i];
	# the seqname of the exported feature is default - "SEQ" -- annoying?
	my $feature = Bio::SeqFeature::Generic->new(-start=> $$annotation{start}->[$i],
						    -end =>  $$annotation{stop}->[$i],
						    -strand => $strand,
						    -primary => $primary_tag,
						    -source => $source_tag,
						    -tag => {
							comment => $$annotation{comment}->[$i],
							note => $$annotation{note}->[$i] }
						    );
	$gbkseq->add_SeqFeature($feature);
	
	# write out
	# $gbkout->write_feature($feature);
    }

    $gbkout->write_seq($gbkseq);
        
    $gbkout->close();

    return $selectedFile;
}

sub export_gff {
    my $main=shift;
    my $seqName=shift; # ref 
    my $seq=shift;
    my $annotation=shift; # ref
    my $annotator_selection = shift; # ref

    # Let the user choose a file for export, and open it
    my @filetypes=(['General Feature Format 2 file','.gff'],["All files",'*']);
    
    my $selectedFile = $$main->getSaveFile(-filetypes=>\@filetypes,-initialfile => "$$seqName.gff",-defaultextension => '.gff');
    
    if(!$selectedFile) {
	print "Cancel export.\n";
	return;
    } else {
	print "Export $selectedFile.\n";
    }

    my $gffout = Bio::Tools::GFF->new(-file => ">$selectedFile", -gff_version => 2);

    ANNOT: for(my $i=0;$i<$$annotation{nr};$i++) {

	if(!$$annotator_selection{$$annotation{type}->[$i]}) { # Check type of annotation against export selection
	    # and, seriously, as cdregion only protein types of annotation should go... we leave that to the user at the moment
	    next ANNOT; 
	}

	# create a temporary feature to write
	my $strand = ($$annotation{frame}->[$i] < 3 || $$annotation{frame}->[$i] > 5) ? 1 : -1;
	my $primary_tag = "CDS"; #$$annotation{type}->[$i];
	my $source_tag = $$annotation{type}->[$i];
	# the seqname of the exported feature is default - "SEQ" -- annoying?
	my $feature = Bio::SeqFeature::Generic->new(-start=> $$annotation{start}->[$i],
						    -end =>  $$annotation{stop}->[$i],
						    -strand => $strand,
						    -primary => $primary_tag,
						    -source => $source_tag,
						    -tag => {
							comment => $$annotation{comment}->[$i],
							note => $$annotation{note}->[$i] }
						    );
	# attach?

	# write out
	$gffout->write_feature($feature);
    }
    
    $gffout->close();
}

sub export_sequin {
    my $main=shift;
    my $seqName=shift; # ref 
    my $seq=shift;
    my $annotation=shift; # ref
    my $annotator_selection = shift; # ref

    my $seq_length=length($$seq);

    # Let the user choose a file for export, and open it
    my @filetypes=(['Sequin ASN.1 file','.sqn'],["All files",'*']);
    
    my $selectedFile = $$main->getSaveFile(-filetypes=>\@filetypes,-initialfile => "$$seqName.sqn",-defaultextension => '.sqn');
    
    if(!$selectedFile) {
	print "Cancel export.\n";
	return;
    } else {
	print "Export $selectedFile.\n";
    }

    open SEQUIN, ">$selectedFile";

    # Write a sequin (part of ASN.1) style Seq-entry. Here goes.
    
    # Header?
    print SEQUIN "Seq-entry ::= set {\n";
    print SEQUIN "\tclass nuc-prot ,\n"; # Ignoring some OPTIONALs
    print SEQUIN "\tseq-set {\n";

    my $line_length=80;
    
    # Beginning sequence-set by writing the actual sequence. 
    print SEQUIN "\t\tseq {\n";
    print SEQUIN "\t\t\tid {\n";
    print SEQUIN "\t\t\t\tlocal\n";
    print SEQUIN "\t\t\t\t\tstr \"$$seqName\" } ,\n";
    print SEQUIN "\t\t\tdescr {\n";
    print SEQUIN "\t\t\t\ttitle \"$$seqName\" ,\n";
    print SEQUIN "\t\t\t\tmolinfo {\n";
    print SEQUIN "\t\t\t\t\tbiomol genomic } } ,\n"; # skipping date
    print SEQUIN "\t\t\tinst {\n";
    print SEQUIN "\t\t\t\trepr raw ,\n";
    print SEQUIN "\t\t\t\tmol dna ,\n";
    print SEQUIN "\t\t\t\tlength $seq_length ,\n";
    print SEQUIN "\t\t\t\tseq-data \n";
    my $iupacseq=$$seq;
    $iupacseq=~tr/atgcnxX/ATGCNNN/; # Make sure sequence format is ok. 
    my $iupacseq_entry="\t\t\t\t\tiupacna \"$iupacseq\" } ,";

    for (my $j=0;$j<POSIX::ceil(length($iupacseq_entry)/$line_length);$j++) {
	print SEQUIN substr($iupacseq_entry,$j*$line_length,$line_length) ,"\n";
    }

    # First add annotations as annot to this sequence
    print SEQUIN "\t\t\tannot {\n";
    print SEQUIN "\t\t\t\t{\n"; # SET-OF Seq-annot
    print SEQUIN "\t\t\t\t\tdata\n"; # ignore some OPTIONALs
    print SEQUIN "\t\t\t\t\t\tftable {\n";
    
    my $firstloop=1;

  ANNOT: for(my $i=0;$i<$$annotation{nr};$i++) {
      # DEBUG!!!
      if(!$$annotator_selection{$$annotation{type}->[$i]}) { # Check type of annotation against export selection
	  # and, seriously, as cdregion only protein types of annotation should go... we leave that to the user at the moment
	  next ANNOT; 
      }

      if($firstloop == 1) {
	  $firstloop=0;
      } else {
	  print SEQUIN "\t\t\t\t\t\t\t\texp-ev not-experimental } ,\n";
      }

      print SEQUIN "\t\t\t\t\t\t\t{\n"; # in a SET-OF Seq-feat
      print SEQUIN "\t\t\t\t\t\t\t\tdata \n";
      print SEQUIN "\t\t\t\t\t\t\t\t\tcdregion {\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tcode {\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tid 1 } } ,\n"; # T. cruzi uses standard code,
      print SEQUIN "\t\t\t\t\t\t\t\tproduct\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\twhole\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tlocal\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\t\tstr \"$$annotation{id}->[$i]\" ,\n"; # IMPROVE: Change to name
      print SEQUIN "\t\t\t\t\t\t\t\tlocation\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\tint {\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tfrom ", $$annotation{start}->[$i]-1, " ,\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tto ", $$annotation{stop}->[$i]-1, " ,\n";
      if($$annotation{frame}->[$i] < 3) {
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\tstrand plus ,\n";
      } else {
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\tstrand minus ,\n";
      }
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tid\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\t\tlocal\n";
      
      print SEQUIN "\t\t\t\t\t\t\t\t\t\t\t\tstr \"$$seqName\" } ,\n"; # location ON genomic strand!
  }
    print SEQUIN "\t\t\t\t\t\t\t\texp-ev not-experimental } } } } } ,\n";
    
    $firstloop=1;
    my $last_nr;
    
  SEQ: for(my $i=0;$i<$$annotation{nr};$i++) {
      # DEBUG!!!
      if($$annotation{type}->[$i] ne "merge") { # Check type of annotation
	  next SEQ;
      }

      if($firstloop == 1) {
	  $firstloop=0;
      } elsif ($firstloop == 0) {
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\t\t\tstr \"$$annotation{id}->[$last_nr]\" } } } } } } ,\n"; 
      }

      print SEQUIN "\t\tseq {\n";
      print SEQUIN "\t\t\tid {\n";
      print SEQUIN "\t\t\t\tlocal\n";
      print SEQUIN "\t\t\t\t\tstr \"$$annotation{id}->[$i]\" } ,\n";
      print SEQUIN "\t\t\tdescr {\n";
      print SEQUIN "\t\t\t\tmolinfo {\n";
      print SEQUIN "\t\t\t\t\tbiomol peptide ,\n";
      print SEQUIN "\t\t\t\t\ttech concept-trans-a } } ,\n";
      print SEQUIN "\t\t\tinst {\n";
      print SEQUIN "\t\t\t\trepr raw ,\n";
      print SEQUIN "\t\t\t\tmol aa ,\n";
      
      my $annotated_sequence=substr($$seq,$$annotation{start}->[$i]-1,$$annotation{stop}->[$i]-$$annotation{start}->[$i]+1);

      my $annotated_protein;
      if($$annotation{frame}->[$i]<3) {
	  # Forward reading frame
	  $annotated_protein=nt2aa($annotated_sequence); # We could remove M and + if we'd like..
      } else {
	  # Reverse reading frame
	  $_=reverse(split(/ */,$annotated_sequence));
	  tr/atgcATGC/tacgTACG/;
	  $annotated_protein=nt2aa($_); # We could remove M and + if we'd like.. 
      }
      my $stop_char = chop($annotated_protein);
      if ($stop_char ne "+") {
	  $WARNING && print "WARNING: Tried to chop a + away from $$annotation{id}->[$i], but got a $stop_char instead! This should only happen for ORFs overlapping contig ends. Reinserting $stop_char at end of protein before search...\n";
	  $annotated_protein .= $stopchar;
      }
      my $annotated_protein_length=length($annotated_protein);
      
      print SEQUIN "\t\t\t\tlength $annotated_protein_length ,\n";
      print SEQUIN "\t\t\t\tseq-data\n";

      my $annotated_protein_entry="\t\t\t\t\tncbieaa \"$annotated_protein\" } ,";

      for (my $j=0;$j<POSIX::ceil(length($annotated_protein_entry)/$line_length);$j++) {
	  print SEQUIN substr($annotated_protein_entry,$j*$line_length,$line_length) ,"\n";
      }

      print SEQUIN "\t\t\tannot {\n";
      print SEQUIN "\t\t\t\t{\n";
      print SEQUIN "\t\t\t\t\tdata\n";
      print SEQUIN "\t\t\t\t\t\tftable {\n";
      print SEQUIN "\t\t\t\t\t\t\t{\n";
      print SEQUIN "\t\t\t\t\t\t\t\tdata\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\tprot {\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tname {\n";
      if(defined($$annotation{name}->[$i]) and $$annotation{name}->[$i] ne "") {
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\t\t\"$$annotation{name}->[$i]\" } ,\n";
      } else { 
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\t\t\"unknown\" } ,\n";
      }
      if($$annotation{note}->[$i]) {
	  # experimenting to see if we could get something like what Esteban wants.
	  # this is an irregular operation, so we can afford to spend quite a few cycles..
	  my $revised_note = "";
	  my $first_hit=1;

	  my $in_note=1;

	  my ($merged_from, $this_note);
	  my @these_notelines = split(/\n/,$$annotation{note}->[$i]);
	  
	  while($_ = shift(@these_notelines)) {
	      chop;
	      if(/\*\*\* NOTE FROM/) {		
		  ($merged_from) = /\*\*\* NOTE FROM\s*(\w+)\s*\w*\s*\*/;
		  $DEBUG && print "DEBUG: note from $merged_from (extracted from $_) for annotation $i.\n";
		  if($first_hit == 1) {
		      $revised_note .= "Hit by $annotatorName{$merged_from} (";
		      $first_hit = 0;
		  } else {
		      $revised_note .= " and also hit by $annotatorName{$merged_from} (";
		  }
		  #    $this_note
	      } elsif ($in_note) { # duh?
		  $revised_note .= $_;
	      }
	  }
	  # print SEQUIN "\t\t\t\t\t\t\t\t\t\tdesc \"$$annotation{note}->[$i]\" } ,\n";
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\tdesc \"$revised_note\" } ,\n";
      } else {
	  print SEQUIN "\t\t\t\t\t\t\t\t\t\tdesc \"unknown function\" } ,\n";
      }
      print SEQUIN "\t\t\t\t\t\t\t\tlocation\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\tint {\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tfrom 0 ,\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tto ",$annotated_protein_length-1," ,\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tstrand plus ,\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\tid\n";
      print SEQUIN "\t\t\t\t\t\t\t\t\t\t\tlocal\n";
      $last_nr=$i;
  }

    if(defined($last_nr) && ($last_nr > 0) ) {
	print SEQUIN "\t\t\t\t\t\t\t\t\t\t\t\tstr \"$$annotation{id}->[$last_nr]\" } } } } } } } }\n"; 
    } else {
	print SEQUIN "\t\t} }\n";
    }
    close SEQUIN;
}

sub testcode {
    # Testcode algortihm by JW Fickett, 1982.
    # I'm grateful to Gautam Aggarwal (SBRI, Seattle, USA) who provided me with his C implementation, 
    # to which this one bears much likeness.
    
    my $seq = shift;
    
    my $seqlen = length($seq);
    
    my @seq = split(/ */, lc $seq);

    my %fnuc;
    my %f;

    my $phase = 0;
    for (my $pos = 0; $pos < $seqlen; $pos++) {
	$phase = $pos % 3;
	$f{$seq[$pos]}[$phase]++;
    }
    
    foreach my $nuc ( 'a', 't', 'g', 'c' ) {
	$fnuc{$nuc} = $f{$nuc}[0] + $f{$nuc}[1] + $f{$nuc}[2];
    }
    
    $DEVEL && $DEBUG && print "The a counts are ", join(" ",@{$f{'a'}}), " t ",join(" ",@{$f{'t'}}), " g ", join(" ",@{$f{'g'}}), " c ", join(" ",@{$f{'c'}}),".\n";
    $DEVEL && $DEBUG && print "Max and min of a are ".maxref($f{'a'})." and ".minref($f{'a'}).".\n";
    $DEVEL && $DEBUG && print "The total frequency of a is ".$fnuc{'a'}." and the sequence (ORF) is $seqlen bp long.\n"; 
    
    my $posia = maxref($f{'a'})/(minref($f{'a'})+1);
    my $posit = maxref($f{'t'})/(minref($f{'t'})+1);
    my $posig = maxref($f{'g'})/(minref($f{'g'})+1);
    my $posic = maxref($f{'c'})/(minref($f{'c'})+1);

    my $conta = $fnuc{'a'} / $seqlen ;
    my $contt = $fnuc{'t'} / $seqlen ;
    my $contg = $fnuc{'g'} / $seqlen ;
    my $contc = $fnuc{'c'} / $seqlen ;

    my $Pmj = 0;
    
    my $tind;
    
    if($posia >= 0 && $posia < 1.1) {
	$tind = .22;
    } elsif($posia >= 1.1  && $posia < 1.2) {
	$tind = .20;
    } elsif($posia >= 1.2 && $posia < 1.3) {
	$tind = .34;
    } elsif($posia >= 1.3 && $posia < 1.4) {
	$tind = .45;
    } elsif($posia >=1.4 && $posia < 1.5) {
	$tind = .68;
    } elsif($posia >=1.5 && $posia < 1.6) {
	$tind = .58;
    } elsif($posia >=1.6 && $posia < 1.7) {
	$tind = .93;
    } elsif($posia >= 1.7 && $posia < 1.8) {
	$tind = .84;
    } elsif($posia >= 1.8 && $posia < 1.9) {
	$tind = .68;
    } elsif($posia >=1.9) {
	$tind=.94;
    }
    $Pmj += .26*$tind;

    if($conta >= 0 && $conta < .17) {
	$tind=.21;
    } elsif($conta >= .17  && $conta < .19) {
	$tind=.81;
    } elsif($conta >= .19 && $conta < .21) {
	$tind= .65;
    } elsif($conta >= .21 && $conta < .23) {
	$tind = .67;
    } elsif($conta >= .23 && $conta < .25) {
	$tind= .49;
    } elsif($conta >= .25 && $conta < .27) {
	$tind = .62;
    } elsif($conta >= .27 && $conta < .29) {
	$tind = .55;
    } elsif($conta >= .29 && $conta < .31) {
	$tind = .44;
    } elsif($conta >= .31 && $conta < .33) {
	$tind = .49;
    } elsif($conta >= .33 && $conta < .99) {
	$tind = .28;
    }
    $Pmj += .11*$tind;

    if($posit >= 0 && $posit < 1.1) {
	$tind=.09;
    } elsif($posit >= 1.1  && $posit < 1.2) {
	$tind = .09;
    } elsif($posit >= 1.2 && $posit < 1.3) {
	$tind = .20;
    } elsif($posit >= 1.3 && $posit < 1.4) {
	$tind = .54;
    } elsif($posit >=1.4 && $posit < 1.5) {
	$tind = .44;
    } elsif($posit >=1.5 && $posit < 1.6) {
	$tind = .69;
    } elsif($posit >=1.6 && $posit < 1.7) {
	$tind = .68;
    } elsif($posit >= 1.7 && $posit < 1.8) {
	$tind = .91;
    } elsif($posit >= 1.8 && $posit < 1.9) {
	$tind = .97;
    } elsif($posit >=1.9) {
	$tind = .97;
    }
    $Pmj += $tind*.33;

    if($contt >= 0. && $contt < .17) {
	$tind = .58;
    } elsif($contt >= .17  && $contt < .19) {
	$tind = .51;
    } elsif($contt >= .19 && $contt < .21) {
	$tind = .69;
    } elsif($contt >= .21 && $contt < .23) {
	$tind = .56;
    } elsif($contt >= .23 && $contt < .25) {
	$tind = .75;
    } elsif($contt >= .25 && $contt < .27) {
	$tind = .55;
    } elsif($contt >= .27 && $contt < .29) {
	$tind = .40;
    } elsif($contt >= .29 && $contt < .31) {
	$tind = .39;
    } elsif($contt >= .31 && $contt < .33) {
	$tind = .24;
    } elsif($contt >= .33 && $contt < .99) {
	$tind = .28;
    }
    $Pmj += $tind*.14;

    if($posig >= 0 && $posig < 1.1) {
	$tind = .08;
    } elsif($posig >= 1.1  && $posig < 1.2) {
	$tind = .08;
    } elsif($posig >= 1.2 && $posig < 1.3) {
	$tind = .16;
    } elsif($posig >= 1.3 && $posig < 1.4) {
	$tind = .27;
    } elsif($posig >=1.4 && $posig < 1.5) {
	$tind = .48;
    } elsif($posig >=1.5 && $posig < 1.6) {
	$tind = .53;
    } elsif($posig >= 1.6 && $posig < 1.7) {
	$tind = .64;
    } elsif($posig >= 1.7 && $posig < 1.8) {
	$tind = .74;
    } elsif($posig >= 1.8 && $posig < 1.9) {
	$tind = .88;
    } elsif($posig >=1.9) {
	$tind = .90;
    }
    $Pmj += $tind*.31;

    if($contg >= 0. && $contg < .17) {
	$tind = .29;
    } elsif($contg >= .17  && $contg < .19) {
	$tind = .33;
    } elsif($contg >= .19 && $contg < .21) {
	$tind = .41;
    } elsif($contg >= .21 && $contg < .23) {
	$tind = .41;
    } elsif($contg >= .23 && $contg < .25) {
	$tind = .73;
    } elsif($contg >= .25 && $contg < .27) {
	$tind = .64;
    } elsif($contg >= .27 && $contg < .29) {
	$tind = .64;
    } elsif($contg >= .29 && $contg < .31) {
	$tind = .47;
    } elsif($contg >= .31 && $contg < .33) {
	$tind = .54;
    } elsif($contg >= .33 && $contg < .99) {
	$tind = .40;
    }
    $Pmj += $tind*.15;

    if($posic >= 0 && $posic < 1.1) {
	$tind=.23;
    } elsif($posic >= 1.1  && $posic < 1.2) {
	$tind=.30;
    } elsif($posic >= 1.2 && $posic < 1.3) {
	$tind = .33;
    } elsif($posic >= 1.3 && $posic < 1.4) {
	$tind = .51;
    } elsif($posic >=1.4 && $posic < 1.5) {
	$tind = .48;
    } elsif($posic >=1.5 && $posic < 1.6) {
	$tind = .66;
    } elsif($posic >=1.6 && $posic < 1.7) {
	$tind = .81;
    } elsif($posic >= 1.7 && $posic < 1.8) {
	$tind = .7;
    } elsif($posic >= 1.8 && $posic < 1.9) {
	$tind = .7;
    } elsif($posic >=1.9) {
	$tind = .8;
    }
    $Pmj += $tind*.18;

    if($contc >= 0. && $contc < .17) {
	$tind = .31;
    } elsif($contc >= .17  && $contc < .19) {
	$tind = .39;
    } elsif($contc >= .19 && $contc < .21) {
	$tind = .44;
    } elsif($contc >= .21 && $contc < .23) {
	$tind = .43;
    } elsif($contc >= .23 && $contc < .25) {
	$tind = .59;
    } elsif($contc >= .25 && $contc < .27) {
	$tind = .59;
    } elsif($contc >= .27 && $contc < .29) {
	$tind = .64;
    } elsif($contc >= .29 && $contc < .31) {
	$tind = .51;
    } elsif($contc >= .31 && $contc < .33) {
	$tind = .64;
    } elsif($contc >= .33 && $contc < .99) {
	$tind = .82;
    }
    $Pmj += $tind*.12;
#    $Pmj *= 10;

    return $Pmj;
}


# qblast at ncbi

sub qblast {
    my $main=shift;
    my $canvas=shift;
    my $sheet = shift;
    my $annotation=shift;
    my $seqName=shift;
    my $sequenceName = $$seqName;
    my $sequence=shift;
    my $status_text=shift;
    
    my $program = shift; # what blast routine?

    my $fullSeq = shift;

    !defined $fullSeq && ($fullSeq = 0);

    # should be possible to select entire sequence - or mayhap even the seq. selection?
    # Get tags for current annotation 

    my ($annotation_id, $annotatorType );
    if(!$fullSeq) {
	($annotation_id, $annotatorType ) = main::annotation_id_type_selected($annotation, $sheet);
    } else {
	$annotation_id = "";
	$annotatorType = "";
    }

    my $annotated_sequence;
    my $annotation_nr;
    if (defined($annotation_id) && ($annotation_id ne "")) {
	print "Annotation $annotation_id made with $annotatorName{$annotatorType} was selected.\n";
	$annotation_nr = annotation_what_nr($annotation,$annotation_id);
	if($annotation_nr == -1) {
	    $WARNING && print "WARNING: no annotation_nr found for annotation_id $annotation_id! Aborting qblast..\n";
	    return -1;
	}
	
	$annotated_sequence = substr($$sequence,$$annotation{start}->[$annotation_nr]-1,$$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1);
	$DEBUG && print "Annotated sequence $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] is $annotated_sequence\n";
    } else {    
	print "Submitting entire sequence $sequenceName.";
	
	$annotated_sequence = $$sequence;
    }

    # If blastp, translate nt into aa sequence for the query

    my $annotated_protein;

    if($program eq "blastp") {
	# Will not work 
	if(!defined($annotation_id) or ($annotation_id eq "") or !defined $$annotation{frame}->[$annotation_nr]) {
	    $WARNING && print "WARNING: Attempted blastp query on a non-oriented annotation. This is not supposed to happen. Using blastx query instead.\n";
	    $program = "blastx";
	} else {
	    if($$annotation{frame}->[$annotation_nr]<3) {
		# Forward reading frame
		$annotated_protein=nt2aa($annotated_sequence); # We could remove M and + if we'd like..
	    } else {
		# Reverse reading frame
		$_=reverse(split(/ */,$annotated_sequence));
		tr/atgcATGC/tacgTACG/;
		$annotated_protein=nt2aa($_); # We could remove M and + if we'd like..
	    }
	    if (chop($annotated_protein) ne "+") {
		$WARNING && print "ERROR: Tried to chop a + away from $$annotation{id}->[$annotation_nr], but got something else!\nThis is not supposed to happen.\n";
	    }
	}
    }

    $$status_text="Connecting to NCBI...";
    $$main->update();

    # Open a socket to ncbi..
    my $remoteHostName="www.ncbi.nlm.nih.gov";
    my $remoteiaddr=inet_aton($remoteHostName);
    $remoteiaddr || die "Unable to resolve host $remoteHostName: $!\n";
    my $proto=getprotobyname('tcp');
    # my $remoteport=getservbyname('http','tcp');
    my $remoteport=80;
    my $remotepaddr=sockaddr_in($remoteport,$remoteiaddr);
    socket(SOCKET,PF_INET,SOCK_STREAM, $proto) || die "Could not create socket: $!\n";
    connect(SOCKET, $remotepaddr) || die "Could not connect: $!\n";
    
    # Construct quer
    my $postus="CMD=Put&PROGRAM=$program&DATABASE=nr"; # perhaps we wish COMPOSITION_BASED_STATISTICS=yes? That is not default in the blast proggy, but on the new NCBI-pages..
    
    ($annotation_id eq "") && ($postus.="&QUERY=>".$sequenceName."\015\012");
    ($annotation_id eq "") || ($postus.="&QUERY=>".$sequenceName."_"."$annotation_nr\015\012");

    if ($program eq "blastp") {
	$postus .= $annotated_protein;
    } else {
	$postus .= $annotated_sequence;
    }
    my %escapedCodes =();
    for (my $i=0;$i<=255;$i++) {
	$escapedCodes{chr($i)} = sprintf("%%%02X", $i);
    }
    $postus =~ s/([^A-Za-z0-9\-_.*&=])/$escapedCodes{$1}/g;
    
    my $contentlength = length($postus);

    # And POST it..
    $DEBUG && print "Attempting a POST\n"; # DEBUG
    $DEBUG && print "POST /blast/Blast.cgi HTTP/1.1\015\012Host: $remoteHostName\015\012Content-Type: application/x-www-form-urlencoded\015\012Content-Length: $contentlength\015\012\015\012$postus\015\012\n";
    
    select(SOCKET);
    print "POST /blast/Blast.cgi HTTP/1.1\015\012Host: $remoteHostName\015\012Content-Type:application/x-www-form-urlencoded\015\012Content-Length:$contentlength\015\012\015\012$postus\015\012";
    $|=1;
    select(STDOUT);
    
    my $ncbireply; # Save the output from ncbi as a string. Right now, we don't use it for anything, though..
    my $rid;
    my $waitingTime; # initial waiting time

    # Get Request id from output after post...
    # IMPROVEMENT: Do as background job with aid of fileevent?
    # BUG: failure to connect, but with some odd replys hangs program...
    while(<SOCKET>) {
	if(m/RID\ =/) {
	    ($rid)=m/RID\s+=\s+([\d\-\.\w]+)/; # RIDs now include ".BLASTQ\d"
	} elsif(m/RTOE = (\d+)/) {
	    $waitingTime=$1;
	}
	# IMPROVEMENT: print ncbi-estimate of computation time...
	$ncbireply.=$_;
	$DEBUG && print;
    }
    close SOCKET;
    
    $DEBUG && print "The request id was $rid\n";
    $$status_text="Submitted $program query to NCBI on sequence $sequenceName";
    if(!defined($annotation_id) || ($annotation_id eq "")) {
    } elsif(defined($$annotation{name}->[$annotation_nr]) && $$annotation{name}->[$annotation_nr] ne "") {
	$$status_text .= ".$$annotation{name}->[$annotation_nr]";
    } else {
	$$status_text .= $annotation_id;
    }
    $$status_text .= " and got request id $rid...";

    $$main->update();
    
    $DEBUG && print "Waiting...";
    $$main->after($waitingTime*1000,sub{check_waiting_time($main,$canvas,$sheet,$sequence,$annotation,$annotation_id,$rid,$program)});

}

sub check_waiting_time {
    # Note that check_waiting_time calls itself!
    my $main=shift;
    my $canvas=shift;
    my $sheet=shift;
    my $seq=shift;
    my $annotation=shift;
    my $annotation_id=shift;
    my $rid=shift;
    my $program=shift;

    # Connect to NCBI web server...
    my $remoteHostName="www.ncbi.nlm.nih.gov";
    my $remoteiaddr=inet_aton($remoteHostName);
    $remoteiaddr || die "Unable to resolve host $remoteHostName: $!\n";
    my $proto=getprotobyname('tcp');
    my $remoteport=80;
    my $remotepaddr=sockaddr_in($remoteport,$remoteiaddr);
    socket(NEWSOCKET,PF_INET,SOCK_STREAM, $proto) || die "Could not create socket: $!\n";
    connect(NEWSOCKET, $remotepaddr) || die "Could not connect: $!\n";
    
    # Set RID in query. Encoding not neccesary since RID only [0-9\-]
    my $postus="CMD=get&RID=$rid";
    my $contentlength=length($postus);

    # POST 
    select(NEWSOCKET);
    print "POST /blast/Blast.cgi HTTP/1.1\015\012Host: $remoteHostName\015\012Content-Type:application/x-www-form-urlencoded\015\012Content-Length:$contentlength\015\012\015\012$postus\015\012";
    $|=1;
    select(STDOUT);
    
    my $waitingTime=0; 
    my $ncbi_reply="";

    # Now check the waiting time..
    while(<NEWSOCKET>) { 
	# if(m/automatically updated in \<b\>(\d+)\<\/b\> seconds/) {
	if(m/RTOE = (\d+)/) {
	    $waitingTime=$1;
	} elsif (m/Status=WAITING/) {
	    $waitingTime=40;
	}
	$ncbi_reply.=$_; # Retain reply for potential analysis.
    }
    close NEWSOCKET;
    
    if($waitingTime) {
	# If a waiting time was announced, schedule this routine to run again after this 
	# time has elapsed..
	$$main->after($waitingTime*1000,sub{check_waiting_time($main,$canvas,$sheet,$seq,$annotation,$annotation_id,$rid,$program)});
	print ".";
	$DEBUG && print "DEBUG: Wait for ". ($waitingTime*1000) . " time units.\n";
	$$main->update();
	return 1;
    }

    $DEBUG && print "got reply!\n$ncbi_reply\n"; # DEBUG

    # Ok, we got a "real" answer this time. Parse it!
    if($program eq "blastn") {
	parse_blastn_reply($main, $canvas,$sheet,$seq,$annotation,$annotation_id,1e-5,$ncbi_reply);
    } elsif ($program eq "blastx") {
	parse_blastx_reply($main, $canvas,$sheet,$seq,$annotation,$annotation_id,0.001,$ncbi_reply);
    } elsif ($program eq "blastp") {
	parse_blastp_reply($main, $canvas,$sheet,$seq,$annotation,$annotation_id,0.001,$ncbi_reply);
    }

    return 0;
}

sub parse_blastn_reply {
    my $main=shift;
    my $canvas=shift;
    my $sheet=shift;
    my $seq=shift;
    my $annotation=shift;
    my $annotation_id=shift;
    my $blastpvalthreshold = shift;
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

    #$ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    $ncbi_reply =~ s/\x0d\x0a.+?\x0d\x0a//mg; # weird lines showing up in ncbi replys nowadays..
    $ncbi_reply=~s/\<.+?\>//mg;

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

    my @queryHitFrame;  # uninteresting (?), but used in latter steps..
    my @queryHitStrand;
    my $queryName;
    my @identities;
    my @alignLength;
    my @score;
    my @expect;
    
    my @alignment;
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
		    # next row for wu-blast result syntax
		    defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/;
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
		# next row for wu-blast result syntax
		defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/;
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


    for(my $i=0; $i<$nHits;$i++) {
	if($annotation_id ne "") {
	    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
	    $queryHitFrame[$i] = $$annotation{frame}->[$annotation_nr];
	} else {
	    $queryHitFrame[$i] = 0;
	}
	
	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i] ($queryHitFrame[$i])) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";
    }

    blastresult_interact($main,$canvas,$annotation,$annotation_id,$sheet,$seq,'blastn',\$queryName,\@queryHitBegin,\@queryHitEnd,\@queryHitStrand,\@queryHitFrame,\@subjectName,\@subjectHitBegin,\@subjectHitEnd,\@subjectHitStrand,\@identities,\@alignLength,\@score,\@expect,\@alignment,$blastpvalthreshold);

    # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
    # Suggest name accordning to Gen. nomenclature "whitepaper"?
}

sub parse_blastx_reply {
    my $main = shift;
    my $canvas = shift;
    my $sheet = shift;
    my $seq = shift;
    my $annotation = shift;
    my $annotation_id = shift;
    my $blastpvalthreshold = shift;
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

    #$ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    $ncbi_reply =~ s/\x0d\x0a.+?\x0d\x0a//mg; # weird lines showing up in ncbi replys nowadays..
    $ncbi_reply=~s/\<.+?\>//mg;
    $WARNING && $ncbi_reply eq "" && print "WARNING: empty blast response recieved in parse_blastx_reply.\n";

    # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

    my $nHits=0;
    my $inHit=0;
    my $inDetail=0;
    my $detailField=0;
    my $inHitName = 0;

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

    my @alignment;
    my $next_row_is_alignment = 0;
    
    foreach $_ (split(/\n{1}/,$ncbi_reply)) {    

	# TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
	
	# Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
	# The rest of the output is kindly enough different.. =)
	if($inHit) {
	    if($inHitName) {
		if( /Length\s+\=\s+\d+/ ) {
		    $inHitName = 0;
		} else {
		    my ($subjectNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		    $subjectName[$nHits-1]  .= " $subjectNamePart";
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
		    # next row for wu-blast result syntax
		    defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/;
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
		    ($expect[$nHits-1])=/Expect[\(\)\d\s]*\=\s+([0-9\.e\-\+]+)/; # how about those mulit-seg hits?
		    (($expect[$nHits-1]) =~ /^e/) and ($expect[$nHits-1] = "1$expect[$nHits-1]");
		} elsif(/^\>/) {
		    # Parse of alignment rows found hit on new subject sequence.
		    $inDetail=0;
		    $inHitName = 1;
		    # Hits supposedly begin with a row "> FASTA_NAME"
		    ($subjectName[$nHits])=/^\>(.+)/;
		    $nHits++;
		    $alignment[$nHits-1]="";
		} 
		# in detail ends
		#Parse a hit header..     
	    } elsif(/Score\s{1}\=/) {
		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		# next row for wu-blast result syntax
		defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/;

		# Syntax of the P value varies btw runs in blastn.. *sigh*
		($expect[$nHits-1])=/Expect[\(\)\d\s]*=\s+([0-9\.e\-\+]+)/;
		(($expect[$nHits-1]) =~ /^e/) and ($expect[$nHits-1] = "1$expect[$nHits-1]");
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
	    $inHit=1;
	    $inHitName = 1;
	    $nHits++;
	    $alignment[$nHits-1]="";
	} elsif(/^Query\=/) {
	    # Actually just evaluated once to get query name..
	    ($queryName)=/^Query\=\s*(.+)/;
#      print "Query $queryName.\n";
	}
    }

    for(my $i=0; $i<$nHits;$i++) {
	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";
    }

    blastresult_interact($main,$canvas,$annotation,$annotation_id,$sheet,$seq,'blastx',\$queryName,\@queryHitBegin,\@queryHitEnd,\@queryHitFrame,\@subjectName, \@subjectHitBegin, \@subjectHitEnd,\@identities,\@alignLength,\@score,\@expect,\@alignment,$blastpvalthreshold);
    # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
    # Suggest name accordning to Gen. nomenclature "whitepaper"?

}

sub blastresult_interact {
    my $main = shift;
    my $canvas = shift;
    my $annotation =shift;
    my $annotation_id= shift;
    my $sheet = shift;
    my $seq = shift;
    my $program = shift; # non-ref
    my $queryName=shift; 
    my $queryHitBegin = shift;
    my $queryHitEnd=shift;

    my $queryHitStrand;
    my $queryHitFrame;
    if($program eq "blastn") {
	$queryHitStrand=shift;
    }

    # elsif($program eq "blastx" or $program eq "blastp" or $program eq "cdblast") {
    $queryHitFrame=shift;

    my $subjectName=shift;
    my $subjectHitBegin=shift;
    my $subjectHitEnd=shift;

    my $subjectHitStrand;
    if($program eq "blastn") {
	$subjectHitStrand = shift;
    }
    
    my $cdLength;
    my $percentAligned;
    my $identities;
    my $alignLength;
    if($program eq "cdblast") {
	$cdLength=shift;
	$percentAligned=shift;
    } else {
	$identities=shift;
	$alignLength=shift;
    }

    my $score=shift;
    my $expect=shift;
    my $alignment=shift;

    my $blastpvalthreshold = shift;

    my $nHits = @$score; # anyone would do

    if($$sheet{blastpopup}) {
	# Then pop up a user interaction window for choosing hits for use in annotation.
	my $qblast_win=$$main->Toplevel;
	$qblast_win->title("$program results");
	# $qblast_win->geometry('+20+0');
	$qblast_win->configure(-background=>'linen',-width=>'600');
	
	my $qblast_main_frame=$qblast_win->Frame(-background=>$$sheet{default_win_background})->pack(-fill => 'both', -expand=> 'yes');
	my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
	my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
	for(my $i=0; $i<$nHits;$i++) {
	    my $hitentry;
	    if($program eq "blastn") {
		$hitentry="$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitStrand[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ($$subjectHitStrand[$i]) ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "blastx") {
		$hitentry="$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitFrame[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "blastp") {
		$hitentry = "$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitFrame[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "cdblast") {
		$hitentry = "$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i]  hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] cd-length=$$cdLength[$i]($$percentAligned[$i]% aligned) score=$$score[$i] Expect=$$expect[$i]";
	    }
	    $DEBUG && print "DEBUG: $hitentry\n";
	    $qblast_hits_list->insert('end',$hitentry);
	}  
	
	my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
	my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
	$qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
	# scrollbars... connect commands...
	
	my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
	if($annotation_id ne "") {
	    my $qblast_annotate=$qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-command=>sub { # Q&D for finding our annotation id...
		#$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
		my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
		my $note = $$annotation{note}->[$annotation_nr];
		foreach my $selected ($qblast_hits_list->curselection) {
		    if($program eq "blastx") {
			$note .= "$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] ($$queryHitFrame[$selected]) hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		    } elsif($program eq "blastp") {
			$note.="$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] ($$queryHitFrame[$selected]) hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		    } elsif($program eq "blastn") {
			$note.="$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] ($$queryHitStrand[$selected]) hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ($$subjectHitStrand[$selected]) ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		    } elsif($program eq "cdblast") {
			$note.="$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] cd-length=$$cdLength[$selected]($$percentAligned[$selected]% aligned) score=$$score[$selected] Expect=$$expect[$selected]\n";
		    }
		}
	      main::annotation_edit($main,$canvas,$annotation,$sheet,$annotation_nr,$note);
	    },-text=>"Annotate query")->pack(-side=>'left');
	}
	my $qblast_manual=$qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-text=>"Add as blasthits",-command=>sub {
	    foreach my $selected ($qblast_hits_list->curselection) {
		my $new_start;
		my $new_stop;

		if($annotation_id ne "") {
		    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);

		    my $color;
		    if($$expect[$selected] < 1e-90) {
			$color = 'DarkOrange1';
		    } elsif($$expect[$selected] < 1e-60) {
			$color = 'orange1';
		    } elsif($$expect[$selected] < 1e-40) {
			$color = 'orange2';
		    } elsif($$expect[$selected] < 1e-12) {
			$color = 'orange3';
		    } elsif($$expect[$selected] < 1e-5) {
			$color = 'DarkOrange3';
		    }

		    if($program eq "blastx") {
			if($$queryHitFrame[$selected] > 0) { # frame is only dep on the hit, since the actual query seq is always "forward"
			    $new_start =  $$annotation{start}->[$annotation_nr]+$$queryHitBegin[$selected]-1; # query coords in nt
			    $new_stop =  $$annotation{start}->[$annotation_nr]+$$queryHitEnd[$selected]-1;
			    $DEBUG && print "DEBUG: Set new stop to $new_stop and start to $new_start since $program on annotation $annotation_nr and query hit frame $$queryHitFrame[$selected].\n";
			} else {
			    $new_stop = $$annotation{stop}->[$annotation_nr]-$$queryHitEnd[$selected]+1; # query coords in nt
			    $new_start = $$annotation{stop}->[$annotation_nr]-$$queryHitBegin[$selected]+1;
			    $DEBUG && print "DEBUG: Set new stop to $new_stop and start to $new_start since $program on annotation $annotation_nr and query hit frame $$queryHitFrame[$selected].\n";
			}
			# add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$annotation{frame}->[$annotation_nr],$$subjectName[$selected],$$subjectHitBegin[$selected],$$subjectHitEnd[$selected],$$cdLength[$selected],$$percentAligned[$selected],$$score[$selected],$$expect[$selected],$$alignment[$selected], $color, 1);
		    } elsif ($program eq "blastp") {

			if($$annotation{frame}->[$annotation_nr] < 3) {
			    $new_start = $$annotation{start}->[$annotation_nr]+3*($$queryHitBegin[$selected]-1);
			    $new_stop = $$annotation{start}->[$annotation_nr]+3*$$queryHitEnd[$selected];
			} else {
			    # $DEBUG && print "Rev frame annotation start is $$annotation{start}->[$annotation_nr], stop $$annotation{stop}->[$annotation_nr] query hit begin $queryHitBegin[$selected] and end $queryHitEnd[$selected].\n";
			    $new_stop = $$annotation{stop}->[$annotation_nr]-3*($$queryHitBegin[$selected]-1);
			    $new_start = $$annotation{stop}->[$annotation_nr]-3*($$queryHitEnd[$selected]);
			}
			# add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$annotation{frame}->[$annotation_nr],$$subjectName[$selected],$$subjectHitBegin[$selected],$$subjectHitEnd[$selected],$$cdLength[$selected],$$percentAligned[$selected],$$score[$selected],$$expect[$selected],$$alignment[$selected], $color, 1);
		    } elsif ($program eq "blastn") {

			if($$annotation{frame}->[$annotation_nr] < 3) {
			    $new_start = $$annotation{start}->[$annotation_nr]+$$queryHitBegin[$selected]-1;
			    $new_stop = $$annotation{start}->[$annotation_nr]+$$queryHitEnd[$selected]-1;
			} else {
			    $new_stop = $$annotation{stop}->[$annotation_nr]-$$queryHitBegin[$selected]+1;
			    $new_start = $$annotation{stop}->[$annotation_nr]-$$queryHitEnd[$selected]+1;
			} 

		    } elsif ($program eq "cdblast") {

			if($$annotation{frame}->[$annotation_nr] < 3) { 
			    $new_start = $$annotation{start}->[$annotation_nr]+3*($$queryHitBegin[$selected]-1);
			    $new_stop = $$annotation{start}->[$annotation_nr]+3*$$queryHitEnd[$selected];
			} else {
			    $new_stop = $$annotation{stop}->[$annotation_nr]-3*($$queryHitBegin[$selected]-1);
			    $new_start = $$annotation{stop}->[$annotation_nr]-3*($$queryHitEnd[$selected]);
			}
		    }
		    $DEBUG && print "DEBUG: Set new stop to $new_stop and start to $new_start (for annotation $annotation_nr with $program).\n";
		} else { # if annotation_id is blank, the query is the entire sequence
		    if($program eq "blastx") {
			if($$queryHitFrame[$selected] > 0) {
			    $new_start = $$queryHitBegin[$selected]-1; # query coords in nt
			    $new_stop = $$queryHitEnd[$selected];
			    $DEBUG && print "DEBUG: Set new stop to $new_stop and start to $new_start since $program on whole seq (length ".length($$seq).").\n";
			} else {
			    $new_stop = $$queryHitBegin[$selected]-1; # query coords in nt
			    $new_start =$$queryHitEnd[$selected];
			    $DEBUG && print "DEBUG: Set new stop to $new_stop and start to $new_start since $program on whole seq (length ".length($$seq).").\n";
			}
		    } elsif ($program eq "blastn") {
			if($$subjectHitStrand[$selected] eq "Plus") {
			    $new_start =$$queryHitBegin[$selected]-1;
			    $new_stop = $$queryHitEnd[$selected]-1;
			} elsif($$subjectHitStrand[$selected] eq "Minus")  {
			    $new_stop = $$queryHitBegin[$selected]+1;
			    $new_start = $$queryHitEnd[$selected]+1;
			} else {
			    print "Could not parse subjectHitStrand (\"$subjectHitStrand[$selected]\") of blasthit with annotation_nr $selected.";
			}
			#add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitStrand[$selected],$$queryHitFrame[$selected],$$subjectName[$selected],$$subjectHitBegin[$selected],$$subjectHitEnd[$selected],$$subjectHitStrand[$selected],$$identities[$selected],$$alignLength[$selected],$$score[$selected],$$expect[$selected],$$alignment[$selected]);
		    }
		}

		if($program eq "blastn") {
		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitStrand[$selected],$$queryHitFrame[$selected],$$subjectName[$selected],$$subjectHitBegin[$selected],$$subjectHitEnd[$selected],$$subjectHitStrand[$selected],$$identities[$selected],$$alignLength[$selected],$$score[$selected],$$expect[$selected],$$alignment[$selected]);
		} else { # blastx, blastp 
		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitFrame[$selected],$$subjectName[$selected], $$subjectHitBegin[$selected],$$subjectHitEnd[$selected],$$identities[$selected],$$alignLength[$selected],$$score[$selected],$$expect[$selected],$$alignment[$selected]);
		}
	    }
	  main::level_layout($annotation, 'blast');
	  main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..
	})->pack(-side=>'left');

	my $qblast_viewalign=$qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-text=>"Display alignments",-command=>sub {
	    # fancy view?
	    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
	    
	    print "Alignments ($program) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";
	    
	    foreach my $selected ($qblast_hits_list->curselection) {
		if($program eq "blastx") {
		    print "$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] ($$queryHitFrame[$selected]) hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		} elsif($program eq "blastp") {
		    print "$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		} elsif($program eq "blastn") {
		    print "$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] ($$queryHitStrand[$selected]) hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] ($$subjectHitStrand[$selected]) ids=$$identities[$selected]/$$alignLength[$selected] score=$$score[$selected] Expect=$$expect[$selected]\n";
		} elsif ($program eq "cdblast") {
               #		foreach my $selected ($qblast_hits_list->curselection) { # duhh???
		    print "$$queryName $$queryHitBegin[$selected]-$$queryHitEnd[$selected] hit $$subjectName[$selected] $$subjectHitBegin[$selected]-$$subjectHitEnd[$selected] CD-length=$$cdLength[$selected]($$percentAligned[$selected]% aligned) Score=$$score[$selected] Expect=$$expect[$selected]\n";
		}
		print "\n$$alignment[$selected]\n";
	    }
	})->pack(-side=>'left');

	my $qblast_cancel=$qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');
    } else {
	# No popup - reconfigure display to a full blast align and autoadd.
	# my %orderSave = %annotatorOrder;
	my %heightSave = %annotatorHeight;
	
	my %sheetSave = ();

	# %annotatorOrder = (ORF => 1, ST => 1, blast => 2);
	%annotatorHeight = (ORF => 3, ST => 3, blast => 27);

	# parse through the annotation levels, save old settings and then disable everything xcept orf/st and blast
	foreach my $sheetkey (keys %$sheet) {
	    if($sheetkey =~ m/^display_\w+_level/) {
		$sheetSave{$sheetkey} = $$sheet{$sheetkey};
		$$sheet{$sheetkey} = 0;
	    }
	}
	
	$$sheet{display_ORF_level} = 1;
	$$sheet{display_ST_level} = 1;
	$$sheet{display_blast_level} = 1;
	
	for(my $i=0; $i<$nHits;$i++) {

	    # program specific!
	    my $hitentry;
	    if($program eq "blastn") {
		$hitentry="$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitStrand[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ($$subjectHitStrand[$i]) ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "blastx") {
		$hitentry="$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitFrame[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "blastp") {
		$hitentry = "$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    } elsif ($program eq "cdblast") {
		$hitentry = "$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] cd-length=$$cdLength[$i]($$percentAligned[$i]% aligned) score=$$score[$i] Expect=$$expect[$i]";
	    }

	    # my $hitentry = "$$queryName $$queryHitBegin[$i]-$$queryHitEnd[$i] ($$queryHitFrame[$i]) hit $$subjectName[$i] $$subjectHitBegin[$i]-$$subjectHitEnd[$i] ids=$$identities[$i]/$$alignLength[$i] score=$$score[$i] Expect=$$expect[$i]";
	    if ($$expect[$i] =~ /^e/) {
		$$expect[$i] = "1".$$expect[$i];
	    }

	    if($$expect[$i] < $blastpvalthreshold) {
		$DEBUG && print "DEBUG: Keep $hitentry\n";
		
		my $color;

		if($$expect[$i] < 1e-90) {
		    $color = 'DarkOrange1';
		} elsif($$expect[$i] < 1e-60) {
		    $color = 'orange1';
		} elsif($$expect[$i] < 1e-40) {
		    $color = 'orange2';
		} elsif($$expect[$i] < 1e-12) {
		    $color = 'orange3';
		} elsif($$expect[$i] < 1e-5) {
		    $color = 'DarkOrange3';
		}

		my $new_start;
		my $new_stop;

		if($program eq "blastx") {
		    if($annotation_id ne "") {
			#$WARNING && print "blastx of a protein type annotation..";
			my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
		    
			if($$queryHitFrame[$i] > 0) { #
			    $new_start = $$annotation{start}->[$annotation_nr]+$$queryHitBegin[$i]-1;
			    $new_stop = $$annotation{start}->[$annotation_nr]+$$queryHitEnd[$i]-1;
			} else {
			    $new_stop = $$annotation{stop}->[$annotation_nr]-$$queryHitEnd[$i]+1;
			    $new_start = $$annotation{stop}->[$annotation_nr]-$$queryHitBegin[$i]+1;
			    $DEBUG && print "Rev frame ($$queryHitFrame[$i]) annotation start is ", $$annotation{start}->[$annotation_nr], ", stop ", $$annotation{stop}->[$annotation_nr], " query hit begin $$queryHitBegin[$i] and end $$queryHitEnd[$i] so new start $new_start and new stop $new_stop.\n";
			}
		    } else {
			if($$queryHitFrame[$i] > 0) { # any digits on the blastx replys??
			    $new_start = $$queryHitBegin[$i];  # query coords in nt 
			    $new_stop = $$queryHitEnd[$i];
			} else {
			    $new_stop = $$queryHitBegin[$i] ;  # query coords in nt
			    $new_start = $$queryHitEnd[$i];
			}
		    }
		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitFrame[$i],$$subjectName[$i], $$subjectHitBegin[$i],$$subjectHitEnd[$i],$$identities[$i],$$alignLength[$i],$$score[$i],$$expect[$i],$$alignment[$i],$color,1);

		} elsif ($program eq "blastp") {
		    if($annotation_id ne "") {
			my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
		    
			if($$annotation{frame}->[$annotation_nr] < 3) {
			    $new_start = $$annotation{start}->[$annotation_nr]+3*($$queryHitBegin[$i]-1);
			    $new_stop = $$annotation{start}->[$annotation_nr]+3*$$queryHitEnd[$i];
			} else {
			    # $DEBUG && print "Rev frame annotation start is $$annotation{start}->[$annotation_nr], stop $$annotation{stop}->[$annotation_nr] query hit begin $queryHitBegin[$i] and end $queryHitEnd[$i].\n";
			    $new_stop = $$annotation{stop}->[$annotation_nr]-3*($$queryHitBegin[$i]-1);
			    $new_start = $$annotation{stop}->[$annotation_nr]-3*($$queryHitEnd[$i]);
			}
		    } else {
			$WARNING && print "WARNING: blank annotation id encountered for program $program in blastresult_interact. Running $program like this seems a rather dubious thing to do, unless the entire sequence is a protein.\n";
			$new_start = 3*($$queryHitBegin[$i]-1);
			$new_stop = 3*$$queryHitEnd[$i];
		    }

		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitFrame[$i],$$subjectName[$i], $$subjectHitBegin[$i],$$subjectHitEnd[$i],$$identities[$i],$$alignLength[$i],$$score[$i],$$expect[$i],$$alignment[$i],$color,1);
		} elsif ($program eq "blastn") { # blastn could possibly be invoked for the entire sequence
		    if($annotation_id ne "") {
			my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
			if($$annotation{frame}->[$annotation_nr] < 3) {
			    $new_start = $$annotation{start}->[$annotation_nr]+$$queryHitBegin[$i]-1;
			    $new_stop = $$annotation{start}->[$annotation_nr]+$$queryHitEnd[$i]-1;
			} else {
			    $new_stop = $$annotation{stop}->[$annotation_nr]-$$queryHitBegin[$i]+1;
			    $new_start = $$annotation{stop}->[$annotation_nr]-$$queryHitEnd[$i]+1;
			}
		    } else {
			if($$subjectHitStrand[$i] eq "Plus") {
			    $new_start =$$queryHitBegin[$i]-1;
			    $new_stop = $$queryHitEnd[$i]-1;
			} elsif($$subjectHitStrand[$i] eq "Minus")  {
			    $new_stop = $$queryHitBegin[$i]+1;
			    $new_start = $$queryHitEnd[$i]+1;
			} else {
			    print "Could not parse subjectHitStrand (\"$subjectHitStrand[$i]\") of blasthit with annotation_nr $annotation_nr.";
			}
		    }

		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$queryHitStrand[$i],$$queryHitFrame[$i],$$subjectName[$i],$$subjectHitBegin[$i],$$subjectHitEnd[$i],$$subjectHitStrand[$i],$$identities[$i],$$alignLength[$i],$$score[$i],$$expect[$i],$$alignment[$i], $color, 1);
		} elsif ($program eq "cdblast") {
		    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);

		    my $new_start;
		    my $new_stop;
		    if($$annotation{frame}->[$annotation_nr] < 3) { 
			$new_start = $$annotation{start}->[$annotation_nr]+3*($$queryHitBegin[$i]-1);
			$new_stop = $$annotation{start}->[$annotation_nr]+3*$$queryHitEnd[$i];
		    } else {
			$new_stop = $$annotation{stop}->[$annotation_nr]-3*($$queryHitBegin[$i]-1);
			$new_start = $$annotation{stop}->[$annotation_nr]-3*($$queryHitEnd[$i]);
		    }

		    add_blasthit($canvas,$annotation,$sheet,$seq,$program,$$queryName,$new_start,$new_stop,$$annotation{frame}->[$annotation_nr],$$subjectName[$i],$$subjectHitBegin[$i],$$subjectHitEnd[$i],$$cdLength[$i],$$percentAligned[$i],$$score[$i],$$expect[$i],$$alignment[$i], $color, 1);
		}
	    } else {
		$DEBUG && print "DEBUG: Drop (pval) $hitentry\n";
	    }
	}

      main::draw_seq($canvas,$main::seqName,$seq,$sheet); # using GLOBAL seqName
	$$canvas->configure(-height=>$$sheet{canvas_seq_height});
      main::level_layout($annotation, 'blast');
	# $DEBUG && print "DEBUG && redrawing annotations with annotation levels as follows: \n";
	# foreach $order (keys %annotatorOrder) { # DEBUG!!
#	  print "$order $annotatorOrder{$order}\n";
	#     }

      main::redraw_annotations($canvas,$annotation,$sheet,$seq);


	if(!Exists($qblast_norm_win)) {
	    # popup "view_normalise" window
	    $qblast_norm_win = $$main->Toplevel;
	    $qblast_norm_win->title("$program results");
	    #  $qblast_norm_win->geometry('+20+0');
	    $qblast_norm_win->configure(-background=>'linen',-width=>'300');
	    
	    my $qblast_frame = $qblast_norm_win->Frame(-background=>$$sheet{default_win_background})->pack(-fill => 'both', -expand=> 'yes');
	    $qblast_frame->Label(-background=>$$sheet{default_win_background},-text=>"Main display now showing $program results.\nPlease remove bad alignments and then hit Ok to keep the remaining hits.\nCancel will remove all of the current hits.")->pack(-fill => 'x',-expand=> 'yes',-side=>'top',-anchor=>'w');
	    my $qblast_action_frame = $qblast_frame->Frame(-background=>$$sheet{default_win_background})->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
	    my $qblast_keep = $qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-text=>"Ok",-command=>sub {
		permanent_temporary($annotation, 'blast');
		$qblast_norm_win->destroy;
	    })->pack(-side=>'left');
	    my $qblast_cancel = $qblast_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub { $qblast_norm_win->destroy; })->pack(-side=>'right');

	    # bind remove of any remaining temp annotations & restoral of original view to close of qblast window
	    $qblast_norm_win->waitWindow;
	    remove_temporary($annotation, 'blast', $sheet);
	    
	    # %annotatorOrder = %orderSave;
	    %annotatorHeight = %heightSave;

	    foreach my $sheetkey (keys %sheetSave) {
		$$sheet{$sheetkey} = $sheetSave{$sheetkey};
	    }

	  main::draw_seq($canvas,$main::seqName,$seq,$sheet); # using GLOBAL seqName
	    $$canvas->configure(-height=>$$sheet{canvas_seq_height});
	  main::level_layout($annotation, 'blast');
	  main::redraw_annotations($canvas,$annotation,$sheet,$seq);
	}
    }
}

sub permanent_temporary {
    my $annotation = shift;
    my $type = shift;    
    
    my $nr = $$annotation{nr};
    
    for (my $i = 0; $i < $nr; $i++) {
	defined($$annotation{temp}->[$i]) && $$annotation{temp}->[$i] && $$annotation{type}->[$i] eq $type && ($$annotation{temp}->[$i]=0);
    }

    return;
}

sub remove_temporary {
    my $annotation = shift;
    my $type = shift;
    my $sheet = shift;

    my $nr = $$annotation{nr};
    my @removeus;

    for (my $i = 0; $i < $nr; $i++) {
	defined($$annotation{temp}->[$i]) && $$annotation{temp}->[$i] && $$annotation{type}->[$i] eq $type && push @removeus, $i;
    }

  main::remove_annotations($annotation, \@removeus,$sheet);

    return;
}

sub parse_blastp_reply {
    my $main=shift;
    my $canvas=shift;
    my $sheet=shift;
    my $seq=shift;
    my $annotation=shift;
    my $annotation_id=shift;
    my $blastpvalthreshold = shift;
    my $ncbi_reply=shift;

    # $ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    # Nowadays, a > is already there.. 
    # Ok, for a crude first attempt: remove other HTML-tags. We lose information this way, 
    # but on the other hand it is possible to use previously written WU-blast parsing code..
    $ncbi_reply =~ s/\x0d\x0a.+?\x0d\x0a//mg; # weird lines showing up in ncbi replys nowadays..
    $ncbi_reply=~s/\<.+?\>//mg;

    $DEBUG && print "DEBUG: $ncbi_reply\n";

    # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

    my $annotation_nr = annotation_what_nr($annotation,$annotation_id); # annotation_id *should* be non-empty..

    my $nHits=0;
    my $inHit=0;
    my $inDetail=0;
    my $detailField=0;

    my @subjectHitBegin;
    my @subjectHitEnd;
    my @subjectName;
    my @queryHitBegin;
    my @queryHitEnd;
    my @queryHitFrame; #uninteresting (methinks), but used in latter steps..
    my $queryName;
    my @identities;
    my @alignLength;
    my @score;
    my @expect;

    my @alignment;
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
		    $nHits++;
		    $queryHitFrame[$nHits-1]=$$annotation{frame}->[$annotation_nr];
		    $inDetail=0;
		    ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		    # next row for wu-blast result syntax
		    defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/; 
#	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
		    ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
		} elsif(/^\>/) {
		    # Parse of alignment rows found hit on new subject sequence.
		    $inDetail=0;
		    # Hits supposedly begin with a row "> FASTA_NAME"
		    ($subjectName[$nHits])=/^\>(.+)/;
		    $nHits++;
		    $queryHitFrame[$nHits-1]=$$annotation{frame}->[$annotation_nr];
		    $alignment[$nHits-1]="";
		} 
		# in detail ends
		#Parse a hit header..     
	    } elsif(/Score\s{1}\=/) {
		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		defined $score[$nHits-1] or ($score[$nHits-1])=/Score\s+=\s+(\d+)\s*\([\d\.]+ bits\)/; 
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
	    $queryHitFrame[$nHits-1]=$$annotation{frame}->[$annotation_nr];
	    $inHit=1;
	} elsif(/^Query\=/) {
	    # Actually just evaluated once to get query name..
	    ($queryName)=/^Query\=\s*(.+)/;
#      print "Query $queryName.\n";
	}
    }

    # blastresult_interact($canvas,$annotation,$sheet,$seq,'blastp',\$queryName,\@queryHitBegin,\@queryHitEnd,\@queryHitFrame,\@subjectName, \@subjectHitBegin, \@subjectHitEnd,\@identities,\@alignLength,\@score,\@expect,\@alignment);
#	    add_blasthit($canvas,$annotation,$sheet,$seq,'blastp',$queryName,$new_start,$new_stop,$$annotation{frame}->[$annotation_nr],$subjectName[$selected],$subjectHitBegin[$selected],$subjectHitEnd[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);

    for(my $i=0; $i<$nHits;$i++) {
	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";
    }

    blastresult_interact($main,$canvas,$annotation,$annotation_id,$sheet,$seq,'blastp',\$queryName,\@queryHitBegin,\@queryHitEnd,\@queryHitFrame,\@subjectName, \@subjectHitBegin, \@subjectHitEnd,\@identities,\@alignLength,\@score,\@expect,\@alignment,$blastpvalthreshold);
    # Suggest name accordning to Gen. nomenclature "whitepaper"?
}

sub qblastcd {
    my $main=shift;
    my $canvas=shift;
    my $sheet = shift;
    my $annotation=shift;
    my $seqName=shift;
    my $sequenceName = $$seqName;
    my $sequence=shift;
    my $status_text=shift;

    # Conserved Domain blast is not yet under the NCBI - Q-systqem...
    my $blastpvalthreshold = shift;
    !defined $blastpvalthreshold && ($blastpvalthreshold = 1e-4);

    # get id and type for currently  selected annotation
    my ($annotation_id, $annotatorType) = main::annotation_id_type_selected($annotation, $sheet);

    # blastcd will not work for entire sequence (ok, could submit all orfs...)
    $DEBUG && print "Annotation $annotation_id made with $annotatorName{$annotatorType} was selected.\n";
    
    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
    
    my $annotated_sequence = substr($$sequence,$$annotation{start}->[$annotation_nr]-1,$$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1);

    print "Annotated sequence $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] is $annotated_sequence\n";
    
    # Translate nt into aa sequence for the query

    my $annotated_protein;

    if(!defined($$annotation{frame}->[$annotation_nr])) {
	$WARNING && print "WARNING: Attempted Conserved Domain query on a non-oriented annotation. This is not supposed to happen. Translating as a forward strand annotation instead..\n";
	# Forward reading frame
	$annotated_protein=nt2aa($annotated_sequence);
    } elsif( $$annotation{frame}->[$annotation_nr]<3 ) {
	# Forward reading frame
	$annotated_protein=nt2aa($annotated_sequence);
    } else {
	# Reverse reading frame
	$_=reverse(split(/ */,$annotated_sequence));
	tr/atgcATGC/tacgTACG/;
	$annotated_protein=nt2aa($_);
    }
    
    if (chop($annotated_protein) ne "+") {
	print "ERROR: Tried to chop a + away from the predicted translation of $$annotation{id}->[$annotation_nr], but got something else!\n";
    }
    
    $$status_text="Connecting to NCBI...";
    $$main->update();

    # Open a socket to ncbi..
    my $remoteHostName="www.ncbi.nlm.nih.gov";
    my $remoteiaddr=inet_aton($remoteHostName);
    $remoteiaddr || die "Unable to resolve host $remoteHostName: $!\n";
    my $proto=getprotobyname('tcp');
    # my $remoteport=getservbyname('http','tcp');
    my $remoteport=80;
    my $remotepaddr=sockaddr_in($remoteport,$remoteiaddr);
    socket(SOCKET,PF_INET,SOCK_STREAM, $proto) || die "Could not create socket: $!\n";
    connect(SOCKET, $remotepaddr) || die "Could not connect: $!\n";

    # Construct query
    my $postus = "DATALIB=oasis_sap"; # All PSSM; both PFAM and SMART 
    $postus .= "&EXPECT=0.01&SMODE=0&INPUT_TYPE=\"fasta\"&GRAPH=0&PAIR=0"; # Parameters..
    
    $postus .= "&SEQUENCE=>".$sequenceName."_"."$annotation_nr\015\012";  
    $postus .= $annotated_protein;
    my %escapedCodes = ();
    for (my $i=0;$i<=255;$i++) {
	$escapedCodes{chr($i)} = sprintf("%%%02X", $i);
    }
    $postus =~ s/([^A-Za-z0-9\-_.*&=])/$escapedCodes{$1}/g;
    
    my $contentlength=length($postus);

    # And POST it..
    $DEBUG && print "Attempting a POST\n"; # DEBUG
    $DEBUG && print "POST /Structure/cdd/wrpsb.cgi HTTP/1.1\015\012Host: $remoteHostName\015\012Content-Type: application/x-www-form-urlencoded\015\012Content-Length: $contentlength\015\012\015\012$postus\015\012\n";
    
    select(SOCKET);
    print "POST /Structure/cdd/wrpsb.cgi HTTP/1.1\015\012Host: $remoteHostName\015\012Content-Type:application/x-www-form-urlencoded\015\012Content-Length:$contentlength\015\012\015\012$postus\015\012";
    $|=1;
    select(STDOUT);

    $$main->update();

    my $ncbireply; # Save the output from ncbi as a string. Right now, we don't use it for anything, though..

    # Get Request id from output after post...
    # IMPROVEMENT: Do as background job with aid of fileevent?
    # BUG: failure to connect, but with some odd replys hangs program...
    while(<SOCKET>) { 

	$ncbireply.=$_;
	
    }
    close SOCKET;

    $DEVEL && print "DEVEL: NCBI reply (unmodified):\n$ncbireply\n";

    parse_blastcd_reply($main, $canvas,$sheet,$sequence,$annotation,$annotation_id,$blastpvalthreshold,$ncbireply);  
}

sub parse_blastcd_reply {
    my $main=shift;
    my $canvas=shift;
    my $sheet=shift;
    my $seq=shift;
    my $annotation=shift;
    my $annotation_id=shift;
    my $blastpvalthreshold = shift;
    my $ncbi_reply=shift;

    # remove annoying tags in reply..
#    $ncbi_reply =~ s$</?FONT.*?>$$g;
#    $ncbi_reply =~ s$</?[pP][rR][eE]>$$g;
#    $ncbi_reply =~ s$</?TD.*?>$$g;
#    $ncbi_reply =~ s$</?TR.*?>$$g;
#    $ncbi_reply =~ s$</?TABLE.*?>$$g;
#    $ncbi_reply =~ s$</?FORM.*?>$$g;
#    $ncbi_reply =~ s$</?INPUT.*?>$$g;
#    $ncbi_reply =~ s$</?SELECT.*?>$$g;
#    $ncbi_reply =~ s$</?OPTION.*?>$$g;
#    $ncbi_reply =~ s$</?IMG.*?>$$g;
#    $ncbi_reply =~ s$</?[bB]>$$g;

    #$ncbi_reply =~ s$</?[aA]{1}.*?>$$g; # substitute a tags for >
    $ncbi_reply =~ s/<a\s+href\=\"http\:\/\/www\.ncbi\.nlm\.nih\.gov\/Structure\/cdd\/.*?>/>/mg;
    $ncbi_reply =~ s/\<.+?\>//mg;
    $ncbi_reply =~ s/\x0d\x0a.+?\x0d\x0a//mg; # weird lines showing up in ncbi replys nowadays..
    
    $DEBUG && print "Got (and pre-processed) reply from NCBI: $ncbi_reply\n";

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
    my @cdLength;
    my @percentAligned;
    my @score;
    my @expect;

    my @alignment;
    my $next_row_is_alignment = 0;
    
    foreach $_ (split(/\n{1}/,$ncbi_reply)) {    

	if(m/...No hits found!/) {
	    print "No RPS BLAST CD hits found for annotation (id $annotation_id).\n";
	    return;
	} 

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
		    ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		    ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
		} elsif(/CD-Length\s+\=/) {
		    # Parse of alignment rows found a new hit on the same subject sequence.
		    $DEBUG && print "DEBUG: new hit - found CD-Length row..\n";
		    $subjectName[$nHits]=$subjectName[$nHits-1];
		    $inDetail=0;
		    $nHits++;
		    ($cdLength[$nHits-1], $percentAligned[$nHits-1])=/CD-Length\s+\=\s+(\d+)\s+residues\,.*?\s+([\d\.]+)\%\s+aligned/;
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
	    } elsif(/CD-Length\s+\=/) {
		# Parse of alignment rows found a new hit on the same subject sequence.
		$DEBUG && print "DEBUG: new hit - found CD-Length row..\n";
		($cdLength[$nHits-1], $percentAligned[$nHits-1])=/CD-Length\s+\=\s+(\d+)\s+residues\,.*?([\d\.]+)\%\s+aligned/;		    
	    } elsif(/Score\s{1}\=/) {
		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
		# Syntax of the P value varies btw runs in blastn.. *sigh*
		($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
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
	    ($queryName)=/^Query\=\s*local sequence:\s*(.+)/; 
#      print "Query $queryName.\n";
	}
    }

    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);

    for(my $i=0; $i<$nHits;$i++) {
	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($$annotation{frame}->[$annotation_nr]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] cd-length=$cdLength[$i]($percentAligned[$i]% aligned) score=$score[$i] Expect=$expect[$i]\n";
    }

    blastresult_interact($main,$canvas,$annotation,$annotation_id,$sheet,$seq,'cdblast',\$queryName,\@queryHitBegin,\@queryHitEnd,$$annotation{frame}->[$annotation_nr],\@subjectName, \@subjectHitBegin, \@subjectHitEnd,\@cdLength,\@percentAligned,\@score,\@expect,\@alignment,$blastpvalthreshold);
}

sub add_blasthit {
    my $canvas = shift;
    my $annotation =shift;
    my $sheet = shift;
    my $seq = shift;
    my $program = shift;
    my $queryName=shift; 
    my $hitBegin=shift; 
    my $hitEnd=shift;

    my $queryHitStrand;
    my $queryHitFrame;
    if($program eq "blastn") {
	$queryHitStrand=shift;
    }

    # elsif($program eq "blastx" or $program eq "blastp" or $program eq "cdblast") {
    $queryHitFrame=shift;

    my $subjectName=shift; 
    my $subjectHitBegin=shift;
    my $subjectHitEnd=shift;

    my $subjectHitStrand;
    if($program eq "blastn") {
	$subjectHitStrand=shift;
    }
    
    my $cdLength;
    my $percentAligned;
    my $identities;
    my $alignLength;
    if($program eq "cdblast") {
	$cdLength=shift;
	$percentAligned=shift;
    } else {
	$identities=shift;
	$alignLength=shift;
    }

    my $score=shift; 
    my $expect=shift; 
    my $alignment=shift;

    my $color = shift;
    !defined($color) && ($color = "DarkOrange1");

    my $temporary = shift; # in case of e g a temporary addition we want to be able to locate these
    !defined($temporary) && ($temporary = 0);


    my $nr;
    if($$annotation{nr}) {
	$nr=$$annotation{nr};
    } else {
	$nr=0;
	$$annotation{nr}=0;
    }

    $$annotation{uid}->[$nr] = $$annotation{unr};
    $$annotation{unr}++;

    $$annotation{id}->[$nr]= "ab_$nr";

    $main::annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;
    $$annotation{seqName}->[$nr]= "";  #$$annotation{seqName}->[$annotation_nr]; get from the parsing routines.. 
    $$annotation{type}->[$nr]="blast";

    $$annotation{temp}->[$nr] = $temporary;

    $$annotation{start}->[$nr]=$hitBegin; # uh-huh? Nah, that would be the coordinate relative to the annotation that was searched!
    $$annotation{stop}->[$nr]=$hitEnd;

    if($program eq "blastn") {
	if($queryHitFrame == 6) {
	    $$annotation{frame}->[$nr] = 6;
	} elsif($queryHitStrand eq 'Minus') {
	    if($queryHitFrame < 3) {
		$$annotation{frame}->[$nr] = 3;
	    } elsif ($queryHitFrame > 3) {
		$$annotation{frame}->[$nr] = 0;
	    }
	} elsif($queryHitStrand eq 'Plus') {
	    $$annotation{frame}->[$nr] = $queryHitFrame;
	} else {
	    $WARNING && print "WARNING: A blast hit ($queryName hit $subjectName $subjectHitBegin-$subjectHitEnd) with unknown query hit orientation $queryHitStrand-$subjectHitStrand!\n";
	    $$annotation{frame}->[$nr] = 6;
	}
    } elsif ($program eq "blastx") {
	if($queryHitFrame =~ /[+]+/ )  {	   
	    $$annotation{frame}->[$nr] = 0;
	} elsif($queryHitFrame =~ /[-]+/ )  {
	    $$annotation{frame}->[$nr] = 3;
	} else {
	    $WARNING && print "WARNING: A blast hit ($queryName hit $subjectName $subjectHitBegin-$subjectHitEnd) with unknown query hit orientation $queryHitStrand!\n";
	    $$annotation{frame}->[$nr] = 6;
	}
	$DEVEL && $DEBUG && print "DEVEL DEBUG: Setting annotation frame for $nr to ", $$annotation{frame}->[$nr], " since query hit frame is $queryHitFrame (for hit from program $program).\n";
    } elsif ($program eq "blastp" or $program eq "cdblast") {
	# ok, should keep track of the orientation of the submitted annotatation here, I guess.. 
	$DEVEL && $DEBUG && print "DEVEL DEBUG: Setting annotation frame for $nr to $queryHitFrame (for hit from program $program).\n";
	$$annotation{frame}->[$nr] = $queryHitFrame;
    }

    my ($namepart) =  $subjectName =~ /(?:.+\|)*\s*(.+?)$/; # gi|blah|blah|blah| product [species]
    $$annotation{name}->[$nr]= $namepart;
    
    $$annotation{comment}->[$nr]="$program hit: $queryName $hitBegin-$hitEnd ";

    if(defined($queryHitStrand)) {
	$$annotation{comment}->[$nr] .= "($queryHitStrand) ";
    }
    $$annotation{comment}->[$nr] .= "hit $subjectName $subjectHitBegin-$subjectHitEnd ";
    if(defined($subjectHitStrand)) {
	$$annotation{comment}->[$nr] .= "($subjectHitStrand) ";
    }
    if($program eq "cdblast") {
	$$annotation{comment}->[$nr] .= "CD-length=$cdLength($percentAligned% aligned) ";
    } else {
	$$annotation{comment}->[$nr] .= "ids=$identities/$alignLength ";
    }
    $$annotation{comment}->[$nr] .= "score=$score Expect=$expect ";

    $$annotation{note}->[$nr]=$alignment;
    if(defined($color) && $color ne "") {
	$$annotation{color}->[$nr]=$color;
    } else {
	$$annotation{color}->[$nr]="DarkOrange1";
    }
    if($program eq 'blastn') {
	$$annotation{level}->[$nr] = 0;   
    } else {
	if($$annotation{frame}->[$nr] < 3) {
	    $$annotation{level}->[$nr] = $$annotation{frame}->[$nr];
	} else {
	    $$annotation{level}->[$nr] = $$annotation{frame}->[$nr] - 5;
	}
    }

    # $DEVEL && print "DEVEL: adding a blasthit $$annotation{id}->[$nr]:\n$$annotation{comment}->[$nr]\n$$annotation{note}->[$nr]\n.\n";

    $$annotation{nr}++;

    $$sheet{status} = "dirty";
}

sub run_external {
    my $main=shift;
    my $canvas=shift;
    my $status_text=shift;
    my $run_win=shift;
    my $runstate=shift;
    my $annotation=shift;
    my $seqName=shift; # ref
    my $seq=shift; # ref
    my $sheet=shift;
    
    # read from preferences
    my $tmppath = $$sheet{tmppath};
    my $binpath = $$sheet{binpath};
    my $blast2path = $$sheet{blast2path};
    my $ncbiblastpath = $$sheet{ncbiblastpath};
    my $testcodepath = $$sheet{testcodepath};
    
    my $session_rand=$$.int(rand 1000000); # The current pid and a random number to make the temp files a bit more unique

    # Running glimmer
    
    my $doGlim=0;
    my $glim="";

    # Is a Glimmer-run requested?
    if ( ($$runstate{glim_run_seq} or $$runstate{glim_run} ) and ($$runstate{glim_train} or $$runstate{glim_train_seq} or $$runstate{glim_model}) ) {
	$doGlim=1;
    } elsif ($$runstate{glim_run_seq} or $$runstate{glim_run} or $$runstate{glim_train} or $$runstate{glim_train_seq} or $$runstate{glim_model} ) { # xor, but using elsif.. :)
	$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"You must specify either both target sequence and training procedure, or none.");
	$doGlim=0;
    } else {
	# Glimmer not selected.
	$doGlim=0;
    }

    if($doGlim) {
	my $glimmer2path;
	my $glim;
	if($$runstate{glim_altstart}) {
	    # Hm, actually, it is glim.pl that should know about that.. Flag, or two different scripts?
	    $glimmer2path=$$sheet{glimmer2altstartpath};
	    $glim=$$sheet{binpath}."/".$$sheet{glimalstartbin};
	} else {
	    $glimmer2path=$$sheet{glimmer2path};
	    $glim=$$sheet{binpath}."/".$$sheet{glimbin};
	}
	
	if($$runstate{glim_train_seq}) {
	    # Save current seq to a temporary fasta file, and use that for training.
	    $$runstate{glim_train_file}="$tmppath/temp.$$seqName.$session_rand.fasta";
	    open SEQOUT,">$$runstate{glim_train_file}" || print "ERROR: could not unlink temp file $$runstate{glim_train_file}.\n";
	    print SEQOUT ">$$seqName\n$$seq\n";
	    close SEQOUT;
	    my $error=0;
	    $DEBUG && print "DEBUG: $glimmer2path/long-orfs \"$$runstate{glim_train_file}\" |$glimmer2path/get-putative > $tmppath/tmp.$session_rand.coord\n";
	    system("$glimmer2path/long-orfs \"$$runstate{glim_train_file}\" |$glimmer2path/get-putative > $tmppath/tmp.$session_rand.coord") == 0 or ($error=1);
	    $DEBUG && print "DEBUG: $glimmer2path/extract \"$$runstate{glim_train_file}\" $tmppath/tmp.$session_rand.coord >$tmppath/tmp.$session_rand.train\n";
	    system("$glimmer2path/extract \"$$runstate{glim_train_file}\" $tmppath/tmp.$session_rand.coord >$tmppath/tmp.$session_rand.train") == 0 or ($error=1);
	    unlink($$runstate{glim_train_file}) || print "ERROR: could not unlink temp file $$runstate{glim_train_file}.\n";
	    $DEBUG && print "DEBUG: $glimmer2path/build-icm <$tmppath/tmp.$session_rand.train >$tmppath/tmp.$session_rand.model\n";
	    system("$glimmer2path/build-icm <$tmppath/tmp.$session_rand.train >$tmppath/tmp.$session_rand.model") == 0 or ($error=1);
	    unlink("$tmppath/tmp.$session_rand.train", "$tmppath/tmp.$session_rand.coord");
	    $$runstate{glim_model_file}="$tmppath/tmp.$session_rand.model";
	    if($error) {
		$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"Glimmer training failed.");
		return 0;
	    }
	} elsif ($$runstate{glim_train}) {
	    # Train on selected file..
	    my $error=0;
	    $DEBUG && print "DEBUG: $glimmer2path/long-orfs \"$$runstate{glim_train_file}\" | $glimmer2path/get-putative > $tmppath/tmp.$session_rand.coord\n";
	    system("$glimmer2path/long-orfs \"$$runstate{glim_train_file}\" | $glimmer2path/get-putative > $tmppath/tmp.$session_rand.coord") == 0 or ($error=1);

	    system("$glimmer2path/extract \"$$runstate{glim_train_file}\" $tmppath/tmp.$session_rand.coord >$tmppath/tmp.$session_rand.train") == 0 or ($error=1);
	    $DEBUG && print "DEBUG: $glimmer2path/build-icm <$tmppath/tmp.$session_rand.train >$tmppath/tmp.$session_rand.model\n";
	    system("$glimmer2path/build-icm <$tmppath/tmp.$session_rand.train >$tmppath/tmp.$session_rand.model") == 0 or ($error=1);
	    if($error) {
		$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"Glimmer training failed.");
		return 0;
	    }
	    $$runstate{glim_model_file}="$tmppath/tmp.$session_rand.model";
	}

	if($$runstate{glim_run_seq}) {
	    # Save current seq to a temporary fasta file, and use that as prediction data.
	    open SEQOUT,">$tmppath/temp.$$seqName.$session_rand.fasta";
	    print SEQOUT ">$$seqName\n$$seq\n";
	    close SEQOUT;
	    $$runstate{glim_run_file}="$tmppath/temp.$$seqName.$session_rand.fasta";
	}

	$DEBUG && print "DEBUG: trying a system $glim \"$$runstate{glim_model_file}\" \"$$runstate{glim_run_file}\" $$runstate{glim_shortest_orf} \"$glimmer2path\" > $tmppath/temp.$session_rand.gui.glim\n";
	# Usage: glim.pl icmFile seqFastaFile minExonLengthq
	system("$glim \"$$runstate{glim_model_file}\" \"$$runstate{glim_run_file}\" $$runstate{glim_shortest_orf} \"$glimmer2path\" > $tmppath/temp.$session_rand.gui.glim") == 0 or ($$main->messageBox(-icon=>'error',-type=>'OK',-message=>"Glimmer run failed.") and (return 0)); # Obfuscated perl candidate. =) Note operator short-circuit properties.
	
	if($$runstate{glim_run_seq}) {
	  main::import_glimmer($main,"$tmppath/temp.$session_rand.gui.glim",$annotation,$$seqName,$seq,$sheet);
	    $DEBUG || unlink("$tmppath/temp.$session_rand.gui.glim","$tmppath/temp.$$seqName.$session_rand.fasta" );
	}

	if($$runstate{glim_train} || $$runstate{glim_train_seq}) {
	    unlink("$tmppath/tmp.$session_rand.model");
	}
    }

    # Running EST-processing

    # Running ESTORF
    my $doEstOrf;

    if(($$runstate{estorf_seq} or $$runstate{estorf_currentseq}) and ( $$runstate{estorf_est} or $$runstate{estorf_db})) {
	$doEstOrf=1;
    } elsif ( ($$runstate{estorf_seq} or $$runstate{estorf_currentseq}) or ($$runstate{estorf_est} or $$runstate{estorf_db})) {
	$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"You must specify either both genomic sequence and est sequence, or none.");
	$doEstOrf=0;
    } else {
	# ESTORF nor requested.
	$doEstOrf=0;
    }
    
    if($doEstOrf) {
	# Usage: pmestLean.pl genomicSeqFastaFileName estSeqFastaFileName 
	# internally calls:
	# Usage: matchEstAutoLean.pl estQuerySeqencesFasta subject(genomic)DbFileName (assume seqname eq dbfilebasename)
	$run_win->title('A GUI - writing sequence file');
	$run_win->update();
	if($$runstate{estorf_currentseq}) {
	    open SEQOUT,">$tmppath/temp.$$seqName.$session_rand.fasta";
	    print SEQOUT ">$$seqName\n$$seq\n";
	    close SEQOUT;
	    $$runstate{estorf_seq_file}="$tmppath/temp.$$seqName.$session_rand.fasta";
	}
	$run_win->title('A GUI - building blastdb...');
	$run_win->update();
	# build a blastable GENOMIC "db".
	
	if($$runstate{estorf_wu}) {
	    system("$blast2path/pressdb \"$$runstate{estorf_seq_file}\"");
	} elsif ($$runstate{estorf_ncbi}) {
	    system("$ncbiblastpath/formatdb -i \"$$runstate{estorf_seq_file}\" -p F -o T");
	}
	# Show a please-wait window
	$run_win->title('Please wait - running ESTORF match...');

	my $run_win_wait_dialog=$$main->Toplevel;
	$run_win_wait_dialog->title("A GUI - Running ESTORF match");
	# $run_win_wait_dialog->geometry('+0+45');
	$run_win_wait_dialog->configure(-background=>'linen');
	my $run_win_wait_dialog_frame=$run_win_wait_dialog->Frame->pack(-expand=>'yes',-fill=>'both');
	my $run_win_wait_dialog_label=$run_win_wait_dialog_frame->Label(-text => "Please wait - running ESTORF match...\n\nIf the EST file is large, this may require serious runtime.", -justify => 'left', -background=>$$sheet{default_win_background})->pack(-expand=>'yes',-fill=>"both");

	$$main->update();
	$run_win_wait_dialog->update(); # Could use this, background the run and produce a progressbar - that might be nice

	my $error = 0;
	if($$runstate{estorf_wu}) {
	    $DEBUG && print "DEBUG: System $binpath/pmestLean.pl \"$$runstate{estorf_seq_file}\" \"$$runstate{estorf_est_file}\" \"$binpath\" \"$blast2path\" \"$tmppath\" > \"$tmppath/$$seqName.$session_rand.pmest\"\n";
	    system("$binpath/pmestLean.pl \"$$runstate{estorf_seq_file}\" \"$$runstate{estorf_est_file}\" \"$binpath\" \"$blast2path\" \"$tmppath\" > \"$tmppath/$$seqName.$session_rand.pmest\"") == 0 or ($error = 1);
	    unlink("$tmppath/temp.$$seqName.$session_rand.fasta", "$tmppath/temp.$$seqName.$session_rand.fasta.nhd", "$tmppath/temp.$$seqName.$session_rand.fasta.ntb","$tmppath/temp.$$seqName.$session_rand.fasta.csq");
	} elsif($$runstate{estorf_ncbi}) {
	    $DEBUG && print "DEBUG: System $binpath/pmestMegablast.pl \"$$runstate{estorf_seq_file}\" \"$$runstate{estorf_est_file}\" \"$binpath\" \"$ncbiblastpath\" > \"$tmppath/$$seqName.$session_rand.pmest\"\n";
	    system("$binpath/pmestMegablast.pl \"$$runstate{estorf_seq_file}\" \"$$runstate{estorf_est_file}\" \"$binpath\" \"$ncbiblastpath\" > \"$tmppath/$$seqName.$session_rand.pmest\"") == 0 or ($error = 1);
	    unlink("$tmppath/temp.$$seqName.$session_rand.fasta", "$tmppath/temp.$$seqName.$session_rand.fasta.nhr", "$tmppath/temp.$$seqName.$session_rand.fasta.nin","$tmppath/temp.$$seqName.$session_rand.fasta.nsd","$tmppath/temp.$$seqName.$session_rand.fasta.nsi","$tmppath/temp.$$seqName.$session_rand.fasta.nsq","$tmppath/formatdb.log");
	}

	if($error == 1) {
	    $$main->messageBox(-icon=>'error',-type=>'OK',-message=>"ESTORF run failed $!.");
	    return 0;
	}

	if($$runstate{estorf_currentseq}) {
	    $run_win->title('A GUI - importing results...');
	    $DEBUG && print "DEBUG: importing results from $tmppath/$$seqName.$session_rand.pmest..\n";
	    $run_win->update();
	  main::import_estorf($main,"$tmppath/$$seqName.$session_rand.pmest",$annotation,$$seqName,$seq,$sheet); # seqName is a ref here, but import expects a string..
	  main::import_est($main,"$tmppath/$$seqName.$session_rand.pmest",$annotation,$$seqName,$seq,$sheet);
	    # unlink("$tmppath/$$seqName.$session_rand.pmest");
	    # LEAVE IT HERE FOR DEBUGGING!
	}
	
	# Close please-wait window
	$run_win_wait_dialog->destroy(); 

	# remove blast-"db"?
	$run_win->title('A GUI - run');
	$run_win->update();
    }

    my $doTestcode;

    if(($$runstate{testcode_run} or $$runstate{testcode_current_seq})) {
	$doTestcode=1;
    } else {
	# Testcode not requested.
	$doTestcode=0;
    }

    if($doTestcode) {

	$run_win->title('A GUI - writing sequence file');
	$run_win->update();

	if($$runstate{testcode_current_seq}) {
	    open SEQOUT,">$tmppath/temp.$$seqName.$session_rand.fasta";
	    print SEQOUT ">$$seqName\n$$seq\n";
	    close SEQOUT;
	    $$runstate{testcode_seq_file}="$tmppath/temp.$$seqName.$session_rand.fasta";
	}

	my ($testcode_forward_out, $testcode_reverse_out); 
	
	if($$runstate{testcode_current_seq}) {
	    ($testcode_forward_out, $testcode_reverse_out) = ("$tmppath/temp.$$seqName.$session_rand.testcode.f","$tmppath/temp.$$seqName.$session_rand.testcode.r");
	} else {
	    my $fasta_name = basename($$runstate{testcode_seq_file});
	    ($testcode_forward_out, $testcode_reverse_out) = ("$tmppath/$fasta_name.$session_rand.testcode.f","$tmppath/$fasta_name.$session_rand.testcode.r");
	}

	my $args = "-i \"$$runstate{testcode_seq_file}\" -f \"$testcode_forward_out\" -r \"$testcode_reverse_out\"";
	if($$runstate{testcode_mwinl}) {
	    $args .= " -m $$runstate{testcode_mwinl}";
	}
	if($$runstate{testcode_increment}) {
	    $args .= " -d $$runstate{testcode_increment}";
	}

	my $error = 0;

	$DEBUG && print "System: $testcodepath/testcode $args.\n";
	system("$testcodepath/testcode $args") == 0 or ($error=1);

	if($error == 1) {
	    $$main->messageBox(-icon=>'error',-type=>'OK',-message=>"TESTCODE run failed $!.");
	    return 0;
	}

	if($$runstate{testcode_current_seq}) {
	  main::import_testcode($main, "$testcode_forward_out", "$testcode_reverse_out",$annotation,$$seqName,$seq,$sheet);
	    # cleanup
	    # unlink("$$runstate{testcode_seq_file}", "$testcode_reverse_out", "$testcode_forward_out");
	}
    }

    my $doSplicemodel = 0;

    if($$runstate{sm_currentseq}) {
	$doSplicemodel = 1;
    }

    # check options
    if($$runstate{sm_report} && (!defined($$runstate{sm_report_file}) ||$$runstate{sm_report_file} eq "")) {
	$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"If you select a report, you must specify a report file.");
	$doSplicemodel = 0;
    }

    if($$runstate{sm_stats} && (!defined($$runstate{sm_stats_file}) ||$$runstate{sm_stats_file} eq "")) {
	$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"If you select statistics output, you must specify a statistics file.");
	$doSplicemodel = 0;
    }
    
    if($doSplicemodel) {
	my $options = "";

	if($$runstate{sm_report}) { 
	    $options .= " --report --reportfile ".$$runstate{sm_report_file};	    
	}

	if($$runstate{sm_stats}) { 
	    $options .= " --rout --statistics ".$$runstate{sm_stats_file};
	}

	if($$runstate{sm_short}) { 
	    $options .= " --short";
	}
	
	if($$runstate{sm_long}) { 
	    $options .= " --long";
	}

	if($$runstate{sm_polypyag}) {
	    $options .= " --polypy-ag ".$$runstate{sm_polypyag};
	}

	if($$runstate{sm_uorf}) {
	    $options .= " --uorf --uorf-file $tmppath/$$seqName.$session_rand.sm.uorf.gff";
	}
	

	# export genbank file - select levels, perhaps? popup? filename...
	my $gbkfile = main::file_export_genbank($main, $seqName, $seq, $annotation, $sheet);

	# add to rc-file..
	my $smpath = "/home/daniel/work/reggae";

	my $error = 0;
	$DEBUG && print STDERR "DEBUG: attempting $smpath/splicemodel.pl --infile $gbkfile --gff $tmppath/$$seqName.$session_rand.sm.gff --uorf --uorf-file $tmppath/$$seqName.$session_rand.sm.uorf.gff $options --polypy-ag 120";
	system("$smpath/splicemodel.pl --infile $gbkfile --gff $tmppath/$$seqName.$session_rand.sm.gff $options");

	# import hits & splicesite features
	my $bioframed = 0; # stacked gff annotations
        main::import_gff($main,"$tmppath/$$seqName.$session_rand.sm.gff",$bioframed,$annotation,$$seqName,$seq,$sheet);

	# import uORFs
        $$runstate{sm_uorf} && main::import_gff($main,"$tmppath/$$seqName.$session_rand.sm.uorf.gff",$bioframed,$annotation,$$seqName,$seq,$sheet);
    }



    my $doBlast;

    if(($$runstate{blast_seq} or $$runstate{blast_annotation} or $$runstate{blast_currentseq}) and ( $$runstate{blast_db_seq_file} or $$runstate{blast_db} or $$runstate{blast_db_remote})) {
	$doBlast=1;
    } elsif ( ($$runstate{blast_seq} or $$runstate{blast_currentseq}) or ($$runstate{blast_db_seq_file} or $$runstate{blast_db} or $$runstate{blast_db_remote})) {
	$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"You must specify either both genomic sequence and database, or none.");
	$doBlast=0;
    } else {
	# BLAST not requested.
	$doBlast=0;
    }

    my $error = 0;
    if($doBlast) {
	if($$runstate{blastx}) {
	    $program = "blastx";
	} elsif($$runstate{blastn}) {
	    $program = "blastn";
	} elsif ($$runstate{blastp}) {
	    $program = "blastp";
	} 

	if($$runstate{blast_db_remote}) {
	    # qblastsubmission -- now with choice of hitlist, and with optional pvalue filter
	    if($run_win != 0) {
		my $title = "A GUI - attempting qblast (";
		$title .= $program.")";	    
		$run_win->title($title);
		$run_win->update();
	    }
	    qblast($main,$canvas,$sheet,$annotation,$seqName,$seq,$status_text,$program,1); # final 1 for full seq

	    return; # qblast will call blast_parse itself
	} elsif($$runstate{blast_annotation}) { # write temp-file
	    $DEBUG && print "DEBUG: blast annotation with annotation id $$runstate{blast_annotationid} against default $program db.\n";
	    my $annotation_nr = annotation_what_nr($annotation, $$runstate{blast_annotationid});
	    if (defined $$annotation{frame}->[$annotation_nr]) {
		$DEBUG && print "DEBUG: frame is $$annotation{frame}->[$annotation_nr] for annotation $$runstate{blast_annotationid}.\n";
		$start = $$annotation{start}->[$annotation_nr];
		$stop = $$annotation{stop}->[$annotation_nr];
		if($$annotation{frame}->[$annotation_nr] < 3) {
		    $subseq = substr($$seq, $start-1, $stop-$start+1);
		} elsif($$annotation{frame}->[$annotation_nr] < 6) {
		    # revcomp
		    $subseq = substr($$seq, $start-1, $stop-$start+1);
		    $DEVEL && $DEBUG && print "DEBUG: subseq is $subseq.\n";
		    $subseq =~ tr/ATGCatgc/TACGtacg/;
		    $subseq = join('',reverse(split(/ */, $subseq)));
		}
		if($program eq "blastp") {
		    $subseq = nt2aa($subseq);
		    
		    my $plus = chop($subseq);
		    ($plus eq "+") || ($WARNING && print "WARNING: Last character of conceptual protein translation was not plus (rather $plus).\n");
		}
		$DEBUG && print "DEBUG: start is $start, stop at $stop (length ".($stop-$start+1).") so subseq is $subseq.\n";
	    } elsif($program eq "blastp") { # frame undef and blastp selected
		$WARNING && print "WARNING: run blastp selected for a non-frame annotation box. This is not supposed to happen.\n";
		return -1;
	    } else {
		$subseq = substr($$seq, $start-1, $stop-$start+1);

		$DEBUG && print "DEBUG: No frame defined. Start is $start, stop at $stop.\n";
	    }
	} elsif($$runstate{blast_currentseq}) {
	    $subseq = $$seq;
#	} elsif($$runstate{blast_db} or $$runstate{blast_db_seq}) {
#	    $subseq = $$seq;
	}
	if($run_win != 0) {
	    $run_win->title('A GUI - writing sequence file');
	    $run_win->update();
	}
	if($$runstate{blast_currentseq} || $$runstate{blast_annotation}) {
	    open SEQOUT,">$tmppath/temp.$$seqName.$session_rand.fasta";
	    print SEQOUT ">$$seqName\n$subseq\n";
	    close SEQOUT;
	    $$runstate{blast_seq_file}="$tmppath/temp.$$seqName.$session_rand.fasta";
	}

	if($run_win != 0) {
	    $run_win->title('A GUI - building blastdb...');
	    $run_win->update();
	}

	if($$runstate{blast_db_seq}) {

	    if($$runstate{blast_wu}) {
		system("$blast2path/pressdb \"$$runstate{blast_db_seq_file}\"") == 0 or ($error=1);
	    } else { # default to ncbi
		system("$ncbiblastpath/formatdb -i \"$$runstate{blast_db_seq_file}\" -o T") == 0 or ($error=1);
	    }
	    
	    $$runstate{blast_db_file} = $$runstate{blast_db_seq_file};
	}

	($run_win != 0) && $run_win->title('Please wait - running blast search...');
	# $$runstate{blast_seq} or $$runstate{blast_currentseq}) or ($$runstate{blast_db_seq_file} or $$runstate{blast_db} or $$runstate{blast_db_remote}
	if($$runstate{blast_db_remote}) {
	    # nop -- q-blast is already done at this stage!
	} else {
	    if($$runstate{blast_wu}) {
		#$DEBUG && print "DEBUG: System $binpath/pmestLean.pl \"$$runstate{estorf_seq_file}\" \"$$runstate{estorf_est_file}\" > \"$tmppath/$$seqName.$session_rand.pmest\"\n";
		#system("$blast2path/ -i \"$$runstate{blast_db_seq_file}\" -o T") == 0 or ($error=1);
		my $filt;
		if($$runstate{blastcomplexityfilter}) {
		    $filt = "-filter dust";
		} else {
		    $filt = "";
		}

		$DEBUG && print "DEBUG: export BLASTMAT=$$sheet{blastmatrixdir}&& $blast2path/$program $$sheet{blastlocaldbdir}/".$$sheet{$program."defaultdb"}." \"$$runstate{blast_seq_file}\" $filt > \"$session_rand.blastres\" &\n";
		system("export BLASTMAT=$$sheet{blastmatrixdir} && $blast2path/$program $$sheet{blastlocaldbdir}/".$$sheet{$program."defaultdb"}." \"$$runstate{blast_seq_file}\" $filt > \"$session_rand.blastres\" &");
	    } else {
		# defaulting to ncbiblast
		my $filt;
		if($$runstate{blastcomplexityfilter}) {
		    $filt = "";
		} else {
		    $filt = "-F F";
		}
		$DEBUG && print "DEBUG: $ncbiblastpath/blastall -p $program -d $$sheet{blastlocaldbdir}/".$$sheet{$program."defaultdb"}." -i \"$$runstate{blast_seq_file}\" -o \"$session_rand.blastres\" -T T $filt &\n";
		system("$ncbiblastpath/blastall -p $program -d $$sheet{blastlocaldbdir}/".$$sheet{$program."defaultdb"}." -i \"$$runstate{blast_s.eq_file}\" -o \"$session_rand.blastres\" -T T $filt &") == 0 or ($error=1);
	    }

	    if($error != 0) {
		$$main->messageBox(-icon=>'error',-type=>'OK',-message=>"An error occured when attempting a blast run. Temporary results-file not parsed.");
		return;
	    }
	    
	    if(! $$runstate{blast_hitlist}) {
		$$sheet{blastpopup} = 0; 		# autoadd as temporaries
	    } elsif ($$runstate{blast_hitlist}) {
		$$sheet{blastpopup} = 1;
	    }
	    
	    # need non-blocking mechanism to handle multiple searches!
	    my $resultsfile = "$session_rand.blastres";
	    my $firstwait = 500; # ms; first wait often terminates on error...
	    local_blast_wait($main, $canvas, $status_text, $runstate, $annotation, $seqName, $seq, $sheet, $resultsfile, 400);
	}
    }
}

sub local_blast_wait {
    my $main=shift;
    my $canvas=shift;
    my $status_text=shift;
    my $runstate=shift;
    my $annotation=shift;
    my $seqName=shift; # ref
    my $seq=shift; # ref
    my $sheet=shift;
    my $resultsfile = shift;
    my $lastwait = shift;
      
    $DEVEL && $DEBUG && print STDERR "DEVEL DEBUG: checking local bast results file $resultsfile...\n";

    # read in the temporary results file
    my $blastresults = "";
    my $htmlend = 0;
    my $ok = open BLASTFILE, $resultsfile;
    if($ok) { # could be too early! The results file may not have been created by the client program at this stage!
	use Fcntl ':flock';
	flock(BLASTFILE, LOCK_SH|LOCK_NB);
	while(<BLASTFILE>) {
	    $blastresults .= $_;	    
	    if(m/<\/HTML>/) {
		$htmlend = 1;
	    }
	}
	close BLASTFILE;	
    }
    $blastresults eq "" && ($ok = 0); # what about partially complete files? Would they ever exist? Nah, probably not given slightly sloppy blast implementations.. =) 
    # DUH! They DO exist! Searching.... will be in there for a long time in ncbi-results. Added html-end flag since we are parsing html results.
    
    if($ok && $htmlend) {
      main::import_blast($main,$blastresults,$annotation,$seqName,$seq,$sheet,$$runstate{blast_annotationid});
    } else {
	$lastwait += 500; # increment waiting time w/ 500ms
	$$main->after($lastwait, sub { local_blast_wait($main, $canvas, $status_text, $runstate, $annotation, $seqName, $seq, $sheet, $resultsfile, $lastwait) });
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
	chomp;
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

sub getFastaQual {
    my $nQual=0;
    
    my $fastaQualFileName=shift;
    
    my $qual=shift; # ref
    open(QUALFILE, "<$fastaQualFileName") || die "Quality fasta input file $fastaQualFileName open failed\n";
    
#  print "Reading quality values...\n";
    while(<QUALFILE>) {
	chomp;
	if(/^\>/) {
	    # On description line.. Name field should be pretty much equal to the fasta file, so ignore it.
	    $nQual++;
	    $$qual[$nQual-1]="";
	} else {
	    $$qual[$nQual-1].=join(',',split(/\s+/));
	    $$qual[$nQual-1].=",";
	}
    }
    
    # Done processing quality file
    close(QUALFILE);    
}

sub select_sequence{
    my $thatName=shift;
    my $seq=shift;
    my $seqName=shift;		# get seq references
    my $nSeq=shift(@_);
    my @Seq=@_[0..($nSeq-1)];
    my @Name=@_[($nSeq)..(2*$nSeq-1)];
    
    # What number was that named sequence?
    
    my $seqNr=0;
    
    foreach my $thisName (@Name) {
	# if($thisName=~m/^$thatName/) {
	if($thisName eq $thatName) {
	    last; 
	}  
	$seqNr++;
    }
    
    # Was it found?
    if($seqNr>=$nSeq) {
	print "Warning: sequence $thatName not found.\n";
	return 0;
    } 
    
    #    $$seqName=$Name[$seqNr];
    $$seqName=$Name[$seqNr];
    print "select_sequence changed seqName ($seqName) to $$seqName\n";
    $$seq=$Seq[$seqNr];
    return 1;
} 

sub open_sequence_fasta_file {
    my $main=shift;
    my $seqFileName=shift;
    my $seqName=shift;		# get seq references
    my $seq=shift;
    
    $$seqName="DEBUG this is about to change..";
    
    # First, get all the EST sequences
#  print "DEBUG: Reading sequences from $seqFileName.\n";

    my @fileSeq;
    my @fileSeqName;
    getFastaSeq($seqFileName, \@fileSeq, \@fileSeqName);
    my $nSeq=@fileSeq;
    
 #   main::add_last_opened_file($sheet, $selected_file);

    # If there is only one sequence, update sequence and return success.
    if($nSeq==1) { 
	$$seq=$fileSeq[0];
	$$seqName=$fileSeqName[0];
	return 1;
    }
    
    # If there are several sequences, show names in sequence listbox...
    my $select_sequence_win=$$main->Toplevel;
    $select_sequence_win->title('A GUI - Select sequence');
    #    $select_sequence_win->geometry('+0+35');
    # $select_sequence_win->OnDestroy(\sequence_selected);
    my $select_sequence_win_frame=$select_sequence_win->Frame->pack(-expand => 'yes', -fill => 'both',-side=>'top');
    
    my $sequence_listbox=$select_sequence_win_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'single'); 
    $sequence_listbox->pack(qw/-side left -expand yes -fill both/);
    foreach my $nameInFile (sort @fileSeqName) {
	$sequence_listbox->insert(0,$nameInFile);
    }
    my $sequence_listbox_sb=$select_sequence_win_frame->Scrollbar(-command => ['yview', $sequence_listbox])->pack(-side=>'right',-fill=>'y');
    $sequence_listbox->configure(-yscrollcommand => ['set', $sequence_listbox_sb]);
    
    my $selectOk=0;
    
    my $sequence_select_button_frame=$select_sequence_win->Frame->pack(-side=>'bottom',-fill=>'x');  
#  print "DEBUG: Outside button command nSeq=$nSeq.\n";
    my $sequence_select_button=$sequence_select_button_frame->Button(-text=>'Select', -command=>sub {
	# Check for sequence by that name, and update sequence..
	# print "DEBUG: Inside button command nSeq=$nSeq\n";
	if(select_sequence($sequence_listbox->get('active'),$seq,$seqName,$nSeq,@fileSeq,@fileSeqName)) {
	    $selectOk=1;
	    # print "seqName ($seqName) inside button command is now $$seqName.\n";
	    $select_sequence_win->destroy;    
	}; } )->pack(-side=>'bottom');
    
    print "DEBUG: Waiting for select_seq_win destruction...\n";
    # wait until window destruction...
    $select_sequence_win->waitWindow;
    
    print "DEBUG: Selected sequence is: ",$$seqName,"\n";
    
    print "DEBUG: select_seq_win destroyed! Return status $selectOk.\n";
    return $selectOk;
}


sub heur_programme {
    my $main = shift;
    my $annotation = shift;
    my $seqName = shift;
    my $seq = shift;
    my $sheet = shift;

    my %runstate = ();
    my %merge_status = ();

# Experimental function!

    # add orf   
    $$sheet{heur_prog_orf} && add_orf($main,$annotation,$seqName,$seq,$sheet);

    # run testcode -- default options, current sequence
    $$sheet{heur_prog_testcode} && ($runstate{testcode_current_seq} = 1);
    
    # run glimmer
    if(defined($$sheet{glimdefaultmodel})) {
	$runstate{glim_run} = 1;
	
	$runstate{glim_model} = 1;
	$runstate{glim_model_file} = $$sheet{glimdefaultmodel};
	$runstate{glim_shortest_orf} = 100;
    }

    # run estorf
    if(defined($$sheet{estorfdefaultestfile})) {
	$runstate{estorf_currentseq} = 1;
	$runstate{estorf_est} = 1;
	$runstate{estorf_est_file} = $$sheet{estorfdefaultestfile};
    }

    # * run full blastx w/ fairly strict thresh

    # actually run those externals.
    my $run_win = $$main->Toplevel;
    my $run_win_label = $run_win->Label(-text=>"Running requested external analysis programs. Please wait.");
    ($runstate{glim_run} || $runstate{estorf_currentseq} || $runstate{testcode_current_seq}) && run_external($main,$run_win,\%runstate,$annotation,$seqName,$seq,$sheet);

    # run polypy
    $$sheet{heur_prog_polypy} && add_polypy($main, $annotation, $seqName, $seq, $sheet);
    
    # * run full repeat check

    # remove short orfs from all (say < 150 to begin with)

    # (?) strandclean -- olap handling??
    # instead of a strandclean, one might chose to remove glimmer, testcode & estorf whith reasonably strict criteria
    # merge

    %merge_status = (glimmer2 => 1, testcode => 1, ESTORF => 1, comp => 1, shortest_olap =>100, framex => 1);
    merge($main,$canvas,$annotation,\%merge_status,$sheet,$seqName);

    # %merge_status = (blast => 1, manual => 1, or => 1);
    # merge($main,$canvas,$annotation,\%merge_status,$sheet,$seqName);

    # but this will leave us with a shitload of overlapping "genes" and otw unlikely phenomena
    # --> extend the strand cleaning routine w/ knowledge of glimmer/testcode/blasthit

    $DEBUG && print "DEBUG: done heuristic annotation programme.\n";
}




sub nt2aa {
    my $tri=$_[0];
    my $aa="";
    # Make nucleotides lowercase, and have all AAs uppercase.
    $tri=~tr/ATCGNX/atcgnx/;

    # DEBUG: Strict??
    # print $DEBUG && "DEBUG: Length ", length($tri), "and the third is ", length($tri)/3 ,".\n";
    for(my $k=0;$k<(length($tri)/3);$k++) {
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
	s/[atgc]*[nx]+[atgc]*[nx]*/U/g; # only substrings of three at the time are processed
	
	# Append the aa the codon encoded to the virtual protein.
	$aa.=$_;
    }
    # Remove the extra nucleotides...

    $aa=~s/[atgcnx]//g;
    return $aa;
}

sub dinucleotide_descriptor {
    my $seq = shift;
    my $dinuc = shift; # BS, PT, PD
    my $winsize = shift;
    my $xAxis = shift;

    @$xAxis = ();

    my $windows = POSIX::floor( (length($$seq) - 2) / $winsize );

    my $pos;
    for(my $j=0;$j < $windows; $j++) {
	$pos = $winsize*$j;
	push @$xAxis, $pos;

	for(my $i=0; $i < $winsize; $i++) {
	    $pos = $i + $winsize*$j;
	    $_ = substr($$seq, $pos, 2);
	    tr/atgc/ATGC/;

	  SWITCH: {
	      (/AA/ || /TT/) && do {
		  $$dinuc{BS}->[$j] += -5.37;
		  $$dinuc{PT}->[$j] += -18.66;
		  $$dinuc{PD}->[$j] += 2.9;
		  last SWITCH;
	      };
	      (/AC/ || /GT/) && do {
		  $$dinuc{BS}->[$j] += -10.51;
		  $$dinuc{PT}->[$j] += -13.10;
		  $$dinuc{PD}->[$j] += 2.3;
		  last SWITCH;
	      };
	      (/AG/ || /CT/) && do {
		  $$dinuc{BS}->[$j] += -6.78;
		  $$dinuc{PT}->[$j] += -14;
		  $$dinuc{PD}->[$j] += 2.1;
		  last SWITCH;
	      };
	      /AT/ && do {
		  $$dinuc{BS}->[$j] += -6.57;
		  $$dinuc{PT}->[$j] += -15.01;
		  $$dinuc{PD}->[$j] += 1.6;
		  last SWITCH;
	      };
	      (/CA/ || /TG/) && do {
		  $$dinuc{BS}->[$j] += -6.57;
		  $$dinuc{PT}->[$j] += -9.45;
		  $$dinuc{PD}->[$j] += 9.8;
		  last SWITCH;
	      };
	      (/CC/ || /GG/) && do {
		  $$dinuc{BS}->[$j] += -8.26;
		  $$dinuc{PT}->[$j] += -8.11;
		  $$dinuc{PD}->[$j] += 6.1;
		  last SWITCH;
	      };
	      /CG/ && do {
		  $$dinuc{BS}->[$j] += -9.69;
		  $$dinuc{PT}->[$j] += -10.03;
		  $$dinuc{PD}->[$j] += 12.1;
		  last SWITCH;
	      };
	      (/GA/ || /TC/) && do {
		  $$dinuc{BS}->[$j] += -9.81;
		  $$dinuc{PT}->[$j] += -13.48;
		  $$dinuc{PD}->[$j] += 4.5;
		  last SWITCH;
	      };
	      /GC/ && do {
		  $$dinuc{BS}->[$j] += -14.59;
		  $$dinuc{PT}->[$j] += -11.08;
		  $$dinuc{PD}->[$j] += 4;
		  last SWITCH;
	      };
	      /TA/ && do {
		  $$dinuc{BS}->[$j] += -3.82;
		  $$dinuc{PT}->[$j] += -11.85;
		  $$dinuc{PD}->[$j] += 6.3;
		  last SWITCH;
	      };
	      
	      # N, X
	      $$dinuc{BS}->[$j] += 0;
	      $$dinuc{PT}->[$j] += 0;
	      $$dinuc{PD}->[$j] += 0;	    
	  }
	    #$DEVEL && print "DEVEL: Substr $j is $_ and has BS = $$dinuc{BS}->[$j], PT = $$dinuc{PT}->[$j], PD = $$dinuc{PD}->[$j].\n";	
	    
	    # lazy normalize to 0 < x < 1
	    $$dinuc{BS}->[$j] /= -14.59;
	    $$dinuc{PT}->[$j] /= -18.66;
	    $$dinuc{PD}->[$j] /= 12.1;
	    
	    #$DEVEL && print "DEVEL: After a lazy normalisation, BS = $$dinuc{BS}->[$j], PT = $$dinuc{PT}->[$j], PD = $$dinuc{PD}->[$j].\n";
	}
    }
}


sub trinucleotide_descriptor {
    my $seq = shift;
    my $trinuc = shift; # B, PP
    my $winsize = shift;
    my $xAxis = shift;

    my $windows = POSIX::floor( (length($$seq) - 3)/$winsize );

    my $pos;

    push @$xAxis, $pos;

    for(my $i=0;$i < $windows; $i++) {
	$pos = $winsize*$i;
	push @$xAxis, $pos;
	
	for(my $j=0;$j < $winsize ; $j++) {
	    $pos = $winsize*$i + $j;

	    $_ = substr($$seq, $pos, 3);
	    tr/atgc/ATGC/;
	    
	  SWITCH: {
	      (/AAA/ || /TTT/) && do {
		  $$trinuc{B}->[$i] += -0.274;
		  $$trinuc{PP}->[$i] += 36;
		  last SWITCH;
	      };
	      (/AAC/ || /GTT/) && do {		
		  $$trinuc{B}->[$i] += -0.205;
		  $$trinuc{PP}->[$i] += 6;
		  last SWITCH;
	      };
	      (/AAG/ || /CTT/) && do {	
		  $$trinuc{B}->[$i] += -0.081;
		  $$trinuc{PP}->[$i] += 6;
		  last SWITCH;
	      };
	      (/AAT/ || /ATT/) && do {	
		  $$trinuc{B}->[$i] += -0.280;
		  $$trinuc{PP}->[$i] += 30;
		  last SWITCH;
	      };
	      (/ACA/ || /TGT/) && do {	
		  $$trinuc{B}->[$i] += -0.006;
		  $$trinuc{PP}->[$i] += 6;
		  last SWITCH;
	      };
	      (/ACC/ || /GGT/) && do {	
		  $$trinuc{B}->[$i] += -0.032;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      (/ACG/ || /CGT/) && do {	
		  $$trinuc{B}->[$i] += -0.033;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      (/ACT/ || /AGT/) && do {
		  $$trinuc{B}->[$i] += -0.183;
		  $$trinuc{PP}->[$i] += 11;
		  last SWITCH;
	      };
	      (/AGA/ || /TCT/) && do {
		  $$trinuc{B}->[$i] += 0.027;
		  $$trinuc{PP}->[$i] += 9;
		  last SWITCH;
	      };
	      (/AGC/ || /GCT/) && do {
		  $$trinuc{B}->[$i] += 0.017;
		  $$trinuc{PP}->[$i] += 25;
		  last SWITCH;
	      };
	      (/AGG/ || /CCT/) && do {	
		  $$trinuc{B}->[$i] += -0.057;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      (/ATA/ || /TAT/) && do {
		  $$trinuc{B}->[$i] += 0.182;
		  $$trinuc{PP}->[$i] += 13;
		  last SWITCH;
	      };
	      (/ATC/ || /GAT/) && do {
		  $$trinuc{B}->[$i] += -0.11;
		  $$trinuc{PP}->[$i] += 7;
		  last SWITCH;
	      };
	      (/ATG/ || /CAT/) && do {
		  $$trinuc{B}->[$i] += 0.134;
		  $$trinuc{PP}->[$i] += 18;
		  last SWITCH;
	      };
	      (/CAA/ || /TTG/) && do {
		  $$trinuc{B}->[$i] += 0.015;
		  $$trinuc{PP}->[$i] += 9;
		  last SWITCH;
	      };
	      (/CAC/ || /GTG/) && do {
		  $$trinuc{B}->[$i] += 0.040;
		  $$trinuc{PP}->[$i] += 17;
		  last SWITCH;
	      };
	      (/CAG/ || /CTG/) && do {
		  $$trinuc{B}->[$i] += 0.175;
		  $$trinuc{PP}->[$i] += 2;
		  last SWITCH;
	      };
	      (/CCA/ || /TGG/) && do {
		  $$trinuc{B}->[$i] += -0.246;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      (/CCC/ || /GGG/) && do {
		  $$trinuc{B}->[$i] += -0.012;
		  $$trinuc{PP}->[$i] += 13;
		  last SWITCH;
	      };
	      (/CCG/ || /CGG/) && do {
		  $$trinuc{B}->[$i] += -0.136;
		  $$trinuc{PP}->[$i] += 2;
		  last SWITCH;
	      };
	      (/CGA/ || /TCG/) && do {
		  $$trinuc{B}->[$i] += -0.003;
		  $$trinuc{PP}->[$i] += 31;
		  last SWITCH;
	      };
	      (/CGC/ || /GCG/) && do {
		  $$trinuc{B}->[$i] += -0.077;
		  $$trinuc{PP}->[$i] += 25;
		  last SWITCH;
	      };
	      (/CTA/ || /TAG/) && do {
		  $$trinuc{B}->[$i] += 0.09;
		  $$trinuc{PP}->[$i] += 18;
		  last SWITCH;
	      };
	      (/CTC/ || /GAG/) && do {
		  $$trinuc{B}->[$i] += 0.031;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      (/GAA/ || /TTC/) && do {
		  $$trinuc{B}->[$i] += -0.037;
		  $$trinuc{PP}->[$i] += 12;
		  last SWITCH;
	      };
	      (/GAC/ || /GTC/) && do {
		  $$trinuc{B}->[$i] += -0.013;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      /GCA/ || /TGC/ && do {
		  $$trinuc{B}->[$i] += 0.076;
		  $$trinuc{PP}->[$i] += 13;
		  last SWITCH;
	      };
	      (/GCC/ || /GGC/) && do {
		  $$trinuc{B}->[$i] += 0.107;
		  $$trinuc{PP}->[$i] += 45;
		  last SWITCH;
	      };
	      (/GGA/ || /TCC/) && do {
		  $$trinuc{B}->[$i] += 0.013;
		  $$trinuc{PP}->[$i] += 5;
		  last SWITCH;
	      };
	      (/GTA/ || /TAC/) && do {
		  $$trinuc{B}->[$i] += 0.025;
		  $$trinuc{PP}->[$i] += 6;
		  last SWITCH;
	      };
	      (/TAA/ || /TTA/) && do {
		  $$trinuc{B}->[$i] += 0.068;
		  $$trinuc{PP}->[$i] += 20;
		  last SWITCH;
	      };
	      (/TCA/ || /TGA/) && do {
		  $$trinuc{B}->[$i] += 0.194;
		  $$trinuc{PP}->[$i] += 8;
		  last SWITCH;
	      };
	      
	      $$trinuc{B}->[$i] += 0;
	      $$trinuc{PP}->[$i] += 0;

	  }
	    #$DEVEL && print "DEVEL: Substr $i is $_ and has B = $$trinuc{B}->[$i], PP = $$trinuc{PP}->[$i].\n";
	    
	    # lazy normalize..
	    $$trinuc{B}->[$i] =  $$trinuc{B}->[$i] / (2 * -0.280) + 0.5;
	    $$trinuc{PP}->[$i] /= 45;
	    #$DEVEL && print "DEVEL: After a lazy noramlistion, B = $$trinuc{B}->[$i], PP = $$trinuc{PP}->[$i].\n";
	}
    }
}

sub log2 {
    my $x = shift;
    return log($x)/log(2);
}

sub maxref {
   my $x = shift; # array ref
   my $max=0;
   foreach my $testme (@{$x}) {
       if($testme>$max) {
	   $max=$testme;
       }
   }
   return $max;  
}

sub minref {

    my $x = shift; # array ref

    my $min=$$x[0];

    foreach my $testme (@{$x}) {
	if($testme<$min) {
	    $min=$testme;
	}
    }
    return $min;

}

sub max {
    my @x = @_;
    my $max=0;
    foreach my $testme (@x) {
	if($testme>$max) {
	    $max=$testme;
	}
    }
    return $max;
}

sub min {
    my @x = @_;
    my $min=$x[0];
    foreach my $testme (@x) {
	if($testme<$min) {
	    $min=$testme;
	}
    }
    return $min;
}

sub sum {
    my $sum=0;
    for(my $n=0;$n<@_;$n++) {
	$sum += $_[$n];
    }
    return $sum;
}

sub logo {
    print "-oO0O---O0Oo-\n  | |   | |  \n  | |\"\"\"| |  \n  | |O O| |  \n  | | ^ | |  \n   \\ \\0/ /   \n(c)2000-03 DN\n";
}

# print "Duh?\n";
return 1;
