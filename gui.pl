#!/usr/bin/perl -W 
#-I/home/daniel/work/gui
#
# Experimental GUI for use in genefinding project.
# Expanded to meet Trypanosoma cruzi genome project annotation needs
#
# Daniel Nilsson, 000119
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
# A hint on module paths:
#
# Add a  -I <installation directory> to the first line above if use lib fails, 
# or copy the pm files to your perl lib path, or add the installation directory
# to your lib path or change the "pathToGuiIO" variable a few lines below. 
#

#use strict;
#use strict 'refs';

# Modules
use Tk;
#use Tk::FileSelect;
#use Tk::FBox;
#use Tk::NoteBook;
#use Tk::Dialog;
#use Tk::Popup;

use POSIX qw(ceil);

our $DEBUG;
our $WARNING;
our $DEVEL;
my %sheet;

my @clargs;

BEGIN {
    # Global verbosity flags default
    $DEBUG=0;
    $WARNING=1;
    $DEVEL=0;    

    # Parse command line parameters
    while (my $arg = shift @ARGV) {
	if($arg eq "--debug") {
	    $DEBUG = 1;
	} elsif ($arg eq "--no-warn") {	
	    $WARNING = 0;
	} elsif ($arg eq "--devel") {
	    $DEVEL = 1;
	} else {
	    $DEVEL && $WARNING && print STDERR "WARNING: Unrecognised command line parameter $arg. Continuing anyway.\n";
	    push @clargs, $arg;
	}
    }

    # Get the libpath up before libraries are included

    ######## Change this to point to the installation directory OR copy the *pm files to somewhere on your @INC path. Obviously, you could also change your @INC path to include the installation directory.
    $sheet{pathToGuiIO} = '/usr/remote/bin';
    ########
    
    $DEVEL && ($sheet{pathToGuiIO} = '/home/daniel/work/gui');

    # Apperance - colors
    $sheet{menu_bar_color} = 'midnightblue';
    $sheet{menu_button_bg_color} = 'midnightblue';
    $sheet{menu_button_active_bg_color} = 'cornflowerblue';
    $sheet{menu_button_foreground_color}= 'white';

    $sheet{default_win_background} = 'linen';
   #  $main::main_win_background=$sheet{default_win_background};
    $sheet{status_text_col}= 'black';

    $sheet{dialog_button_active_background} = 'AntiqueWhite2';

    $sheet{faint_guide_color} = 'antiquewhite';
#    $sheet{faint_guide_color} = 'AntiqueWhite3';
}

use lib (@INC,$sheet{pathToGuiIO});
use GuiAnnotator; # Annotations display configuration and console window data dumps
use GuiIO; # Disk IO, network IO and some other helper functions

# Declare subroutines

# canvas interaction
sub click_sequence;
sub middle_click_sequence;
sub click_annotation;
sub middle_click_annotation;
sub right_click_annotation;
sub right_click_guide;

# menu interaction
sub post_lastfiles_menu;
sub open_last_opened;

sub file_new;
sub file_export_sequin;
sub file_export_gff;
sub file_export_genbank;

sub help_about;
sub help_help;

sub annotation_remove;
sub annotation_removeAll;
sub annotation_merge;
sub annotation_edit;

sub view_zoomIn;
sub view_zoomOut;
sub view_center_on;
sub view_graph;
sub view_legend;
sub view_levels;
sub view_stats;

# annotation helper subs
sub import_glimmer;
sub import_testcode;
sub import_estorf;
sub import_est;
sub import_blast;
sub import_gff;

sub level_layout;
sub level_layout_gff;
sub level_layout_all;

sub add_est;
sub add_estorf;
sub add_glimmer;
sub add_testcode;
sub add_blast;
sub add_manual;
sub add_regexp;
sub add_orf;
sub add_polypy;
sub add_st;
sub add_gff;

sub select_levels;

sub annotation_add_new;

sub fill_in_remove_fields;
sub check_selection;
sub remove_all;
sub remove_annotations;
sub redraw_annotations;

sub annotation_id_type_current;
sub annotation_id_type_selected;

sub annotation_what_selected_nr;
sub update_annotation_nr_cache;

# computations
sub calculate_graph;

sub heuristic_likely_orf;

# run menu
sub run;

# canvas graphics
sub plot_axis;
sub plot_curve;

sub set_main_view_height;

sub draw_seq;
sub draw_annotation_arrow;
sub draw_annotation_box;
sub draw_annotation_st;

# some global variables

# my $AnnotationBoxFrame = 6;

my $mainViewMaxHeight = 47;
my $mainViewHeight = $mainViewMaxHeight;    # initially; later updated dynamically with values from level_selection

# Sheet variables

# As this program expands, we'll probably not be too interested in sending tons of status 
# variables to and fro between different subsystems. Treat this as a worksheet object..

$sheet{canvas_seq_height} = 20 * $mainViewHeight + 40;    # Canvas "layout" - sequence field
$sheet{canvas_plot_height} = 100;   # and a plot field below, to show when needed.
$sheet{status} = "closed";	# {closed,clean,dirty}
$sheet{zoom} = "normal"; # { birdseye,overview,normal,sequence }
$sheet{selected} = ""; # seleceted annotation_id
$sheet{seq_sel_start} = 1;
$sheet{seq_sel_stop} = 1;
$sheet{file} = "";

# Read preferences

# Default preferences
$sheet{tmppath} = "/tmp";
$sheet{binpath} = "/home/daniel/work/gui";
$sheet{blast2path} = "/usr/remote/blast2";
$sheet{ncbiblastpath} = "/usr/remote/ncbiBlast21";

$sheet{usebuiltintestcode} = 1;
$sheet{testcodepath} = "/usr/remote/bin";

$sheet{glimmer2altstartpath}="/home/daniel/install/Glimmer2.0-altstart";
$sheet{glimaltstartbin}=$sheet{binpath}."/glimAltStart.pl";

$sheet{glimmer2path}="/home/daniel/install/Glimmer2.0";
$sheet{glimbin}=$sheet{binpath}."/glim.pl";

$sheet{uselastfiles} = 0;
$sheet{nolastfiles} = 0;

foreach my $key (keys %annotatorName) {
    $sheet{"display_${key}_level"} = 1;
}

# systemwide
my $prefs_file = "/usr/remote/etc/aguirc";
read_prefs($prefs_file, \%sheet);

# override with user settings
#$prefs_file = "~/.aguirc";
my $user = `whoami`;
chop $user;
$prefs_file = "/home/$user/.aguirc";
read_prefs($prefs_file, \%sheet);

# recalc those levels..
set_main_view_height(\%sheet);

# Sequence
our $seqName; # GLOBAL seqName
our $seq; # GLOBAL seq

my %annotation;			# $nr=$annotation{nr,unr}; $annotation{seqName, type, start,stop,frame, level, id,uid}->[$nr];
$annotation{unr} = 0; # reset the "unique" id counter
our %annotation_nr_cache = (); # GLOBAL nr cache, indexed on annotation id

# State variables for plot
my %plotstatus = (axis=>0,GC=>0,AT=>0,AG=>0,CT=>0,H=>0,AMI=>0,BS=>0,PT=>0,PD=>0,B=>0,PP=>0,windowSize=>100); # { axis, GC, AT, AG, CT, H, AMI, BS, PT,PD, B, PP, windowSize }

# Status field var

my $status_text="GUI ready to rock.";

# Main window
$DEBUG && print STDERR "main..\n"; 
my $main = MainWindow->new(); 
$main->title("A GUI");
$main->configure(-background=>$sheet{default_win_background});
# $main->geometry('+150+0');
#$main->minsize( qw(250 

my $canvas;			# Will be needed..
my $scrollbar;

# Menu
my $menu_bar = $main->Frame(-relief=>'groove', -borderwidth=>3, -background=>$sheet{menu_bar_color})->pack(-'side'=>'top',-fill=>'x');

# Predefine menubuttons for use in file_mb subs.
my $view_mb;
my $annotation_mb;
my $run_mb;

# Menu - File
my %mbcol = ( '-background', $sheet{menu_button_bg_color}, '-activebackground', $sheet{menu_button_active_bg_color},'-foreground', $sheet{menu_button_foreground_color});

my $file_mb=$menu_bar->Menubutton(-text=>'File',-underline=>0,%mbcol)->pack(-side=>'left');
my $file_m=$file_mb->Menu(-tearoff => 0, %mbcol);
$file_mb->configure(-menu => $file_m);

my $file_new_label= 'New...';
my $file_new_label_ul = 0;
my $file_load_label = 'Load...';
my $file_load_label_ul = 0;
my $file_open_genbank_label = 'Open genbank...';
my $file_open_genbank_label_ul = 0;
my $file_save_label = 'Save';
my $file_save_label_ul = 0;
my $file_save_as_label = 'Save As...';
my $file_save_as_label_ul = 5;
my $file_export_label = 'Export to sequin...';
my $file_export_label_ul = 0;
my $file_export_gff_label = 'Export GFF...';
my $file_export_gff_label_ul = 7;
my $file_export_genbank_label = 'Export genbank...';
my $file_export_genbank_label_ul = 10;

sub disable_menu_state {
    $view_mb->configure(-state=>"disabled");
    $annotation_mb->configure(-state=>"disabled");
    $run_mb->configure(-state=>"disabled");

    $file_m->entryconfigure($file_save_label, -state=>'disabled');
    $file_m->entryconfigure($file_save_as_label, -state=>'disabled');
    $file_m->entryconfigure($file_export_genbank_label, -state=>'disabled');
    $file_m->entryconfigure($file_export_gff_label, -state=>'disabled');
    $file_m->entryconfigure($file_export_label, -state=>'disabled');
}

sub reset_menu_state {

    $view_mb->configure(-state=>"normal");
    $annotation_mb->configure(-state=>"normal");
    $run_mb->configure(-state=>"normal");

    $file_m->entryconfigure($file_export_genbank_label, -state=>"normal");
    $file_m->entryconfigure($file_export_gff_label, -state=>"normal");
    $file_m->entryconfigure($file_export_label, -state=>"normal");
}

$file_m->add('command', -label=>$file_new_label, -underline => $file_new_label_ul, -command=>sub{ 
		 if(file_new(\$main,\$canvas,\%sheet,\$seqName,\$seq,\%annotation,\%plotstatus)) {
		     reset_menu_state;
		     if($seqName) {
			 $main->title("A GUI - $seqName");
			 my $l=length($seq); $status_text="Opened sequence $seqName, $l bp.";
		     }
		 }
	     } );


$file_m->add('command', -label=>$file_load_label, -underline => $file_load_label_ul,
	     -command=>sub {
		 my $status=file_load(\$main,\$canvas,\%annotation,\%sheet,\%plotstatus,\$seqName,\$seq);
		 if($status) {
		     view_refresh(\$main,\$canvas,\%annotation,\$seqName,\$seq,\%plotstatus,\%sheet);
#		     draw_seq(\$canvas,$seqName,\$seq,\%sheet);
#		     redraw_annotations(\$canvas,\%annotation,\%sheet,\$seq);

		     reset_menu_state;

		     if($seqName) {
			 $main->title("A GUI - $seqName [$sheet{file}]");
			 my $l=length($seq); $status_text="Loaded sequence $seqName, $l bp.";
		     }
		 }
	     });


$file_m->add('command', -label=>$file_open_genbank_label,-underline=>$file_open_genbank_label_ul, 
	     -command=>sub {
		 $ok = file_close(\$main,\$canvas,\%sheet,\%annotation,\$seqName,\$seq,\%plotstatus);
		 if($ok) {
		     $WARNING && print "File saved ok, or discarded from file_close.\n";
		 } else {
		     $WARNING && print "Cancel new.\n";
		     return 0;			# Ok, so cancel new.
		 }

		 my $status = open_genbank(\$main, \%annotation, \$seqName, \$seq, \%sheet);
		 if($status) {
		     draw_seq(\$canvas,$seqName,\$seq,\%sheet);
		     $canvas->configure(-height=>$sheet{canvas_seq_height});
		     redraw_annotations(\$canvas,\%annotation,\%sheet,\$seq);
		     if($plotstatus{axis} && ($sheet{zoom} ne 'sequence')) {
			 $DEBUG && print "DEBUG: reintroducing graph...\n";
			 $plotstatus{axis}=0; # Tell calculate_graph to redraw axis...
			 calculate_graph(\$main,\$canvas,\$seq,\%plotstatus,\%sheet);      
		     }
		     reset_menu_state;
		     
		     if($seqName) {
			 $main->title("A GUI - $seqName [$sheet{file}]");
			 my $l=length($seq); $status_text="Loaded sequence $seqName, $l bp.";
		     }
		 } 
	     });



$file_m->add('command',-label=>$file_save_label, -underline=>$file_save_label_ul,
	     -command=>sub { 
		 file_save(\$main,\%annotation,\%sheet,\%plotstatus);
		 if($seqName && $sheet{file}) {
		     $main->title("A GUI - $seqName [$sheet{file}]");
		 }
	     });


$file_m->add('command',-label=>$file_save_as_label,-underline=>$file_save_as_label_ul,
	     -command=> sub { 
		 file_save(\$main,\%annotation,\%sheet,\%plotstatus,1);
		 if($seqName && $sheet{file}) {
		     $main->title("A GUI - $seqName [$sheet{file}]");
		 }
	     });


$file_m->add('command',-label=>$file_export_gff_label, -underline => $file_export_gff_label_ul,
	     -command=> sub {
		 file_export_gff(\$main,\$seqName,\$seq,\%annotation, \%sheet);
	     });


$file_m->add('command',-label=>$file_export_genbank_label, -underline =>$file_export_genbank_label_ul,
	     -command=> sub {
		 file_export_genbank(\$main,\$seqName,\$seq,\%annotation, \%sheet);
	     });

$file_m->add('command',-label=>$file_export_label, -underline=>$file_export_label_ul,
	     -command=> sub {
		 file_export_sequin(\$main,\$seqName,\$seq,\%annotation, \%sheet);
	     });

my $lastfilemenu;

if($sheet{uselastfiles}) {
    $file_mb->separator();
    $lastfilemenu = $file_mb->cascade(-tearoff => 0, -label=>'Open ~recent',
				      -postcommand => sub { post_lastfiles_menu(\$main,\$canvas, \%annotation, \%sheet, \%plotstatus,\$seqName,\$seq,\$lastfilemenu); } );
}

$file_mb->separator();
$file_m->add('command',-label=>'Print', -command=>sub { 
    my $ps_file_name = $main->getSaveFile(-initialfile => "out.ps",-defaultextension => "ps",-title=>"Print to file");

    if( $ps_file_name ) {
	$canvas->postscript(-file=>$ps_file_name);
    }
});
$file_mb->separator();

$file_m->add('command',-label=>'Close', -underline=>0, -command=>sub{ 
		    if(file_close(\$main,\$canvas,\%sheet,\%annotation,\$seqName,\$seq,\%plotstatus)) {
			disable_menu_state;
			$main->title("A GUI");
			$status_text="A GUI ready to rock.";
		    }
		  });
$file_mb->separator();
$file_m->add('command',-label=>'Exit', -underline=>1, -command=>sub{$main->destroy; &logo});

$file_m->configure(-postcommand=>sub { 
		     if($sheet{status} eq "dirty") {
			     $file_m->entryconfigure($file_save_label, -state=>'normal');
			     $file_m->entryconfigure($file_save_as_label, -state=>'normal');
		     } elsif ($sheet{status} eq "clean") {
			     $file_m->entryconfigure($file_save_label, -state=>'disable');
			     $file_m->entryconfigure($file_save_as_label, -state=>'normal');
		     } else {
			     $file_m->entryconfigure($file_save_label, -state=>'disabled');
			     $file_m->entryconfigure($file_save_as_label, -state=>'disabled');
		     }
		   } );

# Menu - View

$view_mb=$menu_bar->Menubutton(-text=>'View',-underline=>0)->pack(-side=>'left');
# Invoking view_m->refresh is pretty convenient, hence "our". 
our $view_m = $view_mb->Menu(-tearoff => 0,%mbcol); 
$view_mb->configure(-menu => $view_m, %mbcol);

my $view_mb_zoomIn_label = 'Zoom In';
my $view_mb_zoomIn_label_ul = 5;
my $view_mb_zoomIn=$view_m->add('command',-label=>$view_mb_zoomIn_label,-underline=>$view_mb_zoomIn_label_ul,
				-command=>sub {
				    my @old_scroll_pos = $canvas->xview(); # save scroll center position
				    my $old_scroll_center = $old_scroll_pos[0]+($old_scroll_pos[1]-$old_scroll_pos[0])/2;
				    
				    view_zoomIn(\$main,\$canvas,\%annotation,\%sheet,\%plotstatus,\$seqName,\$seq);
				    
				    my @current_scroll_pos = $canvas->xview(); # calculate new scroll width
				    my $current_scroll_width = ($current_scroll_pos[1] - $current_scroll_pos[0]);
				    
				    $canvas->xview(moveto=>$old_scroll_center-($current_scroll_width/2)); # scroll the canvas - it controls the scrollbar in turn!
				    
				    $DEVEL && print "DEVEL: Old scroll got $old_scroll_pos[0]-$old_scroll_pos[1], center $old_scroll_center, new scroll got $current_scroll_pos[0]-$current_scroll_pos[1], width $current_scroll_width so the new scroll should be set to ",$old_scroll_center-($current_scroll_width/2),", ",$old_scroll_center+($current_scroll_width/2), ".\n";
				});
my $view_mb_zoomOut_label = 'Zoom Out';
my $view_mb_zoomOut_label_ul = 5;
my $view_mb_zoomOut=$view_m->add('command',-label=>$view_mb_zoomOut_label,-underline=>$view_mb_zoomOut_label_ul,
				 -command=>sub {
				     my @old_scroll_pos = $canvas->xview(); # save scroll center position
				     my $old_scroll_center = $old_scroll_pos[0]+($old_scroll_pos[1]-$old_scroll_pos[0])/2;
				     
				     view_zoomOut(\$main,\$canvas,\%annotation,\%sheet,\%plotstatus,\$seqName,\$seq);
				     
				     my @current_scroll_pos = $canvas->xview(); # calculate new scroll width
				     my $current_scroll_width = ($current_scroll_pos[1] - $current_scroll_pos[0]);
				     
				     $canvas->xview(moveto=>$old_scroll_center-($current_scroll_width/2)); # and set old scroll center, preserving new width
				     
				     $DEVEL && print "DEVEL: Old scroll got $old_scroll_pos[0]-$old_scroll_pos[1], center $old_scroll_center, new scroll got $current_scroll_pos[0]-$current_scroll_pos[1], width $current_scroll_width so the new scroll should be set to ",$old_scroll_center-($current_scroll_width/2),", ",$old_scroll_center+($current_scroll_width/2), ".\n";
				 });
$DEVEL && ($view_m->add('command',-label=>'Scrolltest', -underline=>0, -foreground=>$sheet{menu_button_foreground_color}, -command=>sub { 
    my @current_scroll_pos = $canvas->xview();
    my $current_scroll_width = ($current_scroll_pos[1] - $current_scroll_pos[0]);   
#    $scrollbar->set(0.5-($current_scroll_width/2),0.5+($current_scroll_width/2));
    $canvas->xview(moveto=>0.5-($current_scroll_width/2));
}));
my $view_mb_zoom = $view_m->add('command',-label=>'Zoom...', -underline=>0,-command=>sub {  
    view_zoom(\$main,\$canvas,\%annotation,\%sheet,\%plotstatus,\$seqName,\$seq);});

my $view_mb_goto = $view_m->add('command',-label=>'Center on...', -underline=>0, -command=>sub { view_center_on(\$main, \$canvas,\%annotation,\$seq,\%sheet); });

my $view_mb_levels = $view_m->add('command',-label=>'Annotation levels...', -underline=>0, -command=>sub { 
    view_levels(\$main, \$canvas,\%annotation,\$seqName,\$seq,\%plotstatus,\%sheet); 
});

my $view_mb_refresh = $view_m->add('command',-label => 'Refresh',-underline=>0, -command=>sub {
    view_refresh(\$main,\$canvas,\%annotation,\$seqName,\$seq,\%plotstatus,\%sheet);
});

$view_mb->separator();
$view_mb_graph_label = 'Display Graph...';
$view_mb_graph_label_ul = 8;
$view_m->add('command',-label=>$view_mb_graph_label,-underline=>$view_mb_graph_label_ul, 
	     -command=>sub{view_graph(\$main,\$canvas,\$seq,\%plotstatus,\%sheet)});
$view_m->add('command',-label=>'Display Legend', -underline=>8, -command=>sub{ view_legend(\$main,\%annotation,\%plotstatus,\%sheet) });
$view_mb->separator();
$view_m->add('command',-label=>'Display Stats...', -underline=>8, -command=>sub{view_stats(\%annotation,\$seqName,\$seq)});

my $v_menu=$view_mb->cget("-menu");
$v_menu->configure(-postcommand=>sub {
    if($sheet{zoom} eq "birdseye") {
	$view_m->entryconfigure($view_mb_zoomIn_label, -state=>'normal');
	$view_m->entryconfigure($view_mb_zoomOut_label, -state=>'disabled');
	$view_m->entryconfigure($view_mb_graph_label, -state=>'normal');
    } elsif($sheet{zoom} eq "overview" or $sheet{zoom} eq "normal") { 
	$view_m->entryconfigure($view_mb_zoomIn_label, -state=>'normal');
	$view_m->entryconfigure($view_mb_zoomOut_label, -state=>'normal');
	$view_m->entryconfigure($view_mb_graph_label, -state=>'normal');
    } elsif($sheet{zoom} eq "sequence") {
	$view_m->entryconfigure($view_mb_zoomIn_label,-state=>'disabled');
	$view_m->entryconfigure($view_mb_zoomOut_label,-state=>'normal');
	$view_m->entryconfigure($view_mb_graph_label-state=>'disabled');
    } elsif ($sheet{zoom} eq "factor") {
	$view_m->entryconfigure($view_mb_zoomIn_label, -state=>'disabled');
	$view_m->entryconfigure($view_mb_zoomOut_label, -state=>'disabled');
	$view_m->entryconfigure($view_mb_graph_label, -state=>'normal');
    }});

# Menu - Annotation

$annotation_mb=$menu_bar->Menubutton(-text=>'Annotation',-underline=>0,%mbcol)->pack(-side=>'left');
my $annotation_m = $annotation_mb->Menu(-tearoff => 0, %mbcol);

$annotation_mb->configure(-menu => $annotation_m);
my $annotation_add_menu = $annotation_mb->Menu(-tearoff=>0, %mbcol);

$annotation_m->add('cascade', -label=>'Add', -underline=>0, -menu => $annotation_add_menu);

my $annotation_add_menu_manual_command = $annotation_add_menu->add('command',-label=>$annotatorName{manual}, -command=>sub { add_manual(\$main,\%annotation,$seqName,\$seq,\%sheet) });
my $annotation_add_menu_st_command = $annotation_add_menu->add('command',-label=>$annotatorName{ST}, -command=>sub { add_st(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_orf_command = $annotation_add_menu->add('command',-label=>$annotatorName{ORF}, -command=>sub { add_orf(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_polypy_command = $annotation_add_menu->add('command',-label=>$annotatorName{polypy}, -command=>sub { add_polypy(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_regexp_command = $annotation_add_menu->add('command',-label=>$annotatorName{regexp}, -command=>sub { add_regexp(\$main,\%annotation,$seqName,\$seq,\%sheet)});

#my $addsep = $annotation_add_menu->separator();
#  -activebackground=>$sheet{menu_button_active_bg_color});

my $annotation_add_menu_blast_command = $annotation_add_menu->add('command',-label=>$annotatorName{blast}, -command=>sub { add_blast(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_est_command = $annotation_add_menu->add('command',-label=>$annotatorName{EST}, -command=>sub { add_est(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_estorf_command = $annotation_add_menu->add('command',-label=>$annotatorName{ESTORF}, -command=>sub { add_estorf(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_glimmer_command = $annotation_add_menu->add('command',-label=>$annotatorName{glimmer2}, -command=>sub { add_glimmer(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_testcode_command = $annotation_add_menu->add('command',-label=>$annotatorName{testcode}, -command=>sub { add_testcode(\$main,\%annotation,$seqName,\$seq,\%sheet)});
my $annotation_add_menu_gff_command = $annotation_add_menu->add('command',-label=>$annotatorName{gff}." (stacked)", -command=>sub { add_gff(\$main,\%annotation,$seqName,\$seq,0,\%sheet)});
my $annotation_add_menu_gff_framed_command = $annotation_add_menu->add('command',-label=>$annotatorName{gff}." (framed)", -command=>sub { add_gff(\$main,\%annotation,$seqName,\$seq,1,\%sheet)});

#$annotation_m->add('command',-label=>'Import...',   -command=>sub { annotation_import(\$main,\%annotation) } );
$annotation_m->add('command',-label=>'Edit...', -underline=>0, -command=>sub { 

    my ($annotation_id,$annotation_type) = annotation_id_type_selected(\%annotation, \%sheet); # annotator type unused..
    
    if($annotation_id eq "error") {
	$WARNING && print "WARNING: No annotation selected before calling annotation-edit!\n";
    }
    $DEVEL && print "DEVEL: edit annotation $annotation_id\n";
    
    my $annotation_nr = annotation_what_nr(\%annotation,$annotation_id);
    if(!$annotation{note}->[$annotation_nr]) {
	$annotation{note}->[$annotation_nr]="";
    }
    
    annotation_edit(\$main,\$canvas,\%annotation,\$seq,\%sheet,$annotation_nr,$annotation{note}->[$annotation_nr]);
});

$annotation_m->add('command',-label=>'Merge...', -underline=>0, -command=>sub { annotation_merge(\$main,\$canvas,\%annotation,\$seqName,\$seq,\%sheet) } );
$annotation_m->add('command',-label=>'Remove...', -underline=>0, -command=>sub { annotation_remove(\$main,\$canvas,\%annotation,\%sheet) } );
$annotation_m->add('command',-label=>'Remove All...', -underline=>1, -command=>sub { annotation_removeAll(\$main,\$canvas,\%annotation,\$seq,\%sheet) } );

# Menu - Run 

$run_mb=$menu_bar->Menubutton(-text=>'Run',-underline=>0, %mbcol)->pack(-side=>'left');

my $run_m = $run_mb->Menu(-tearoff => 0, %mbcol);
$run_mb->configure(-menu => $run_m);

$run_m->add('command',-label=>'Glimmer 2.0',-underline=>0, -command=>sub { run(\$main,\$canvas,\$status_text,\%annotation,\$seqName,\$seq,\%sheet,'glimmer2') } );
$run_m->add('command',-label=>'Testcode',  -underline=>0, -command=>sub { run(\$main,\$canvas,\$status_text,\%annotation,\$seqName,\$seq,\%sheet,'testcode') } );
$run_m->add('command',-label=>'BLAST', -underline=>0, -command=>sub { run(\$main,\$canvas,\$status_text,\%annotation,\$seqName,\$seq,\%sheet,'blast') } );
$run_m->add('command',-label=>'ESTORF', -underline=>0, -command=>sub { run(\$main,\$canvas,\$status_text,\%annotation,\$seqName,\$seq,\%sheet,'estorf') } );
#$run_m->add('command',-label=>'~pEST',  -command=>sub { run(\$main,\$canvas,\%annotation,\$seqName,\$seq,\%sheet,'pest') } ); 

# Menu - Help

my $help_mb=$menu_bar->Menubutton(-text=>'Help',-underline=>0, %mbcol)->pack(-side=>'right');
my $help_m = $help_mb->Menu(-tearoff => 0, %mbcol);
$help_mb->configure(-menu => $help_m);

$help_m->add('command',-label=>'Help', -underline=>0, -command=>sub {&help_help(\$main,\%sheet)} );
$help_m->add('command',-label=>'About', -underline=>0, -command=>sub {&help_about(\$main,\%sheet)});

# Create main canvas
$canvas = $main->Canvas(-height=>$sheet{canvas_seq_height},-width=>600,-background=>$sheet{default_win_background})->pack(-expand => 'yes', -fill => 'both', -side=>'top');

# bind commands to canvas
$canvas->bind('annotation','<1>'=>sub {click_annotation(\$main,\$canvas,\%sheet,\%annotation,\$seq,\$status_text)});
$canvas->bind('annotation','<2>'=>sub{middle_click_annotation(\$main,\$canvas,\%annotation,$seqName,\$seq,\$status_text)});
$canvas->bind('annotation','<3>'=>sub {right_click_annotation(\$main,\$canvas,\%sheet,\%annotation,\$seqName,\$seq,\$status_text)});

$canvas->bind('guide','<3>'=>sub {right_click_guide(\$view_m)});

#$canvas->bind('sequence','<1>'=>sub {click_sequence($main,\$canvas,\$status_text)});

# $DEBUG && print "DEBUG: canvas is $canvas\n";
$canvas->bind('sequence','<1>'=> [\&click_sequence,\$main,\$canvas,\%sheet,Tk::Ev('x'),Tk::Ev('y'),\$status_text] );
$canvas->bind('sequence','<2>'=> [\&middle_click_sequence,\$main,\$canvas,\$seq,\%sheet,Tk::Ev('x'),Tk::Ev('y'),\$status_text] );
# add sequence right click, with cut, paste, add manual here, add as manual+blast, etc..

# $canvas->bind('sequence','<ButtonPress-1>'=> sub {$canvas->selectFrom('sequence' , "@".Ev('x').",".Ev('y'))});
#  $canvas->bind('sequence','<ButtonRelease-1>'=> \&click_sequence($main,\$canvas,Ev('x'), Ev('y'), \$status_text));

#$canvas->bind('all','<2>'=> sub { print "Canvas info: xview ",$canvas->xview," scrollregion ",$canvas->cget('scrollregion'),"\n";});

#
# Redraw by parsing the worksheet file, or by other method? There ought to be some sort of 1-1 mapping between graphical apperance and
# a saved work-sheet anyway, but it might slow operation. Use for something like "menu_view_redraw"?
#

# Return height/width after resize: $canvas->Width;

$scrollbar = $main->Scrollbar(-orient=>'horiz',-command=> [xview => $canvas])->pack(-side=>'top',-fill=>'x');
$canvas->configure(-xscrollcommand => [set => $scrollbar]);

# canvas->scan might also be fun..

# Status field frame
my $status_frame = $main->Frame(-relief => 'sunken', -background=>$sheet{default_win_background})->pack(-side=>'bottom', -fill=>'x');
my $status_field = $status_frame->Label(-textvariable => \$status_text, -justify=>'left', -foreground=>$sheet{status_text_col}, -background=>$sheet{default_win_background}, -padx=>5, -pady=>5,-wraplength=>790)->pack(-fill=>'x', -side=>'left');

# Note - it should really not be possible to do much before a worksheet is open and a sequence has been selected 

disable_menu_state; 

# read cl arguments of local scope, needing 
# Parse command line parameters
while (my $arg = shift @clargs) {
    if($arg eq "--new") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --new flag.\n";
	} else {
	    if(open_sequence_fasta_file(\$main,$arg,\$seqName,\$seq)) {
		$canvas->delete('all'); 
		
		$sheet{selected} = '';
		$sheet{file}=''; # done in file_close..
		
		draw_seq(\$canvas,$seqName,\$seq,\%sheet);

		my $unr = $annotation{unr};
		%annotation=();
		$annotation{unr}=$unr;
		
		$sheet{status}="clean"; # Just after opening a new seq, we consider the document clean.

		reset_menu_state;
		
		if($seqName) {
		    $main->title("A GUI - $seqName");
		    my $l=length($seq); $status_text="Opened sequence $seqName, $l bp.";
		}
	    }

	}
    } elsif ($arg eq "--open") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --open flag.\n";
	} else {
	    open_last_opened(\$main,\$canvas,\%annotation,\%sheet,\%plotstatus,\$seqName,\$seq, $arg);
	}
    } else {
	$WARNING && print STDERR "WARNING: Unrecognised command line parameter $arg. Continuing anyway.\n";
    }
}

$DEBUG && print STDERR "going to loop\n"; 
MainLoop();

############################### END MAIN ################################

# Canvas interaction

sub click_annotation {
  my $main=shift;
  my $canvas=shift;
  my $sheet=shift;
  my $annotation=shift;
  my $seq=shift;
  
  my $status_text=shift;
  
  my ($annotation_id, $annotation_type) = annotation_id_type_current($canvas);

  my $last_selection_id = $$sheet{selected};

  $$sheet{selected} = $annotation_id;

  $WARNING && print "Annotation $annotation_id made with $annotatorName{$annotation_type} was selected.\n";
  
  if ($annotation_type eq 'merge') {
      view_merged($main,$annotation,$annotation_id,$status_text);
  } else { # generalised view - blast, regexp, EST, ORF, polypy, ESTORF, genbankcds, manual, testcode, glimmer2
      view_annotation($main,$annotation,$annotation_id,$status_text);
  }

  # Redraw the selected annotation since selected annoatation changes the "focus"-stipple
  redraw_selected_annotation($canvas, $annotation, $sheet, $seq, $last_selection_id);
}

sub click_sequence {
  my $canvas = shift; # The bind call passes the canvas as the first argument, just for the heck of it. Nasty debug.. 
  # $DEBUG && print "The waste sent to click sequence is: $waste\n"; 
  my $main=shift;
  $canvas=shift;
  my $sheet = shift; #ref
  
  my $mouse_x_coord=shift;
  my $mouse_y_coord=shift;

  my $status_text=shift; # ref

  if($$sheet{zoom} ne "sequence") {
      print "You clicked on sequence $seqName.\n"; # access violation $seqName
      $$status_text = "Sequence $seqName (".length($seq)."bp)."; # access violation $seqName, $seq
      return;
  }

  my $xy_coord_string = "@".$$canvas->canvasx($mouse_x_coord).",".$$canvas->canvasy($mouse_y_coord);
  
  # $DEBUG && print "DEBUG: canvas is $canvas\n";
  # Check seq name
  my $seqName;
  my $seqpos;
  my $comp=0;

  my @taglist=$$canvas->gettags('current');
      foreach (@taglist) {	  
	  ($seqpos) = m/^xseqpos(\d+)/ if (m/^xseqpos/);       
	  next if ($_ eq 'current');
	  next if ($_ eq 'selected');
	  next if ($_ eq 'sequence');
	  $comp = 1 if ($_ eq 'comp');
	  $seqName=$_;
      }
  if($comp == 0) {
      $WARNING && print "You clicked on sequence $seqName at mouse coordinates $mouse_x_coord, $mouse_y_coord which is canvas coord $xy_coord_string which corresponds to sub-sequence position ",$$canvas->index('current',$xy_coord_string), ". With the attached xseqpos value $seqpos, the position is ",$seqpos+$$canvas->index('current',$xy_coord_string)," and thus the clicked (forward strand) base is ", substr($seq,$seqpos+$$canvas->index('current',$xy_coord_string)-1,1),"\n";
  } else {
      my $compbase = substr($seq,$seqpos+$$canvas->index('current',$xy_coord_string)-1,1);
      $compbase =~ tr/acgtACGT/tgcaTGCA/;
      $WARNING && print "You clicked on sequence $seqName at mouse coordinates $mouse_x_coord, $mouse_y_coord which is canvas coord $xy_coord_string which corresponds to sub-sequence position ",$$canvas->index('current',$xy_coord_string), ". With the attached xseqpos value $seqpos, the position is ",$seqpos+$$canvas->index('current',$xy_coord_string)," and thus the clicked (reverse strand) base is ",$compbase ,"\n";
  }

  $$sheet{seq_sel_start} = $seqpos+$$canvas->index('current',$xy_coord_string);
  
  $$status_text="Mark from sequence $seqName position ". $$sheet{seq_sel_start} .".";

  #mark text??
#  $$canvas->selectTo('sequence', "@".$dragto_x.);

}

sub middle_click_sequence {
  my $canvas = shift; # The bind call passes the canvas as the first argument, just for the heck of it. Nasty debug.. 
  # $DEBUG && print "The waste sent to click sequence is: $waste\n"; 
  my $main=shift;
  $canvas=shift;
  my $seq =shift; #ref
  my $sheet = shift; #ref
 
  my $mouse_x_coord=shift;
  my $mouse_y_coord=shift;
  
  my $status_text=shift; # ref

  my $xy_coord_string = "@".$$canvas->canvasx($mouse_x_coord).",".$$canvas->canvasy($mouse_y_coord);

  if($$sheet{zoom} ne "sequence") {
      print "You clicked on sequence $seqName.\n"; # access violation $seqName
      $$status_text = "Sequence $seqName (".length($seq). "bp)."; # access violation $seqName, $seq
      return;
  }
    
  # $DEBUG && print "DEBUG: canvas is $canvas\n";
  # Check seq name
  my $seqName;
  my $seqpos;
  my $comp = 0;
  
  my @taglist=$$canvas->gettags('current');
  foreach (@taglist) {	  
      ($seqpos) = m/^xseqpos(\d+)/ if (m/^xseqpos/);       
      next if ($_ eq 'current');
      next if ($_ eq 'selected');
      next if ($_ eq 'sequence');
      $comp = 1 if ($_ eq 'comp');
      $seqName=$_;
  }

  if($comp == 0) {
      $WARNING && print "You right-clicked on sequence $seqName at mouse coordinates $mouse_x_coord, $mouse_y_coord which is canvas coord $xy_coord_string which corresponds to sub-sequence position ",$$canvas->index('current',$xy_coord_string), ". With the attached xseqpos value $seqpos, the position is ",$seqpos+$$canvas->index('current',$xy_coord_string)," and thus the clicked (forward strand) base is ", substr($$seq,$seqpos+$$canvas->index('current',$xy_coord_string)-1,1),"\n";
  } else {
      my $compbase = substr($$seq,$seqpos+$$canvas->index('current',$xy_coord_string)-1,1);
      $compbase =~ tr/acgtACGT/tgcaTGCA/;
      $WARNING && print "You right-clicked on sequence $seqName at mouse coordinates $mouse_x_coord, $mouse_y_coord which is canvas coord $xy_coord_string which corresponds to sub-sequence position ",$$canvas->index('current',$xy_coord_string), ". With the attached xseqpos value $seqpos, the position is ",$seqpos+$$canvas->index('current',$xy_coord_string)," and thus the clicked (reverse strand) base is ",$compbase,"\n";
  }

  $$sheet{seq_sel_stop}=$seqpos+$$canvas->index('current',$xy_coord_string);

  # Could ofcourse do things based on comp/non-comp string here..
  if( ($$sheet{seq_sel_start} < $$sheet{seq_sel_stop})) {
      if($comp != 0) {
	  # complementary strand clicked in forward order...
	  $_ = reverse(split(/ */,substr($$seq, $$sheet{seq_sel_start}-1, $$sheet{seq_sel_stop}-$$sheet{seq_sel_start}+1)));
	  tr/atgcATGC/tacgTACG/;
	  print "The (revcomp) sequence on pos $$sheet{seq_sel_start}-$$sheet{seq_sel_stop} is $_.\n";
      } else {
	  $_ = substr($$seq, $$sheet{seq_sel_start}-1, $$sheet{seq_sel_stop}-$$sheet{seq_sel_start}+1);
	  print "The sequence on pos $$sheet{seq_sel_start}-$$sheet{seq_sel_stop} is: $_.\n";
      }
  } else {
      # reverse order

      $WARNING && print "WARNING: selection stop upstream of selection stop (in the upper strand sense). Ignoring and reversing..\n";

      my $tmp = $$sheet{seq_sel_start};
      $$sheet{seq_sel_start} = $$sheet{seq_sel_stop};
      $$sheet{seq_sel_stop} = $tmp;

      if($comp != 0) {
	  $_ = reverse(split(/ */,substr($$seq, $$sheet{seq_sel_start}-1, $$sheet{seq_sel_stop}-$$sheet{seq_sel_start}+1)));
	  tr/atgcATGC/tacgTACG/;
	  print "The (revcomp) sequence on pos $$sheet{seq_sel_stop}-$$sheet{seq_sel_start} is $_.\n";
      } else {
	  $_ = substr($$seq, $$sheet{seq_sel_start}-1, $$sheet{seq_sel_stop}-$$sheet{seq_sel_start}+1);
	  print "The sequence on pos $$sheet{seq_sel_start}-$$sheet{seq_sel_stop} is: $_.\n";      
      }
  }

  draw_seq_text($canvas,$seqName,$seq,$sheet);
  $$canvas->SelectionClear(-selection=>CLIPBOARD);
  $$canvas->clipboardClear;
  $$canvas->clipboardAppend($_);
  # $$canvas->SelectionNotify();

  $$status_text="Mark to sequence $seqName position ". $$sheet{seq_sel_stop} .".";
}

sub annotation_id_type_selected {
    my $annotation = shift;
    my $sheet = shift;

    if((!defined($$sheet{selected})) or $$sheet{selected} eq '') {
	$DEBUG && print "DEBUG WARNING: annotation_id_type_selected called while no annotation selected ($$sheet{selected}).\n";
	return "error";
    }

    my $selected_nr = $annotation_nr_cache{$$sheet{selected}};
    return ($$sheet{selected}, $$annotation{type}->[$selected_nr]);
}

sub annotation_id_type_current {
    my $canvas = shift;
    
    # Get tags for current annotation 
    my @taglist=$$canvas->gettags('current'); 
    my $annotation_type;
    my $annotation_id;

    TAG: foreach (@taglist) {
	next if ($_ eq 'current');
	next if ($_ eq 'selected');
	next if ($_ eq 'annotation');
	next if ($_ eq 'comp');
	next if m/^xseqpos/;
	foreach my $knownType (keys %annotatorName) {
	    if($_ eq $knownType) {
		$annotation_type=$_;
		next TAG;
	    }
	}
	$annotation_id=$_;		# Gotta be an id then.. 
    }
    return ($annotation_id,$annotation_type);
}

sub middle_click_annotation {
  my $main=shift;
  my $canvas=shift;
  my $annotation=shift;
  my $sequenceName=shift; # ref
  my $sequence=shift; # ref
  my $status_text=shift;
  
  my ($annotation_id,$annotation_type) = annotation_id_type_current($canvas);

  print "Annotation $annotation_id made with $annotatorName{$annotation_type} was selected.\n";

  # Call the ordinary annotation view.. 
  if ($annotation_type eq 'merge') {
      view_merged($main,$annotation,$annotation_id,$status_text);
  } else {
      view_annotation($main,$annotation,$annotation_id,$status_text);
  }
  
  my $annotation_nr = annotation_what_nr($annotation,$annotation_id);

  if($$annotation{note}->[$annotation_nr]) {
    print "Notes: $$annotation{note}->[$annotation_nr]\n";
  }
  
  my $annotated_sequence = substr($$sequence,$$annotation{start}->[$annotation_nr]-1,$$annotation{stop}->[$annotation_nr]-$$annotation{start}->[$annotation_nr]+1);
  
  my $annotated_protein;
  if( (not defined $$annotation{frame}->[$annotation_nr]) || $$annotation{frame}->[$annotation_nr]<3 ) {
      # Forward reading frame
      print "Annotated sequence $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] is $annotated_sequence.\n";

      $annotated_protein=nt2aa($annotated_sequence); # We could remove M and + if we'd like..
      print "Forward strand translation: $annotated_protein\n";

  } else {
      # Reverse reading frame
      my $rev_comp=reverse(split(/ */,$annotated_sequence));
      $rev_comp=~tr/atgcATGC/tacgTACG/;
      my $annotated_protein=nt2aa($rev_comp); # We could remove M and + if we'd like..
      print "Annotated sequence $$annotation{start}->[$annotation_nr] - $$annotation{stop}->[$annotation_nr] is (rev comp) $rev_comp.\n";
      
      print "Reverse strand translation: $annotated_protein\n";
  } 
}

sub right_click_guide {
    my $menu = shift; #
    #  clone main popup instead?

    $$menu->Popup(-popover => 'cursor', -popanchor => 'nw');
    
    # Flush the menu out - don't wanna wait for the annotations to redraw to select from the menu!
    $$menu->update();
    
#    my $main_popup_menu = $$main->Menu(-tearoff => 0, -activebackground=>$sheet{menu_button_active_bg_color});
    # commands to invoke main::mb_zoomIn and mb_zoomOut if not disabled (get config state)
  
}
  
sub right_click_annotation { 
  my $main=shift;
  my $canvas=shift;
  my $sheet=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $status_text=shift;

  # Select right-clicked sequence
  
  my ($annotation_id, $annotation_type) = annotation_id_type_current($canvas);

  my $last_selection_id = $$sheet{selected};

  $$sheet{selected} = $annotation_id;
  
  my $annotation_nr = annotation_what_nr($annotation,$annotation_id);

  # Popup right click menu
  my $annotation_popup_menu = $$main->Menu(-tearoff => 0, %mbcol);
  if(defined($annotation{name}->[$annotation_nr]) && $annotation{name}->[$annotation_nr] ne "")  {
      $annotation_popup_menu->add('command',-state=>'disabled', -label=>$annotation{name}->[$annotation_nr],  -foreground=>$sheet{menu_button_foreground_color});
      $annotation_popup_menu->separator();
  }
  
  $annotation_popup_menu->add('command',-label=>'Edit...', -command=>sub { annotation_edit($main,$canvas,$annotation,$seq,$sheet,$annotation_nr,$annotation{note}->[$annotation_nr]);});
  $annotation_popup_menu->add('command',-label=>'Remove...', -command=>sub { annotation_remove($main,$canvas,$annotation,$sheet);});
$annotation_popup_menu->add('command',-label=>'Add as manual', -command=>sub { add_selected_as_manual($main, $canvas, $annotation, $seqName, $seq, $sheet);});
  
  $annotation_popup_menu->separator();

  $annotation_popup_menu->add('command',-label=>'Blastn', -command=>sub {
      my %runstate = ();
      $runstate{blast_annotation}=1;
      $runstate{blast_annotationid}= $annotation{id}->[annotation_what_selected_nr($annotation,$sheet)];
      $runstate{blast_db} = 1;
      $runstate{blastn} = 1;

      $$sheet{blastlocaldefault} eq "wu" && ($runstate{blast_wu} = 1);
      
      $DEBUG && print "DEBUG: calling blastn for annotation id $runstate{blast_annotationid}.\n";

      run_external($main,$canvas,$status_text,0,\%runstate,$annotation,$seqName,$seq,$sheet);
  });

  if(defined $$annotation{frame}->[$annotation_nr]) {
      # Blastp will only work for oriented (aa) annotations
      $annotation_popup_menu->add('command',-label=>'Blastp', -command=>sub  {
	  my %runstate = ();
	  $runstate{blast_annotation}=1;
	  $runstate{blast_annotationid}=$annotation{id}->[annotation_what_selected_nr($annotation,$sheet)];
	  $runstate{blast_db} = 1;
	  $runstate{blastp} = 1;

	  $$sheet{blastlocaldefault} eq "wu" && ($runstate{blast_wu} = 1);

	  $DEBUG && print "DEBUG: calling blastp for annotation id $runstate{blast_annotationid}.\n";
	  run_external($main,$canvas,$status_text,0,\%runstate,$annotation,$seqName,$seq,$sheet);
      });
  }
  $annotation_popup_menu->separator();

  $annotation_popup_menu->add('command',-label=>'Qblastn', -command=>sub {qblast($main,$canvas,$sheet,$annotation,$seqName,$seq,$status_text,'blastn')});
  if(defined $$annotation{frame}->[$annotation_nr]) {
      # Blastp will only work for oriented (aa) annotations
      $annotation_popup_menu->add('command',-label=>'Qblastp', -command=>sub {qblast($main,$canvas,$sheet,$annotation,$seqName,$seq,$status_text,'blastp')});
      $annotation_popup_menu->add('command',-label=>'Qblast RPS CD', -command=>sub { qblastcd($main,$canvas,$sheet,$annotation,$seqName,$seq,$status_text)});
  }
  $annotation_popup_menu->add('command',-label=>'Qblastx', -command=>sub {qblast($main,$canvas,$sheet,$annotation,$seqName,$seq,$status_text,'blastx')});

  $annotation_popup_menu->Popup(-popover => 'cursor', -popanchor => 'nw');

  # Flush the menu out - don't wanna wait for the annotations to redraw to select from the menu!
  $annotation_popup_menu->update();

  # Redraw the selected annotation since selected annoatation changes the "focus"-stipple
  redraw_selected_annotation($canvas,$annotation,$sheet,$seq, $last_selection_id);
}

# File menu helpers
sub post_lastfiles_menu {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $plotstatus=shift;
    my $seqName=shift;
    my $seq=shift;
    my $lastfilesmenu = shift; 

    # $DEBUG && print "DEBUG: posting lastfiles menu for seqName $$seqName (ref $seqName), seq $$seq (ref $seq)\n";

    my @lastfiles = last_opened_files($sheet);

    if(@lastfiles > 0) {
	($$lastfilesmenu->menu)->delete(0,$$sheet{nolastfiles}); # empty menu
	my $min_entries = ( (@lastfiles < $$sheet{nolastfiles}) ? scalar(@lastfiles) : $$sheet{nolastfiles} );

	my $i;

	for ($i = 0; $i < $min_entries; $i++) {
	    # ignore any extra entries in lastfiles file
	    ($$lastfilesmenu->menu)->add('command', -label => $lastfiles[$i], %mbcol, -command=>[\&open_last_opened,$main,$canvas,$annotation,$sheet,$plotstatus,$seqName,$seq,$lastfiles[$i]]);
	}
    }
    return;
}

# File menu

sub file_new {
  my $main=shift;
  my $canvas=shift;
  my $sheet=shift;	# catch reference to sheet status
  my $seqName=shift; # ref
  my $seq=shift; # ref
  my $annotation=shift; # ref
  my $plotstate=shift; # ref

  # DEBUG
  #    $$sheet_status='dirty';
 
  $ok = file_close($main,$canvas,$sheet,$annotation,$seqName,$seq,$plotstate);
  if($ok) {
      $WARNING && print "File saved ok, or discarded from file_close.\n";
  } else {
      $WARNING && print "Cancel new.\n";
      return 0;			# Ok, so cancel new.
  }

  my $new_dialog=$$main->Toplevel;
  $new_dialog->title("A GUI - New worksheet");
#  $new_dialog->geometry('+10+0');
  $new_dialog->configure(-background=>'linen');
  my $new_dialog_frame=$new_dialog->Frame->pack(-expand=>'yes',-fill=>'both');
  my $new_dialog_label=$new_dialog_frame->Label(-text => "A new worksheet has been created.\nPlease select a sequence to work with.\n", -justify => 'left',-background=>$$sheet{default_win_background})->pack(-expand=>'yes',-fill=>"both");
  
  my $sequence_opened = open_sequence($main,$seqName,$seq,$sheet);
  if (!$sequence_opened) {
    print "No sequence selected - cancel new...\n";
    $new_dialog->destroy();
    return 0;
  }
  $DEVEL && print "DEVEL: Selected sequence ($seqName) now $$seqName.\nMore \"new\" things happen...\n";
  $DEVEL && print "DEVEL: And corresponding seq is $$seq.\n";
  
  # delete old seq first?
  $$canvas->delete('all'); 
  
  $$sheet{selected} = '';
  $$sheet{file}=''; # done in file_close..

#  my $seqLength = length($$seq);
#  draw_seq($canvas,$$seqName,$seqLength,$sheet);
  draw_seq($canvas,$$seqName,$seq,$sheet);

  my $unr = $$annotation{unr};
  %$annotation=();
  $$annotation{unr}=$unr;
  
  $$sheet{status}="clean"; # Just after opening a new seq, we could consider the document clean. This is a little dubious, but hey, there isn't really much of a problem here, right?
  
  $new_dialog->destroy;
  return 1;
}

sub open_last_opened {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $plotstatus=shift;
    my $seqName=shift;
    my $seq=shift;
    my $this_last_file = shift;

    $DEBUG && print "DEBUG: This last file is known as $this_last_file.\n";
    
    $DEBUG && ((ref $seq) || die "Seq is not a ref in open_last_opened!\n");

    $ok = file_close($main,$canvas,$sheet,$annotation,$seqName,$seq,$plotstatus);
    if($ok) {
	$WARNING && print "File saved ok, or discarded from file_close.\n";
    } else {
	$WARNING && print "Cancel new.\n";
	return 0;			# Ok, so cancel new.
    }

    my $status = file_load($main,$canvas,$annotation,$sheet,$plotstatus,$seqName,$seq,$this_last_file);
    if($status) {
	draw_seq($canvas,$$seqName,$seq,$sheet);                       # package following 8 lines in subroutine
	$$canvas->configure(-height=>$$sheet{canvas_seq_height});
	redraw_annotations($canvas,$annotation,$sheet,$seq);

	reset_menu_state;

	if($$seqName) {
	    $$main->title("A GUI - $$seqName [$$sheet{file}]");
	    my $l=length($$seq); $status_text="Loaded sequence $$seqName, $l bp.";
	}
    } else {
	# failed opening file as a worksheet save file - maybe its a genbank file?
	my $status = open_genbank($main, $annotation, $seqName, $seq, $sheet, $this_last_file);
	if($status) {
	    draw_seq($canvas,$$seqName,$seq,$sheet);
	    $$canvas->configure(-height=>$$sheet{canvas_seq_height});
	    redraw_annotations($canvas,$annotation,$sheet,$seq);

	    reset_menu_state;

	    if($seqName) {
		$$main->title("A GUI - $$seqName [$$sheet{file}]");
		my $l=length($$seq); $status_text="Loaded sequence $$seqName, $l bp.";
	    }
	}
    }
}


# Help menu

sub help_about {
  my $main=shift;
  my $sheet = shift;
  
  my $about_win = $$main->Toplevel;
  $about_win->title("A GUI - About");
#  $about_win->geometry('+0+0');
  $about_win->configure(-background=>'linen');
  
  my $about_text_frame=$about_win->Frame(-background=>$$sheet{default_win_background})->pack;
  my $about_text=$about_text_frame->Label(-text => "Gene finding GUI for semi-automagical annotation.\n\nBy Daniel Nilsson\ndaniel.nilsson\@cgb.ki.se\n2000-01-19\n\nAll reservable rights reserved.\n\nUse on your risk, if you have the right to do so - the author cannot be held liable\n for any damage, mental, reputational or otherwise, caused.\nAnd, hey, please try to have fun doing so!\n", -justify => 'center', -background=>$$sheet{default_win_background})->pack(-side => 'top');
  my $about_ok_button=$about_text_frame->Button(-text => 'Ok', -command => sub { $about_win->destroy })->pack(-side => 'bottom');
  
}

sub file_export_sequin {
  my $main=shift;
  my $seqName=shift; # ref
  my $seq=shift; # ref
  my $annotation=shift; # ref
  my $sheet = shift; # ref
  
  # Dialog for chosing what to export, and into what format..
  # 

  my %annotator_selection;
  foreach $annotator (keys %annotatorName) {
      $annotator_selection{$annotator} = 0;      
  }

  $annotator_selection{'merge'} = 1; # default
  select_levels($main, $canvas, \%annotator_selection, $sheet);
  
  export_sequin($main,$seqName,$seq,$annotation, \%annotator_selection);
}

sub file_export_gff {
  my $main=shift;
  my $seqName=shift; # ref
  my $seq=shift; # ref
  my $annotation=shift; # ref
  my $sheet = shift; # ref

  # Dialog for chosing what to export, and into what format..
  # 

  my %annotator_selection;
  foreach $annotator (keys %annotatorName) {
      $annotator_selection{$annotator} = 0;      
  }

  $annotator_selection{'merge'} = 1; # default
  select_levels($main,$canvas, \%annotator_selection, $sheet);
  
  export_gff($main,$seqName,$seq,$annotation, \%annotator_selection);
}

sub file_export_genbank {
  my $main=shift;
  my $seqName=shift; # ref
  my $seq=shift; # ref
  my $annotation=shift; # ref
  my $sheet = shift; # ref

  # Dialog for chosing what to export, and into what format..
  # 

  my %annotator_selection;
  foreach $annotator (keys %annotatorName) {
      $annotator_selection{$annotator} = 0;      
  }

  $annotator_selection{'merge'} = 1; # default
  select_levels($main,$canvas, \%annotator_selection, $sheet);
  
  my $selectedFile = export_genbank($main,$seqName,$seq,$annotation, \%annotator_selection);

  return $selectedFile;
}

sub help_help {
  my $main=shift;
  my $sheet = shift;

  my $help_win = $$main->Toplevel;
  $help_win->title("A GUI - Help");
#  $help_win->geometry('+0+0');
  $help_win->configure(-background=>'linen');
  
  my $help_text_frame=$help_win->Frame(-background=>$$sheet{default_win_background})->pack;
  my $help_text=$help_text_frame->Label(-text => "Gene finding GUI for semi-automagical annotation.\n\nBy Daniel Nilsson\ndaniel.nilsson\@cgb.ki.se\n2000-01-19\n", -justify => 'center', -background=>$$sheet{default_win_background})->pack(-side => 'top');
  my $help_ok_button=$help_text_frame->Button(-text => 'Ok', -command => sub { $help_win->destroy })->pack(-side => 'bottom');
}

# Annotation helpers

sub import_blast {
    my $main = shift;
    my $blastHtmlResults = shift;
    my $annotation = shift;
    my $seqName=shift;
    my $seq = shift;
    my $sheet = shift;
    my $annotation_id = shift;

    my $program = "";

    foreach $_ (split(/\n{1}/,$blastHtmlResults)) {
	if(/BLASTP/) {
	    $program = "blastp";
	} elsif(/BLASTN/) {
	    $program = "blastn";
	} elsif(/BLASTX/) {
	    $program = "blastx";
	}
    }

    defined ($annotation_id) || ($annotation_id = "");

    # open file, and place contents in ncbi_results, if html version..
    if($program eq "blastn") {
	parse_blastn_reply($main, \$canvas,$sheet,$seq,$annotation,$annotation_id,1e-5,$blastHtmlResults); # using global canvas
    } elsif ($program eq "blastx") {
	parse_blastx_reply($main, \$canvas,$sheet,$seq,$annotation,$annotation_id,0.001,$blastHtmlResults); # using global canvas
    } elsif ($program eq "blastp") {
	parse_blastp_reply($main, \$canvas,$sheet,$seq,$annotation,$annotation_id,0.001,$blastHtmlResults); # using global canvas
    } else {
	print "No blast program type designation found in file. If you are certain this is a valid ncbi blast HTML reply, this is a program error.\n";
    }

    level_layout($annotation, 'blast');
}

sub import_gff {
    my $main=shift;
    my $gffFile=shift;
    my $bioframed = shift;
    my $annotation=shift;
    my $seqName=shift; # NOT ref
    my $seq=shift;
    my $sheet=shift;
    
    $DEVEL && $DEBUG && ref $seqName && print "DEBUG: seqName is reference in import_gff.\n";

    use Bio::Tools::GFF;

    my $gffin = Bio::Tools::GFF->new(-file => $gffFile, -gff_version => 2);
    my $feature;

    my $nr=($$annotation{nr})?($$annotation{nr}):0;

    while($feature = $gffin->next_feature()) {
	my $frame;
	if($feature->strand == -1) {
	    $frame = 5-((length($$seq)-$feature->end)%3);
	} else {
	    $frame = ($feature->start%3)?$feature->start%3-1:2;
	}
	my $comment = $feature->primary_tag." ".$feature->source_tag;

	my $note;
	foreach $tag ($feature->all_tags()) {
	    $note .= $tag." = ".(join(' ',$feature->each_tag_value($tag)))."\n";
	 }
	$note .= join(' ', $feature->sub_SeqFeature());

	my $level = $frame - 5;

	# open new annotation type (add req props to Annotator)
	# I dont really have the time to make this a full-fledged oo "annotator->new()" but given time to evolve..
	my $type = "gff_".$feature->primary_tag;
	my $color = 'aquamarine3'; 

	if ($type eq 'gff_megablasthit' ){
	    $color = 'blue';
	    # config option..
#	    $color = $$main->chooseColor(-title => "Choose an annotation color", -initialcolor => $color);
	    my ($name) = ($note =~ /ame\s+=\s+([\w\d._+-\[\]\(\)]+)\s+/);

	    add_gff_annotator($sheet, $type, 3, $bioframed); # height 3 as a default for a gff-derived type (non-bioframed)
	    annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type => $type, start => $feature->start, stop => $feature->end, frame => $frame, comment => $comment, note => $note, color => $color, level => $level, name => $name);
	    $nr++;
	} else {

	    add_gff_annotator($sheet, $type, 3, $bioframed); # height 3 as a default for a gff-derived type (non-bioframed)
	    $color = $annotatorDefaultColor{$type};
	    
	    # $color = $$main->chooseColor(-title => "Choose an annotation color", -initialcolor => $color);
	    # have a list ready to assing to each new gff-type?
	    
	    annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type => $type, start => $feature->start, stop => $feature->end, frame => $frame, comment => $comment, note => $note, color => $color, level => $level);
# 	$level += $annotatorDisplayLevel{$type}	    
	    $nr++;
	}
    }
    
    $gffin->close();

    level_layout_gff($annotation);

    redraw_annotations(\$canvas,$annotation,$sheet,$seq);  # using global canvas
}

sub import_glimmer {
    my $main=shift;
    my $glimFile=shift;
    my $annotation=shift;		# annotation reference..
    my $seqName=shift; # NOT ref
    my $seq=shift;
    my $sheet=shift;
    
    $DEVEL && $DEBUG && ref $seqName && print "DEBUG: seqName is reference in import_glimmer.\n";
    $DEVEL && print "DEVEL import_glimmer: $seqName appears to be current sequence...\n";
    
    my $nr=($$annotation{nr})?($$annotation{nr}):0;

    my $inPutatives=0;

    my $subjectName = "";
    
    open GLIM, "<$glimFile" or die "DEBUG: Could not open $glimFile. Exiting...\n";
    
  GLIM:  while(<GLIM>) {
      chop;

      my $glimStart=0;
      my $glimStop=0;
      my $glimFrame=0;
      my $glimComment="";
      my $glimFound=0;
      
      # Use color to signify the quality...
      my $color;
      
      if(m/^Putatives/) {
	  $inPutatives=1;
	  # Ok, we have found a putative section. Now, first get the subjectName, then move on to next row 
	  # for some actual action.
	  ($subjectName)=m/^Putatives in\s+(.*)\s+\[start end frame comments\]/;
	  $DEBUG && print print "DEBUG: new putative list on subjectName $subjectName in file $glimFile.\n"; 
	  # Deactivated. Take care not to exec on $subjectName without quotes..
	  # $subjectName=~s/ //g; 
	  next GLIM;
      } elsif (m/DEBUG/) {
	  $inPutatives=0;
      }
      
      if($inPutatives) {
	  # Obtain glim.pl output data..
	  ($glimStart,$glimStop,$glimFrame,$glimComment)=m/(\d+)\s+(\d+)\s+(\d+)\s+(.+)/;
	  
	  # Deal with the discrepancy in start/stop position handling.
	  # First, swap start and stop on any >2 orf.
	  if($glimFrame>2) {
	      my $tmp=$glimStart;
	      # Adjust the "stop" -3, since glimmer reports "the position of the last base *before* the stop codon". (Make it 3nt longer i e.)
	      $glimStart=$glimStop-3;
	      $glimStop=$tmp;

	      #if($glimFrame == 3) {
	      #    $glimFrame = 4;
	      #} elsif($glimFrame == 4) {
	      #  $glimFrame = 5;	    
	      #} elsif($glimFrame == 5) {
	      #    $glimFrame = 3;	    
	      #}
	      
	  } else {
	      # Adjust the "stop" +3, since glimmer reports "the position of the last base *before* the stop codon".
	      $glimStop+=3;
	  }

	  # Glimmer frames do not correspond to my definition; conversion follows:
	  if($glimFrame < 3) {
	      # update frame -- forward pass
	      $glimFrame=($glimStart%3)?($glimStart%3)-1:2;
	  } else {
	      $glimFrame = 5-((length($$seq)-$glimStop)%3);
	  }
	  if ($subjectName ne $seqName) {
	      $DEBUG && print "DEBUG: Found glimmer prediction on $subjectName, when expecting prediction to current seq $seqName.\n";
	      next GLIM;
	  } else {
	      # print "Found glimmer2 hit $glimStart-$glimStop ($glimFrame) for $subjectName.\n";
	  }
	  
#      if ($glimComment=~m/^\s*\[\#\d+\]\s*$/ or $glimComment=~m/^\s*\[LowScoreBy\s+\#\d+\s+L\=\d+\s+S\=\d+\]\[\#\d+\]\s*$/ or  ) {		
	  if ($glimComment=~m/\[OlapWith\s+\#\d+\s+L\=\d+\s+S\=\d+\]/ or $glimComment=~m/\[ShadowedBy\s+\#\d+\]/ or $glimComment=~m/\[DelayedBy\s+\#\d+\s+L\=\d+\]/ or $glimComment=~m/\[ShorterThan.+\]/) {
	      # Shadowed (etc.) putative...
	      $color="red3";
	  } else {
	      # Un-shadowed (etc.) putative...
	      $color="red1";
	  }

	  my $level;
	  # Check if the startcodon is standard or alternative and assign level	
	  if($glimFrame<3) {
	      $level=$glimFrame;	
	      if (substr($$seq,$glimStart-1,3) eq "ATG") {
		  $glimComment.=' start: ATG';
	      } elsif (substr($$seq,$glimStart-1,3)=~m/[CGT]TG/) {     
		  $glimComment.= ' alt start: ';
		  $glimComment.=substr($$seq,$glimStart-1,3);
	      }
	  } else {
	      $level=$glimFrame - 5; 

	      my $startcodon=substr($$seq,$glimStop-3,3);
	      $startcodon=~tr/atgcATGC/tacgTACG/;
	      my $revcomp=join('',reverse(split(/ */,$startcodon))); 
	      if($revcomp ne 'ATG') {
		  $glimComment.= " alt";
	      }
	      $glimComment.=" start: ";
	      $glimComment.= $revcomp;  
	  }

	  # add annotation
	  annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type=> 'glimmer2', start => $glimStart, stop => $glimStop, frame => $glimFrame, comment => $glimComment, color => $color, level => $level);
	  $nr = $$annotation{nr};
	  $$sheet{status} = "dirty"; # Sheet is now dirty.
      }
  }
    close(GLIM);

    redraw_annotations(\$canvas,$annotation,$sheet,$seq);  # using global canvas
}

sub do_testcode {
    my $main=shift;
    my $annotation=shift;
    my $seqName=shift;
    my $seq=shift;
    my $sheet=shift;

    my $shortestOrf=200; # testcode statistic not confident below 200nt threshold.
    my $testcode_orf_average_threshold = .74;
    my $testcode_orf_confident_threshold = .95;
#    my $increment = 3;
#    my $mwinl = 200;

    my $nr=($$annotation{nr})?($$annotation{nr}):0;

    # locate orfs

    # Translate trinucleuotides -> aa in all six reading frames.
    my @phase;
    my @phaseOffset=(0,1,2,2,1,0);
    
    $phase[0]=nt2aa($$seq);
    $phase[1]=nt2aa(substr($$seq,1));
    $phase[2]=nt2aa(substr($$seq,2));

    # Same story with the complementary strand.
    my $complementSeq=$$seq;
    # The complement has complementary bases, and in the
    # translation perspective, reverse orientation.
    $complementSeq =~ tr/atgcATGC/tacgTACG/;
    $complementSeq = join('',reverse(split(/ */,$complementSeq)));
    $phase[3]=nt2aa(substr($complementSeq,2));
    $phase[4]=nt2aa(substr($complementSeq,1));
    $phase[5]=nt2aa($complementSeq);

    for(my $i=0;$i<6;$i++) {
	my @orf=();
	while($phase[$i]=~m/((?:M\w*\+)|(?:^\w*\+)|(?:M\w*$))/g) {
	    push @orf,$1;
	    #  print "DEBUG: longest hits from phase $i: $1\n";
	}
	
	my $startAaPos;
	my $stopAaPos;
	my $startNtPos;
	my $stopNtPos;

	foreach my $orf (@orf) {
	    if(3*length($orf) > $shortestOrf) {
		$startAaPos=index($phase[$i],$orf);
		$startAaPos != -1 || ($WARNING && print "WARNING: Index on phase $i failed for $orf.\n");
		$stopAaPos=$startAaPos+length($orf);

		# average testcode over orfs, according to strand
		my $orf_testcode;
		if($i<=2) {
		    $startNtPos=3*$startAaPos+$phaseOffset[$i]+1; # to 1-based nt-pos
		    $stopNtPos=3*$stopAaPos+$phaseOffset[$i]; # to 1-based nt-pos
		    $DEBUG && $DEVEL && print "Calling testcode for forward strand orf $startNtPos to $stopNtPos (".($stopNtPos-$startNtPos+1)." bp):";
		    my $orf_seq = substr($$seq, $startNtPos - 1, $stopNtPos - $startNtPos);
		    $orf_testcode = testcode( $orf_seq );
		    $DEBUG && $DEVEL && print " score $orf_testcode.\n" ;
		    $DEBUG && $DEVEL && print "That would be the following sequence: ".$orf_seq.".\n";
		} else {
		    $startNtPos=3*$startAaPos+$phaseOffset[$i]+1; # to 1-based nt-pos
		    $stopNtPos=3*$stopAaPos+$phaseOffset[$i]; # to 1-based nt-pos
		    $DEBUG && $DEVEL && print "Calling testcode for rev strand orf $startNtPos to $stopNtPos (".($stopNtPos-$startNtPos+1)." bp):";
		    my $orf_seq = substr( $complementSeq, $startNtPos - 1, $stopNtPos-$startNtPos);
		    $orf_testcode = testcode( $orf_seq );
		    $DEBUG && $DEVEL && print " score $orf_testcode.\n" ;
		    $DEBUG && $DEVEL && print "That would be the following sequence: ".$orf_seq.".\n";

		    # switch to forward strand coords before adding annotation..
		    $startNtPos=-3*$stopAaPos+length($$seq)+1-$phaseOffset[$i]; # to 1-based nt-pos
		   $stopNtPos=-3*$startAaPos+length($$seq)-$phaseOffset[$i]; # to 1-based nt-pos, inverted start & stop always to have lowest number as NtStart.
		}

		# add only orfs over a certain threshold
		if($orf_testcode >= $testcode_orf_average_threshold) {

		    my $color = 'SkyBlue3';
		    my $note = "";
		    my $comment = sprintf "ORF Testcode score %.2f", $orf_testcode;
		    if ($orf_testcode >=  $testcode_orf_confident_threshold) {
			$comment .= " (predicted Coding)";
		    } else {
			$comment .= " (borderline score - No Opinion)";
		    }
		    my $frame=$i;

		    my $level;
		    if($i<3) {
			$level = $frame;
		    } else {
			$level = $frame - 5;
		    }

		    $nr = $$annotation{nr}; # update nr

		    if($DEVEL) {
			annotation_add_new($annotation, type=> 'testcode', id => "ab_$nr", seqName => $seqName, start => $startNtPos, stop => $stopNtPos, comment => $comment, frame => $frame, color => $color, level => $level, note => $note, name => "ab_$nr" );
		    } else {
			annotation_add_new($annotation, type=> 'testcode', id => "ab_$nr", seqName => $seqName, start => $startNtPos, stop => $stopNtPos, comment => $comment, frame => $frame, color => $color, level => $level, note => $note);
		    }

		    $$sheet{status} = 'dirty';
		}
	    }
	}
    }
    
    redraw_annotations(\$canvas,$annotation,$sheet,$seq); # using global canvas
}

sub import_testcode {
    my $main=shift;
    my $testcode_file_forward=shift;
    my $testcode_file_reverse=shift;
    my $annotation=shift;
    my $seqName=shift;
    my $seq=shift;
    my $sheet=shift;

    my $shortestOrf=150;
    my $testcode_orf_average_threshold = 7;

    my $increment = 3;
    my $mwinl = 200;

    my @testcode_forward;
    my @testcode_reverse;

    my $nr=($$annotation{nr})?($$annotation{nr}):0;

    my $status = open(TCFORWARD, "<$testcode_file_forward");
    if ($status == 0) {
	print "ERROR opening testcode file $testcode_file_forward!"; # WARNING
	close TCFORWARD;
	return;
    }

    $DEBUG && print "DEBUG: importing testcode results from $testcode_file_forward..\n";
    
    my $flen = 0;
    while(<TCFORWARD>) {
	chop;
	$testcode_forward[$flen] = $_;      
	$testcode_forward[$flen+1] = $_;
	$testcode_forward[$flen+2] = $_;
	$flen+=$increment;
    }

    close TCFORWARD;

    $status = open(TCREVERSE, "<$testcode_file_reverse");
    if ($status == 0) {
	print "ERROR opening testcode file $testcode_file_reverse!"; # WARNING
	close TCREVERSE;
	return;
    }

    $DEBUG && print "DEBUG: importing testcode results from $testcode_file_reverse..\n";

    my $rlen = 0;
    while(<TCREVERSE>) {
	chop;
	$testcode_reverse[$rlen] = $_;
	$testcode_reverse[$rlen+1] = $_;
	$testcode_reverse[$rlen+2] = $_;
	
	$rlen+=$increment;
    }

    close TCREVERSE;
   
    if( $rlen != $flen || $flen+$increment*$mwinl+3 < length($$seq) ) { # AND then they use integer-division...
	$WARNING && print "WARNING: error in testcode sequence length. With increment $increment and mwinl $mwinl, forward is $flen, reverse is $rlen and sequence ",length($seq),".\n";
	return;
    }
    
    for(my $i = 0; $i<(length($seq)-$flen); $i++) { 
	$testcode_forward[$flen+$i] = 0;
	$testcode_reverse[$flen+$i] = 0;
    }
    
    # locate orfs

    # Translate trinucleuotides -> aa in all six reading frames.
    my @phase;
    my @phaseOffset=(0,1,2,2,1,0);
    
    $phase[0]=nt2aa($$seq);
    $phase[1]=nt2aa(substr($$seq,1));
    $phase[2]=nt2aa(substr($$seq,2));

    # Same story with the complementary strand.
    my $complementSeq=$$seq;
    # The complement has complementary bases, and in the
    # translation perspective, reverse orientation.
    $complementSeq =~ tr/atgcATGC/tacgTACG/;
    $complementSeq = join('',reverse(split(/ */,$complementSeq)));
    $phase[3]=nt2aa(substr($complementSeq,2));
    $phase[4]=nt2aa(substr($complementSeq,1));
    $phase[5]=nt2aa($complementSeq);

    for(my $i=0;$i<6;$i++) {
	my @orf=();
	while($phase[$i]=~m/(M\w*\+)/g) {
	    push @orf,$1;
	    #  print "DEBUG: longest hits from phase $i: $1\n";
	}
	
	my $startAaPos;
	my $stopAaPos;
	my $startNtPos;
	my $stopNtPos;

	foreach (@orf) {
	    if(3*length($_) > $shortestOrf) {
		$startAaPos=index($phase[$i],$_);
		$startAaPos != -1 || ($WARNING && print "WARNING: Index on phase $i failed for $_.\n");
		$stopAaPos=$startAaPos+length($_);

		# average testcode over orfs, according to strand
		my $orf_testcode_sum;
		if($i<=2) {
		    $startNtPos=3*$startAaPos+$phaseOffset[$i]+1; # to 1-based nt-pos
		    $stopNtPos=3*$stopAaPos+$phaseOffset[$i]; # to 1-based nt-pos
		    $orf_testcode_sum = sum(@testcode_forward[$startNtPos..$stopNtPos]);
		} else {
		    $startNtPos=-3*$stopAaPos+length($$seq)+1-$phaseOffset[$i]; # to 1-based nt-pos
		    $stopNtPos=-3*$startAaPos+length($$seq)-$phaseOffset[$i]; # to 1-based nt-pos, inverted start & stop always to have lowest number as NtStart.
		    $orf_testcode_sum = sum(@testcode_reverse[$startNtPos..$stopNtPos]);
		}

		# add only orfs over a certain threshold
		if(($orf_testcode_sum/($stopNtPos-$startNtPos+1) ) > $testcode_orf_average_threshold) {

		    my $color = 'SkyBlue3';
		    my $note = "";
		    my $comment = sprintf "ORF Testcode score %.2f", ($orf_testcode_sum/($stopNtPos-$startNtPos+1));
		    my $frame=$i;

		    my $level;
		    if($i<3) {
			$level = $frame;
		    } else {
			$level = $frame - 5;
		    }

		    $nr = $$annotation{nr}; # update nr

		    $DEVEL && print "DEVEL: Adding annotation nr $nr; uid $$annotation{uid}->[$nr], start $$annotation{start}->[$nr], stop $$annotation{stop}->[$nr], frame $$annotation{frame}->[$nr], color $$annotation{color}->[$nr], comment $$annotation{comment}->[$nr], note $$annotation{note}->[$nr].\n";		    
		    
		    if($DEVEL) {
			annotation_add_new($annotation, type=> 'testcode', id => "ab_$nr", seqName => $seqName, start => $startNtPos, stop => $stopNtPos, comment => $comment, frame => $frame, color => $color, level => $level, note => $note, name => "ab_$nr" );
		    } else {
			annotation_add_new($annotation, type=> 'testcode', id => "ab_$nr", seqName => $seqName, start => $startNtPos, stop => $stopNtPos, comment => $comment, frame => $frame, color => $color, level => $level, note => $note);
		    }

		    $$sheet{status} = 'dirty';
		}
	    }
	}
    }
    
    redraw_annotations(\$canvas,$annotation,$sheet,$seq); # using global canvas
}


sub import_estorf {
  my $main=shift;
  my $pmestFile=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift; # ref
  my $sheet=shift;

  my $orfClass;

  open PMEST,"<$pmestFile" or die "ERROR: Could not open $pmestFile for reading: $!\n";
  $DEBUG && print "DEBUG: importing estorf annotation from $pmestFile..\n";

ORF: while(<PMEST>) {
    chop;
    
    my $subjectName;
    my $orfStart;
    my $orfStop;
    my $orfPhase;

    #  if(m/^Hit/) {
    #   # Might wanna catch the subject sequence locations for sanity checks - or should that be incorporated in the pmest program directly?
    # } 
    if(m/Real good hit/) {
      $orfClass=1;
    } elsif (m/Pretty similar hit/) {
      $orfClass=2;
    }
    
    if (m/touches/) {
      # print "$subjectName[$nEST][$hit[$nEST]] touches ORF $orfStartSave[$nEST][$hit[$nEST]][$orfCount] - $orfStopSave[$nEST][$hit[$nEST]][$orfCount] $orfPhase[$nEST][$hit[$nEST]][$orfCount]\n";
      ($subjectName,$orfStart, $orfStop, $orfPhase)=m/^([\w\d\.]+)\s+touches ORF\s+(\d+)\s+\-\s+(\d+)\s+(\d+)/;
      $subjectName=~s/ //g;
      
      # I think orfStart is 0-based, whereas glimStart is 1-based. Since this is for human-readability, we choose 1-based cooridnates.
      $orfStart++;
      $orfStop++; # Then I also would assume that in the name of consitency, stop should 0-based, but that seems not to be the case... POTENTIAL ERROR?
      # Resolved 000131 DN - error in pmest...
      if ($subjectName ne $seqName) {
	next ORF;
      } else {
	# print "Found ESTORF hit $orfStart-$orfStop ($orfPhase) for $subjectName.\n";
      }
      
      #print "DEBUG:$subjectName touches ORF $orfStart - $orfStop $orfPhase\n";

      my $nr;
  
      if($$annotation{nr}) {
	$nr=$$annotation{nr};
      } else {
	$nr=0;
	$$annotation{nr}=0;
      }

      # use add_new instead - but perhaps we should resolve the "redundant ESTORF" issue already?
      $$annotation{uid}->[$nr] = $$annotation{unr};
      $$annotation{unr}++;

      $$annotation{id}->[$nr]="ab_$nr";
      $annotation_nr_cache{$$annotation{id}->[$nr]} = $nr;
      $$annotation{seqName}->[$nr]=$seqName; # Pretty redundant.. Remove?
      $$annotation{type}->[$nr]='ESTORF';
      $$annotation{start}->[$nr]=$orfStart;
      $$annotation{stop}->[$nr]=$orfStop;
      $$annotation{frame}->[$nr]=$orfPhase;
      $$annotation{comment}->[$nr]="Hit class $orfClass";
      if($orfClass==1) {
	$$annotation{color}->[$nr]='blue1';
      } else {
	$$annotation{color}->[$nr]='blue3';
      }

      $$annotation{nr}++;      
      $$sheet{status}="dirty"; # Sheet is now dirty.

      if($orfPhase<3) {
	$$annotation{level}->[$nr]=$$annotation{frame}->[$nr];
      } else {
#	$$annotation{level}->[$nr]=-($annotatorLevel{ESTORF}+($$annotation{frame}->[$nr]-3))
	$$annotation{level}->[$nr] = $$annotation{frame}->[$nr] - 5; 
      }
    }
  }
  close(PMEST);

  my %removestate;
  $removestate{estorfredundant} = 1;
  remove_all($annotation,\%removestate,$seq,$sheet);

  redraw_annotations(\$canvas,$annotation,$sheet,$seq);  # using global canvas

}

sub import_est {
    my $main=shift;
    my $pmestFile=shift;
    my $annotation=shift;
    my $seqName=shift;
    my $seq = shift; 
    my $sheet=shift;

    open PMEST,"<$pmestFile" or die "ERROR: Could not open $pmestFile for reading: $!\n";
    $DEBUG && print "DEBUG: importing est annotation from $pmestFile..\n";

  EST: while(<PMEST>) {
      chop;
      if(m/^Hit/) {
	  my ($subjectName,$subjectHitBegin,$subjectHitEnd,$subjectHitStrand,$hitP,$hitIds,$hitLength,$hitIdFreq,$queryHitBegin,$queryHitEnd,$queryHitStrand)=m/^Hit:\s+(.+)\s+(\d+)\-(\d+)\((.+?)\)\s+P\=\s+([0-9\.e\-\+]+),ids\=\s+(\d+)\/(\d+)\(([.\d]+?)\%\)\s+Query:\s+(\d+)\-(\d+)\((.+?)\)/;
	  
	  if($subjectName ne $seqName) {
	      next EST;
	  }
	 
	  # So, this is an EST that hit the current sequence.
	  # Read next row for EST name etc.
	  $_=<PMEST>;
	  chop;
	  my ($estName,$hitClassText)=m/^(.+)\s+on.+\((.+)\)/;     
	  $DEVEL && print "DEVEL: Found $hitClassText with $estName.\n";
	  
	  my $nr;
	 
	  if($$annotation{nr}) {
	      $nr=$$annotation{nr};
	  } else {
	      $nr=0;
	      $$annotation{nr}=0;
	  }
      
	  my $color;
	  if($hitClassText=~m/good/) {
	      $color = "yellow1";
	  } elsif($hitClassText=~m/similar/) {
	      $color = "yellow3";
	  } else {
	      # Shaky hit - implement adjustable parameter for P-levels?
	      # Current policy - don't even add.
	      next EST;
	      # $$annotation{color}->[$nr]="yellow3";
	  }
	 
	  my $frame;
	  if($estName =~ /5{1}\'{1}/ ) {
	      $frame = 0;
	  } elsif ($estName =~ /3{1}\'{1}/ ) {
	      $frame = 3;
	  } else {
	      $frame = $AnnotationBoxFrame; # Box
	      $DEBUG && print "DEBUG: Could not determine orientation for $estName. Boxing it!\n";
	  }

	  if($DEVEL) {
	      annotation_add_new($annotation, type=> 'EST', id => "ab_$nr", seqName => $seqName, start => $subjectHitBegin, stop => $subjectHitEnd, comment => "$estName\n$queryHitBegin-$queryHitEnd($queryHitStrand) hit $subjectHitBegin-$subjectHitEnd($subjectHitStrand),P=$hitP,ids=$hitIds/$hitLength($hitIdFreq)", frame => $frame, color => $color, name => "ab_$nr");
	  } else {
	      annotation_add_new($annotation, type=> 'EST', id => "ab_$nr", seqName => $seqName, start => $subjectHitBegin, stop => $subjectHitEnd, comment => "$estName\n$queryHitBegin-$queryHitEnd($queryHitStrand) hit $subjectHitBegin-$subjectHitEnd($subjectHitStrand),P=$hitP,ids=$hitIds/$hitLength($hitIdFreq)", color => $color, frame => $frame);
	  }
	  $$sheet{status}="dirty"; # Sheet is now dirty.   
      }
  }
    
    level_layout($annotation, 'EST');
  
    redraw_annotations(\$canvas,$annotation,$sheet,$seq); # using global canvas
}

sub level_layout_all {
    my $annotation = shift;
    
    foreach $annotator (keys %annotatorOrder) {
	if($annotatorBioframed{$annotator} == 0) {
	    level_layout($annotation, $annotator);
	}
    }
} 

sub level_layout_gff {
    my $annotation = shift;
    
    foreach $annotator (keys %annotatorOrder) {
	if($annotator =~ /^gff/) {
	    level_layout($annotation, $annotator);
	}
    }
} 

sub level_layout {
    my $annotation = shift;
    my $type = shift;
    
    # A simplistic layout scheme, that will throw still overlapping sequences in the lowest frame..

#    my $baseLevel = $annotatorLevel{$type};
    my $baseLevel = 0; # Display levels introduced..

    my $stackingHeight = $annotatorHeight{$type} - 1; # nr of levels to stack on - 1

    my @frames;
#    if(($type eq 'EST') or ($type eq 'blast') or ($type eq 'regexp') or ($type eq 'gff')) {
    if($annotatorBioframed{$type} == 0) { # includes polypy -- but never called?!
	$DEVEL && print "DEVEL: type $type for level layout.\n";
	@frames = (0,3);
    } else  {
	@frames = ( 0 ); # no other annotation than EST really uses this feature yet, but in case that changes, better be prepared...
    }

    foreach my $currentFrame (@frames) { 
	$DEVEL && print "DEVEL: layouting on frame $currentFrame.\n";
	my @nrOfType = ();
	my @sortedOfType = ();

	for(my $i = 0; $i < $$annotation{nr}; $i ++) {
	    my $frame;
	    if ($type eq 'EST') { # still needed?!
		$frame = $$annotation{frame}->[$i];
#	    if($type eq 'blast' or $type eq 'regexp' or $type =~ m/^gff/) {
	    } elsif($annotatorBioframed{$type} == 0 && $annotatorDirected{$type} == 2) {
		if($$annotation{frame}->[$i] < 3) {
		    $frame = 0;
		} elsif($$annotation{frame}->[$i] < 6) {
		    $frame = 3;
		} elsif($$annotation{frame}->[$i] == 6) {
		    $frame = $AnnotationBoxFrame;
		}
	    } elsif($annotatorBioframed{$type} == 0 && $annotatorDirected{$type} == 1) {
		if($$annotation{frame}->[$i] < 3) {
		    $frame = 0;
		} elsif($$annotation{frame}->[$i] < 6) {
		    $frame = 3;
		} else {
		    $WARNING && print STDERR "WARNING: boxframe detected for directed non-bioframed annotation ($type).\n";
		}
	    }
	    if( $$annotation{type}->[$i] eq $type && ( $frame == $currentFrame || $frame == $AnnotationBoxFrame ) ) {
		push @nrOfType, $i;
	    }
	}
    
	if (@nrOfType == 0) {
	    $DEVEL && print "DEVEL: no $type annotations found!\n";
	    next;
	}

	@sortedOfType = sort { $$annotation{start}->[$a] <=> $$annotation{start}->[$b] } @nrOfType;

	my @lastOnLevel = split (/ */, "0"x($stackingHeight+1));
       
	$DEVEL && print "DEVEL: Placing $$annotation{id}->[ $sortedOfType[0] ] on level ",$baseLevel,"; it starts on $$annotation{start}->[ $sortedOfType[0] ].\n";
	$lastOnLevel[0] = $$annotation{stop}->[ $sortedOfType[0] ];
	
	if($$annotation{frame}->[$sortedOfType[0]] < 3) {
	    $$annotation{level}->[$sortedOfType[0]] = $baseLevel;
	} elsif ($$annotation{frame}->[$sortedOfType[0]] < 6)  {
	    $$annotation{level}->[$sortedOfType[0]] =  - $baseLevel;
	} elsif ($$annotation{frame}->[$sortedOfType[0]] == $AnnotationBoxFrame)  {
	    $$annotation{level}->[$sortedOfType[0]] = $baseLevel;
	} else {
	    $WARNING && print "WARNING: Annotation $$annotation{id}->[ $sortedOfType[0] ] has incorrect frame $$annotation{frame}->[$sortedOfType[0]], and was not assigned a level in level_layout for $type.\n";
	}

	my $currentIncrement = 0;
	
	for (my $i = 1; $i< @sortedOfType; $i++) {
	    if( $$annotation{start}->[ $sortedOfType[$i] ] < $$annotation{stop}->[ $sortedOfType[$i-1] ] ) {
		$currentIncrement++;
		if($currentIncrement != 1) {
		    # Can we go back to a lower level? As low a level as possible!
		    for(my $j = $stackingHeight; $j>0; $j--) {
			if ( $currentIncrement > $j && $$annotation{start}->[ $sortedOfType[$i] ] > $$annotation{stop}->[ $sortedOfType[$i-$j-1] ] ) {
			    $currentIncrement -= $j+1;
			}
		    }
		}
	    } else {
		$currentIncrement = 0;
	    }
	    
	    # Is this, lowest possible level, still being obscured by a last on level annotation
	    while(($currentIncrement <= $stackingHeight) && ($$annotation{start}->[ $sortedOfType[$i] ] < $lastOnLevel[$currentIncrement]) ) {
	        $DEVEL && print "DEVEL: A last on level obscuring was prevented for $$annotation{id}->[ $sortedOfType[$i] ] on level ", $baseLevel + $currentIncrement,".\n";
		$currentIncrement++;
	    }
	    
	    if($currentIncrement > $stackingHeight) {
		$WARNING && print "WARNING: Max level of $type stacking reached when performing level_layout for annotation $$annotation{id}->[ $sortedOfType[$i] ] between $$annotation{start}->[ $sortedOfType[$i] ] and $$annotation{stop}->[ $sortedOfType[$i] ]. Stacking on last possible level.\n";
		# $currentIncrement = 0;
		$currentIncrement --;
	    }
	    
	    $DEVEL && print "DEVEL: Placing $$annotation{id}->[ $sortedOfType[$i] ] on level ",$baseLevel + $currentIncrement,"; it starts on $$annotation{start}->[ $sortedOfType[$i] ], but last on this level ends on $lastOnLevel[$currentIncrement].\n";
	    $lastOnLevel[$currentIncrement] = $$annotation{stop}->[ $sortedOfType[$i] ];
	
	    if($$annotation{frame}->[$sortedOfType[$i]] < 3) {
		$$annotation{level}->[$sortedOfType[$i]] = $baseLevel + $currentIncrement;
	    } elsif ($$annotation{frame}->[$sortedOfType[$i]] < 6)  {
		$$annotation{level}->[$sortedOfType[$i]] =  - ($baseLevel + $currentIncrement);
	    } elsif ($$annotation{frame}->[$sortedOfType[$i]] == $AnnotationBoxFrame) {
		$$annotation{level}->[$sortedOfType[$i]] = $baseLevel + $currentIncrement;
	    } else {
		$WARNING && print "WARNING: Annotation $$annotation{id}->[ $sortedOfType[$i] ] has incorrect frame $$annotation{frame}->[$sortedOfType[$i]], and was not assigned a level in level_layout for $type.\n";
	    }
	}
    }
}

sub add_glimmer {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

#  my $loadDialog = $$main->FileSelect(-directory=>'.',-filelabel=>'Load file',-filelistlabel=>'Files');
#  $loadDialog->title("A GUI - Add Glimmer annotation from glim file");
  
#  my $selectedFile = $loadDialog->Show;

  my @filetypes=(['glim results','*.glim']);
  my $selectedFile = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.glim',-title=>"Extract Glimmer2 annotations from glim output");

  if(!$selectedFile) {
    $WARNING && print "Cancel load.\n";
  } else {
    $WARNING && print "Load glimmer annotations to $seqName from $selectedFile.\n";
    import_glimmer($main,$selectedFile,$annotation,$seqName,$seq,$sheet); # IMPROVEMENT: ugly ref handling!
  }
}

sub add_gff {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $bioframed=shift;
  my $sheet=shift;

  my @filetypes=(['gff2 results file','*.gff']);
  my $selectedFile = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.gff',-title=>"Extract GFF2 annotations from file");

  if(!$selectedFile) {
    $WARNING && print "Cancel load.\n";
  } else {
    $WARNING && print "Load gff2 annotations to $seqName from $selectedFile.\n";
    import_gff($main,$selectedFile,$bioframed,$annotation,$seqName,$seq,$sheet);
  }
}

sub add_testcode {
      my $main=shift;
      my $annotation=shift;
      my $seqName=shift;
      my $seq=shift;
      my $sheet=shift;

      my @filetypes;
      my $selectedForwardFile;
      my $selectedReverseFile;

      # viewing testcode as graph? postponed.

      if($$sheet{usebuiltintestcode}) {
	  do_testcode($main, $annotation, $seqName, $seq, $sheet);
      } else {
	  @filetypes=(['Testcode results','*'],['forward','*.f'],['reverse','*.r']);
	  $selectedForwardFile=$$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.f',-title=>"Open testcode run forward file");
	  $selectedReverseFile=$$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.r',-title=>"Open testcode run reverse file");

	  if(!($selectedForwardFile && $selectedReverseFile)) {
	      $WARNING && print "Cancel load.\n";
	  } else {
	      $WARNING && print "Load TestCode values for $seqName from $selectedForwardFile and $selectedReverseFile.\n";
	      import_testcode($main,$selectedForwardFile,$selectedReverseFile,$annotation,$seqName,$seq,$sheet);
	  }
      }
}

sub add_estorf {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift; # ref
  my $sheet=shift;
  
#  my $loadDialog = $$main->FileSelect(-directory=>'.',-filelabel=>'Load file',-filelistlabel=>'Files');
#  $loadDialog->title("A GUI - Add ESTORF annotation from pmest output");
  
#  my $selectedFile = $loadDialog->Show;

  my @filetypes=(['pmest results','*.pmest']);
  my $selectedFile = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.pmest',-title=>"Extract ESTORF annotation from pmest output");
  
  if(!$selectedFile) {
    $WARNING && print "Cancel load.\n";
  } else {
    $WARNING && print "Load ESTORF annotations to $seqName from $selectedFile.\n";
    import_estorf($main,$selectedFile,$annotation,$seqName,$seq,$sheet);
  }
}

sub add_est {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

  my @filetypes=(['pmest results','*.pmest']);
  my $selectedFile=$$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.pmest',-title=>"Extract EST hit annotation from pmest output");

  if(!$selectedFile) {
    $WARNING && print "Cancel load.\n";
  } else {
    $WARNING && print "Load EST annotations to $seqName from $selectedFile.\n";
    import_est($main,$selectedFile,$annotation,$seqName,$seq,$sheet);
  }
}

sub add_regexp {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

  my ($start, $stop, $frame);
  my $regexp;
  my $caseinsensitive = 1;

  my $annotation_nr=($$annotation{nr})?($$annotation{nr}):0; # If this is the first annotation, its 0 and not undef.
  my $ar_win=$$main->Toplevel;
  $ar_win->title("Add regexp matches");
 # $ar_win->geometry('+200+0');
  $ar_win->configure(-width=>"500", -background=>$$sheet{default_win_background});
  
  my $ar_main_frame=$ar_win->Frame(-background=>$sheet{default_win_background})->pack(-padx=>6,-pady=>6,-fill=>'both',-expand=>'yes');
  my $ar_type_label=$ar_main_frame->Label(-background=>$sheet{default_win_background}, -text=>"Type: regexp")->pack(-side=>'top',-anchor=>'w');

  my $ar_regexp_frame=$ar_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
  my $ar_regexp_label=$ar_regexp_frame->Label(-background=>$sheet{default_win_background}, -text=>"Find (perl regexp)")->pack(-side=>'left');
  my $ar_regexp=$ar_regexp_frame->Entry(-background=>$sheet{default_win_background}, -textvariable=>\$regexp, -width=>0)->pack(-side=>'left');

  my $ar_case_cb_frame = $ar_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
  my $ar_case_cb = $ar_case_cb_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$caseinsensitive,-text=>"Case insensitive")->pack(-side=>"left");

  # window, regexp search, check for
  
  # add as overlapping  (a la blast)

  my $ar_action_frame=$ar_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-side=>'top',-anchor=>'w',-expand=>'yes');
  my $ar_apply_button;
  my $ar_ok_button=$ar_action_frame->Button(-activebackground=>$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-text=>"Ok",-command=>sub {  
      $ar_apply_button->invoke();
      $ar_win->destroy;
  })->pack(-side=>'left');
  $ar_apply_button=$ar_action_frame->Button(-activebackground=>$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-text=>"Apply",-command=>sub {
      my $color='hot pink';
      # forward
      while( $caseinsensitive ? $$seq=~m/$regexp/gi : $$seq=~m/$regexp/g ) {
	  $start=$-[0]+1;
	  $stop=$+[0];

          # update frame -- forward pass
	  $frame=($start%3)?($start%3)-1:2;
	      
	  $DEBUG && print "DEBUG: hit start $-[0], hit stop $+[0] so setting forward start to $start, stop $stop and frame $frame.\n";

	  # add annotation
	  annotation_add_new($annotation,type => 'regexp', id => "ab_$annotation_nr", start => $start, stop => $stop, frame => $frame, seqName => $seqName, color => $color, name => $&, comment => "match $regexp", note => "$& matches $regexp");
	  $annotation_nr=$$annotation{nr}; # catch the updated annotation_nr - counted up in annotation_add_new

	  $$sheet{status}="dirty"; # Sheet is now dirty.
	  
	  $$main->update();
      }

      # reverse
      my $revseq=reverse(split(/ */,$$seq));
      $revseq=~tr/atgcATGC/tacgTACG/;
      my $seqlen = length($revseq);
      
      while( $caseinsensitive ? $revseq=~m/$regexp/gi : $revseq=~m/$regexp/g ) {
	  # note -- revcoords!!!   sl------*-----+----------1
	  $start = $seqlen - $+[0]+1; 
	  $stop = $seqlen - $-[0];

          # update frame -- reverse pass
	  $frame = 5-(($seqlen-$stop)%3);

	  $DEBUG && print "DEBUG: hit start $-[0], hit stop $+[0] and seqlen is $seqlen so setting reverse start to $start, stop $stop and frame $frame.\n";

	  if($start > $stop) {
	      $DEBUG && print "DEBUG: Oops - error in revcomp coord computation in add_regexp..\n";
	  }

	  annotation_add_new($annotation,type => 'regexp', id => "ab_$annotation_nr", start => $start, stop => $stop, frame => $frame, seqName => $seqName, color => $color, name => $&, comment => "match $regexp", note => "$& matches $regexp");
	  $annotation_nr = $$annotation{nr}; # catch the updated annotation_nr - counted up in annotation_add_new

	  $$sheet{status}="dirty"; # Sheet is now dirty.

	  $$main->update();
      }

      $DEBUG && print "DEBUG: Done regexp searching - layout and redraw.\n";
      level_layout($annotation,'regexp');
      redraw_annotations(\$canvas,$annotation,$sheet,$seq); # using global canvas
  })->pack(-side=>'left');
  my $ar_cancel_button=$ar_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-text=>"Cancel",-command=>sub {
      $ar_win->destroy;
  })->pack(-side=>'right');
}

sub add_selected_as_manual {
  my $main=shift;
  my $canvas = shift;
  my $annotation=shift; 
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

  # pick selected 
  # add manual at selected coords.

  $selected_id = $$sheet{selected};
  
  my $annotation_nr = annotation_what_nr($annotation,$selected_id);

  my $start = $$annotation{start}->[$annotation_nr];
  my $stop  = $$annotation{stop}->[$annotation_nr];
  my $frame = $$annotation{frame}->[$annotation_nr];
  my $name  = $$annotation{name}->[$annotation_nr];
  my $note  = $$annotation{note}->[$annotation_nr];
  my $color = "magenta";
  
  my $level;
  if($frame < 3) {
      $level = $frame;
  } elsif ($frame < $AnnotationBoxFrame) {
      $level = $frame - 5;
  } elsif ($frame == $AnnotationBoxFrame) {
      $level = 0;
  }

  my $new_annotation_nr = ($$annotation{nr})?($$annotation{nr}):0; # If this is the first annotation, its 0 and not undef. Unlikely, btw... =)

  annotation_add_new($annotation, type => 'manual',id => "ab_$new_annotation_nr", start => $start, stop => $stop, seqName => $seqName, frame => $frame, name => $name, note => $note, color => $color, level => $level);
  
  redraw_annotations($canvas,$annotation,$sheet,$seq);
}

sub add_manual {
  my $main=shift;
#  my $canvas = shift;
  my $annotation=shift; 
  my $seqName=shift;  #not ref
  my $seq=shift;
  my $sheet=shift;
  
  my $annotation_nr=($$annotation{nr})?($$annotation{nr}):0; # If this is the first annotation, its 0 and not undef.
  my ($start,$stop,$frame,$comment,$name,$note);
  my $color = "magenta";
  my ($forward,$reverse,$none);
  # edit_annotation, but with editable frame, start and stop.

  my $seq_length=length($$seq);

  my $am_win=$$main->Toplevel;
  $am_win->title("Add annotation");
 #  $am_win->geometry('+200+0');
  $am_win->configure(-width=>"500", -background=>$$sheet{default_win_background});
  
  my $am_main_frame=$am_win->Frame(-background=>$sheet{default_win_background})->pack(-padx=>6,-pady=>6,-fill=>'both',-expand=>'yes');

  my $am_type_label=$am_main_frame->Label(-background=>$sheet{default_win_background}, -text=>"Type: manual")->pack(-side=>'top',-anchor=>'w');

  my $am_start_frame=$am_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
  my $am_start_label=$am_start_frame->Label(-background=>$sheet{default_win_background}, -text=>"Start: ")->pack(-side=>'left');
  my $am_start_entry=$am_start_frame->Entry(-background=>$sheet{default_win_background}, -textvariable=>\$start)->pack(-side=>'left');

  my $am_stop_frame=$am_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
  my $am_stop_label=$am_stop_frame->Label(-background=>$sheet{default_win_background}, -text=>"Stop: ")->pack(-side=>'left');
  my $am_stop_entry=$am_stop_frame->Entry(-background=>$sheet{default_win_background}, -textvariable=>\$stop)->pack(-side=>'left');

  my $am_frame_frame=$am_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
  my $am_frame_label=$am_frame_frame->Label(-background=>$sheet{default_win_background},-text=>"Frame: ")->pack(-side=>'left'); # spam, spam, frame and spam, please
  my $am_frame_radio_frame=$am_frame_frame->Frame(-background=>$sheet{default_win_background})->pack(-side=>'left',-anchor=>'w');
  
  my $am_frame_forward_radio=$am_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background}, -command=>sub {
      if($forward) {
	  $none=$reverse=0;
      }
      $frame=bioframe( (($start%3)?($start%3)-1:2) , 'manual')  ; # 0, 1 or 2.      
  },-variable=>\$forward,-value=>'1',-text=>'Forward')->pack(-side=>'top',-anchor=>'w');
  my $am_frame_reverse_radio=$am_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-command=>sub {
      if($reverse) {
	  $none=$forward=0;
      }
      $frame=bioframe( (5-(($seq_length-$stop)%3)), 'manual'); # 3, 4 or 5.
  },-variable=>\$reverse,-value=>'1',-text=>'Reverse')->pack(-side=>'top',-anchor=>'w');
  my $am_frame_none_radio=$am_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-command=>sub {
      if($none) {
	  $reverse=$forward=0;
      }
      $frame=bioframe($AnnotationBoxFrame,'manual'); # ANY
  },-text=>'None',-variable=>\$none,-value=>'1')->pack(-side=>'top',-anchor=>'w');
  my $am_frame_frame_label=$am_frame_frame->Label(-background=>$sheet{default_win_background},-textvariable=>\$frame)->pack(-side=>'left');
  $am_frame_none_radio->invoke();
  
  my $am_name_heading=$am_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation name")->pack(-pady=>6,-side=>'top',-anchor=>'w');
  my $am_name=$am_main_frame->Entry(-background=>$sheet{default_win_background},-textvariable=>\$name,-justify=>'left')->pack(-side=>'top',-anchor=>'w',-fill=>'x');
  my $am_comment_heading=$am_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotator comment")->pack(-pady=>6,-side=>'top',-anchor=>'w');
  my $am_comment=$am_main_frame->Entry(-background=>$sheet{default_win_background},-textvariable=>\$comment,-justify=>'left')->pack(-side=>'top',-anchor=>'w',-fill=>'x');
  my $am_notes_heading=$am_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation notes")->pack(-pady=>6,-side=>'top',-anchor=>'w');
  my $am_notes_text=$am_main_frame->Text(-background=>$sheet{default_win_background})->pack(-pady=>6,-side=>'top',-anchor=>'w',-fill=>'x');

  my $am_color_frame=$am_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-pady=>6,-side=>'top',-anchor=>'w');
  my $am_color_heading=$am_color_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation color")->pack(-side=>'left');
  my $am_color_button;
  $am_color_button = $am_color_frame->Button(-activebackground=>$color, -background=>$color,-text=>$color, -command=>sub {
      $color = $am_win->chooseColor(-title => "Choose an annotation color", -initialcolor => $color);
      $am_color_button->configure(-activebackground => $color, -background => $color,-text => $color);
  })->pack(-side=>'left');

  my $am_action_frame=$am_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-side=>'top',-anchor=>'w',-expand=>'yes');
  my $am_apply_button;
  my $am_ok_button=$am_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-text=>"Ok",-command=>sub {
      # sanity checks
      if($start>$stop && ($reverse or $none)) {
	  $_=$stop;
	  $stop=$start;
	  $start=$_;
      } elsif ($start>$stop && $forward) {
	  $am_win->messageBox(-icon=>'error',-type=>'OK',-message=>"Start is further upstream than stop, and orientation is set to forward.");
	  return 0;
      }
      if((!$start)||(!$stop)) {
	  $am_win->messageBox(-icon=>'error',-type=>'OK',-message=>"In order to add this annotation, you must specify start and stop positions. (1st pos is 1, last is $seq_length");
	  return 0;
      }
      
      # update frame...
      if($forward) {
	  $frame=($start%3)?($start%3)-1:2;
      } elsif ($reverse) {
	  $frame=5-(($seq_length-$stop)%3);
      } elsif ($none) {
	  $frame=6;
      }
      
      my $level;
      if($frame < 3) {
	  $level = $frame;
      } elsif ($frame < $AnnotationBoxFrame) {
	  $level = $frame - 5;
      } elsif ($frame == $AnnotationBoxFrame) {
	  $level = 0;
      }
      
      $note = $am_notes_text->get('1.0','end');
      
      annotation_add_new($annotation, type => 'manual',id => "ab_$annotation_nr", start => $start, stop => $stop, seqName => $seqName, frame => $frame, name => $name, note => $note, color => $color, level => $level);
      
      # prepare for next annotation
      $annotation_nr = $$annotation{nr};
      
      $$sheet{status}="dirty"; # Sheet is now dirty.
      
      redraw_annotations(\$canvas,$annotation,$sheet,$seq); # canvas is a global..

      $am_win->destroy;
  })->pack(-side=>'left');
  my $am_cancel_button=$am_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-text=>"Cancel",-command=>sub {
      $am_win->destroy;
  })->pack(-side=>'right');
}

sub add_orf {
    my $main=shift;
    my $annotation=shift;
    my $seqName=shift;
    my $seq=shift;
    my $sheet=shift;
    
    # DEBUG: bend, it should.
    my $shortestOrf = 15; #150

    # Translate trinucleuotides -> aa in all six reading frames.
    my @phase;
    my @phaseOffset = (0,1,2,2,1,0);
    
    $phase[0] = nt2aa($$seq);
    $phase[1] = nt2aa(substr($$seq,1));
    $phase[2] = nt2aa(substr($$seq,2));

    # Same story with the complementary strand.
    my $complementSeq = $$seq;
    # The complement has complementary bases, and in the
    # translation perspective, reverse orientation.
    $complementSeq=~tr/atgcATGC/tacgTACG/;
    $complementSeq = join('',reverse(split(/ */,$complementSeq)));
    $phase[3]=nt2aa(substr($complementSeq,2));
    $phase[4]=nt2aa(substr($complementSeq,1));
    $phase[5]=nt2aa($complementSeq);

    for(my $i=0;$i<6;$i++) {
	my @orf=();
	while( $phase[$i] =~ m/((?:M\w*\+)|(?:^\w*\+)|(?:M\w*$))/g ) { # ignoring the (unlikely) event that the entire seq is an orf with no start or stop in - that would be painfully obvious from showing ST..
	    push @orf,$1;
	    $DEVEL && $DEBUG && print "DEBUG: longest hits from phase $i: $1\n";
	}
	my $startAaPos;
	my $stopAaPos;
	my $startNtPos;
	my $stopNtPos;
	       
	foreach (@orf) {
	    my $orflen = length($_);
	    if(3 * $orflen > $shortestOrf) {
		$startAaPos = index($phase[$i],$_);
		$startAaPos != -1 || ($WARNING && print "WARNING: Index on phase $i failed for $_.\n");
		$stopAaPos=$startAaPos+$orflen;
		if($i<=2) {
		    $startNtPos = 3*$startAaPos+$phaseOffset[$i]+1; # to 1-based nt-pos
		    $stopNtPos = 3*$stopAaPos+$phaseOffset[$i]; # to 1-based nt-pos
		} else {
		    $startNtPos = -3*$stopAaPos+length($$seq)+1-$phaseOffset[$i]; # to 1-based nt-pos
		    $stopNtPos = -3*$startAaPos+length($$seq)-$phaseOffset[$i]; # to 1-based nt-pos, inverted start & stop always to have lowest number as NtStart.
		}

		# Get current top annotation nr, or initialize it if this is the first
		# annotation made.
		if(! defined($$annotation{nr})) {
		    $$annotation{nr} = 0;
		}
		
		my $nr = $$annotation{nr};

		my $frame = $i;
		my $color = 'gray';

		my $level;
		if($i < 3) {
		    $level = $frame;
		} else {
		    $level = $frame - 5;
		}

		$DEVEL && ($$annotation{name}->[$nr] = "ab_$nr");
		
		annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type => 'ORF', start => $startNtPos, stop => $stopNtPos, frame => $frame, color => $color, comment => "", note => "", level => $level);

		$$sheet{status} = "dirty"; # Sheet is now dirty.
	    }
	}
	$$main->update();
    }
    redraw_annotations(\$canvas,$annotation,$sheet,$seq); # using global canvas!
}

sub add_polypy {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

  my $polyLength = 7;

  my $nr;
  if($$annotation{nr}) {
      $nr=$$annotation{nr};
  } else {
      $nr=0;
      $$annotation{nr}=0;
  }

  my $frame = 0;
  my $polyPy = 0;
  my $revPolyPy = 0;
  
  $DEBUG && print "DEBUG: PolyPy search on frame $frame\n";
  for(my $i=0;$i<length($$seq) ; $i++) {
      my $thisbase = lc substr($$seq,$i,1);
      if($thisbase eq 'c' || $thisbase eq 't') { # [CcTt]{3,}[AaGg]{,2}[CcTt]{3,}[AaGg]{,2}[CcTt]{3,}
	  $polyPy++;
	  if($revPolyPy) {
	      if( $revPolyPy > $polyLength ) {
		  $DEBUG && print "DEBUG: Found polyPy on reverse strand ",$i-$revPolyPy+1," - ", $i+1,"\n" ;
		  # add annotation
		  my $start = $i-$revPolyPy+1;
		  my $stop = $i; # +1 -1
		  my $add_frame = 3; #frame assignment "wrong" -- could assign correct internal frames instead
		  my $color = "green3";

		  annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type => 'polypy', start => $start, stop => $stop, frame => $add_frame, comment => "", note => "", color => $color, level => 0);
		  $nr = $$annotation{nr};

		  $$sheet{status}="dirty"; # Sheet is now dirty. 
	      }
	      $revPolyPy = 0;
	  }
      } elsif ($thisbase eq 'a' || $thisbase eq 'g') {
	  $revPolyPy++;
	  if($polyPy) {
	      if ( $polyPy > $polyLength ) {
		  $DEBUG && print "DEBUG: Found polyPy on forward strand ",$i-$polyPy+1," - ", $i+1,"\n" ;
		  # add annotation
		  my $start = $i+1-$polyPy;;
		  my $stop = $i; # +1 -1 
		  my $add_frame = 0;
		  my $color = "green3";
		  
		  annotation_add_new($annotation, id => "ab_$nr", seqName => $seqName, type => 'polypy', start => $start, stop => $stop, frame => $add_frame, comment => "", note => "", color => $color, level => 0);
		  $nr = $$annotation{nr};

		  $$sheet{status}="dirty"; # Sheet is now dirty. 
	      }	
	      $polyPy = 0;
	  }
      } elsif ($thisbase eq 'n' || $thisbase eq 'x') { # Count x and n as "hits"
	  $polyPy++;
	  $revPolyPy++;
      }
  }
  redraw_annotations(\$canvas,$annotation,$sheet,$seq);
}

sub add_st {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

  # Design question: do we want a separate entry for each start/stop?
  # One per codon and frame? 

  # Translate trinucleuotides -> aa in all six reading frames.
  my @phase;
  my @phaseOffset=(0,1,2,2,1,0);
  my @phaseEndOffset;
  
  $phase[0] = nt2aa($$seq);
  $phase[1] = nt2aa(substr($$seq,1));
  $phase[2] = nt2aa(substr($$seq,2));
  # Same story with the complementary strand.
  my $complementSeq = $$seq;
  # The complement has complementary bases, and in the 
  # translation perspective, reverse orientation.
  $complementSeq=~tr/atgcATGC/tacgTACG/;
  $complementSeq = join('',reverse(split(/ */,$complementSeq)));
  $phase[3] = nt2aa(substr($complementSeq,2));
  $phase[4] = nt2aa(substr($complementSeq,1));
  $phase[5] = nt2aa($complementSeq);

  my $nr;
  if($$annotation{nr}) {
    $nr=$$annotation{nr};
  } else {
    $nr=0;
    $$annotation{nr}=0;
  }

  my $seqlen = length($$seq);

  for(my $i=0;$i<6;$i++) {

    my $startPos=0;
    my @starts=();
    while($startPos!=-1) {
#     $DEBUG && print "DEBUG: start pos = $startPos\n";
      $startPos=index($phase[$i],'M',$startPos);
      if($startPos!=-1) {
	  
	  # Convert 0-based aa-positions to 1-based nt-positions.
	  if($i<=2) {
	      $ntStartPos=3*$startPos+$phaseOffset[$i] +1;
	  } else {
	      $ntStartPos=-3*$startPos+$seqlen-$phaseOffset[$i];
	  }
	  push @starts,$ntStartPos;
	  # Increase to keep index going..
	  $startPos++;
      }
      $$main->update();
    }

  $DEBUG && print "DEBUG: in frame $i, add_st found starts at @starts\n";

    my @stops=();
    my $stopPos=0;
    while($stopPos!=-1) {
 #     $DEBUG && print "DEBUG: stop pos = $stopPos\n";
	$stopPos=index($phase[$i],'+',$stopPos);
	if($stopPos!=-1) {
	    # Nah, I think we'd better use nt-positions instead... =)
	    # Convert 0-based aa-positions to 1-based nt-positions.
	    if($i<=2) {
		$ntStopPos=3*($stopPos+1)+$phaseOffset[$i];
	    } else {
		$ntStopPos=-3*($stopPos+1)+$seqlen+1-$phaseOffset[$i];
	    }
	    push @stops,$ntStopPos;
	    
	    # Increase to keep index going..
	    $stopPos++;
	}
    }

    $DEBUG && print "DEBUG: in frame $i, add_st found stops at @stops\n";
 
    # Actually add annontation

    $$annotation{uid}->[$nr+$i] = $$annotation{unr};
    $$annotation{unr}++;

    $$annotation{id}->[$nr+$i]="st_".($nr+$i);
    $annotation_nr_cache{$$annotation{id}->[$nr+$i]} = $nr+$i; # update cache
    
    $$annotation{type}->[$nr+$i]='ST';
    $$annotation{frame}->[$nr+$i]=$i;

    # Ugly, but functional..
    my $j=0;
    foreach $start (@starts) {
      $$annotation{start}->[$nr+$i]->[$j]=$start;
      $j++;
    }

    $j=0;
    foreach $stop (@stops) {
      $$annotation{stop}->[$nr+$i]->[$j]=$stop;
      $j++;
    }

    if($$annotation{frame}->[$nr+$i] < 3) {
      $$annotation{level}->[$nr+$i] = $$annotation{frame}->[$nr+$i];    
    } else {
      $$annotation{level}->[$nr+$i] = $$annotation{frame}->[$nr+$i] - 5;
    }

    $$annotation{nr}+=1;
   
    draw_annotation_st(\$canvas,$annotation,$nr+$i,$sheet,$seq); # using global canvas
  }
}

sub add_blast {
  my $main=shift;
  my $annotation=shift;
  my $seqName=shift;
  my $seq=shift;
  my $sheet=shift;

#  my $selectedFile = $loadDialog->Show;

  my @filetypes=(['blast results','*.html']);
  my $selectedFile = $$main->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '*.html',-title=>"Extract blast hits from saved NCBI blast response");

  if(!$selectedFile) {
    print "Cancel load.\n";
  } else {
    $DEBUG && print "DEBUG: Load blasthits to $seqName from $selectedFile.\n";
    my $blastresults = "";
    open BLASTFILE, $selectedFile;
    while(<BLASTFILE>) {
	$blastresults .= $_;
    }
    close BLASTFILE;

    import_blast($main,$blastresults,$annotation,$seqName,$seq,$sheet);
  }
}

sub annotation_add_new {
    my $annotation = shift;  # ref
    my %content = (@_); 

    my $annotation_nr=($$annotation{nr})?($$annotation{nr}):0; # If this is the first annotation, its 0 and not undef.
    $$annotation{uid}->[$annotation_nr] = $$annotation{unr};
    $$annotation{unr}++;     # prepare for next annotation

    $DEBUG && print "Adding new annotation; nr $annotation_nr, uid = $$annotation{uid}->[$annotation_nr], ";

    foreach $key (keys %content) {
	# any content checks?
	$$annotation{$key}->[$annotation_nr]=$content{$key};
	$DEBUG && print "$key = $content{$key}, ";
    }
    if($$annotation{id}->[$annotation_nr]) {
	$annotation_nr_cache{$$annotation{id}->[$annotation_nr]}=$annotation_nr; # cache
    } else {
	$WARNING && print "WARNING: no annotation id given for annotation $annotation_nr - reluctantly adding anyway!\n";
    }
    # check for vital missing fields
    if(!defined($$annotation{color}->[$annotation_nr]) or $$annotation{color}->[$annotation_nr] eq "") {
	if ($DEBUG) {
	    print "NO COLOR SET - default to black, ";
	} elsif ($WARNING) { 
	    print "WARNING: no color chosen for annotation type $$annotation{type}->[$annotation_nr] in annotation_add_new - set to black for now..\n";
	}
	$$annotation{color}->[$annotation_nr]="brown";
    }

    # six/0 frame display 
    if($$annotation{type}->[$annotation_nr] eq 'manual') {
	if($$annotation{frame}->[$annotation_nr] < 3) {
	    $$annotation{level}->[$annotation_nr]=$$annotation{frame}->[$annotation_nr];
	} elsif ($$annotation{frame}->[$annotation_nr] < $AnnotationBoxFrame) {
	    $$annotation{level}->[$annotation_nr]=$$annotation{frame}->[$annotation_nr]-5;
	} elsif ($$annotation{frame}->[$annotation_nr] == $AnnotationBoxFrame) {
	    $$annotation{level}->[$annotation_nr]=0;
	}
    } elsif($$annotation{type}->[$annotation_nr] eq "glimmer2" or $$annotation{type}->[$annotation_nr] eq "testcode" or $$annotation{type}->[$annotation_nr] eq 'polypy' or $$annotation{type}->[$annotation_nr] eq 'ORF') {
	# level already set on call; ESTORF, merge and  does currently not use this routine
    } elsif ($annotatorBioframed{$$annotation{type}->[$annotation_nr]} == 0) {
        # eq "regexp" or $$annotation{type}->[$annotation_nr] eq "blast" or $$annotation{type}->[$annotation_nr] eq "EST"
	$DEBUG && print STDERR "DEBUG: no level for type $$annotation{type}->[$annotation_nr] - awaiting level_layout\n";
    } else {
	$$annotation{level}->[$annotation_nr]=0;
	if ($DEBUG) {
	    print "NO FRAME RULE - setting default frame; ";
	} elsif ($WARNING) {
	    print "WARNING: no frame rule for annotation type $$annotation{type}->[$annotation_nr] in annotation_add_new. Adding annotation nr $annotation_nr anyway..\n";
	}
    }
    
    $DEBUG && print "frame = $$annotation{frame}->[$annotation_nr].\n";
    $$annotation{nr}++;     # prepare for next annotation

    # Please remember:
    # # layout after all annotations of one kind is in place,
    # # redraw annotations if neccessary
    # # update local copy of current annotation_nr
}

sub annotation_edit {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $seq = shift;
    my $sheet=shift;
    my $annotation_nr=shift;
    my $new_note=shift;

    my $seq_length=length($$seq);

    my $ae_win=$$main->Toplevel;
    $ae_win->title("Edit annotation $$annotation{start}->[$annotation_nr]-$$annotation{stop}->[$annotation_nr] ($$annotation{id}->[$annotation_nr])");
 #    $ae_win->geometry('+200+0');
    $ae_win->configure(-background=>$$sheet{default_win_background},-width=>"600");
  
    my $annotation_id = $$annotation{id}->[$annotation_nr];
    my $annotation_uid = $$annotation{uid}->[$annotation_nr];

    my $name = $$annotation{name}->[$annotation_nr];
    my $start = $$annotation{start}->[$annotation_nr];
    my $stop = $$annotation{stop}->[$annotation_nr];
    my $frame = $$annotation{frame}->[$annotation_nr];
    my $color = $$annotation{color}->[$annotation_nr];

    my $ae_main_frame=$ae_win->Frame(-background=>$sheet{default_win_background})->pack(-padx=>6,-pady=>6,-fill=>'both',-expand=>'yes');
  
    my $ae_type_label=$ae_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Type: $$annotation{type}->[$annotation_nr]")->pack(-side=>'top',-anchor=>'w');
    my $ae_loc_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w',-fill=>'x');

    my $ae_start_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
    my $ae_start_label=$ae_start_frame->Label(-background=>$sheet{default_win_background}, -text=>"Start: ")->pack(-side=>'left');
    my $ae_start_entry=$ae_start_frame->Entry(-background=>$sheet{default_win_background}, -textvariable=>\$start)->pack(-side=>'left');

    my $ae_stop_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
    my $ae_stop_label=$ae_stop_frame->Label(-background=>$sheet{default_win_background}, -text=>"Stop: ")->pack(-side=>'left');
    my $ae_stop_entry=$ae_stop_frame->Entry(-background=>$sheet{default_win_background}, -textvariable=>\$stop)->pack(-side=>'left');

    my $ae_frame_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-expand=>'yes',-side=>'top',-anchor=>'w');
    my $ae_frame_label=$ae_frame_frame->Label(-background=>$sheet{default_win_background},-text=>"Frame: ")->pack(-side=>'left'); # spam, spam, frame and spam, please
    my $ae_frame_radio_frame=$ae_frame_frame->Frame(-background=>$sheet{default_win_background})->pack(-side=>'left',-anchor=>'w');
    
    my $ae_frame_forward_radio=$ae_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background}, -command=>sub {
								  if($forward) {
								    $none=$reverse=0;
								  }
								  $frame=bioframe( (($start%3)?($start%3)-1:2) , 'manual')  ; # 0, 1 or 2.
								  
							      },-variable=>\$forward,-value=>'1',-text=>'Forward')->pack(-side=>'top',-anchor=>'w');
    my $ae_frame_reverse_radio=$ae_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-command=>sub {
								  if($reverse) {
								    $none=$forward=0;
								  }
								  $frame=bioframe( (5-(($seq_length-$stop)%3)), 'manual'); # 3, 4 or 5.
								},-variable=>\$reverse,-value=>'1',-text=>'Reverse')->pack(-side=>'top',-anchor=>'w');
  my $ae_frame_none_radio=$ae_frame_radio_frame->Radiobutton(-activebackground=>$$sheet{dialog_button_active_background}, -background=>$sheet{default_win_background},-command=>sub {
								  if($none) {
								    $reverse=$forward=0;
								  }
								  $frame=bioframe($AnnotationBoxFrame,'manual'); # ANY
								},-text=>'None',-variable=>\$none,-value=>'1')->pack(-side=>'top',-anchor=>'w');
    my $ae_frame_frame_label=$ae_frame_frame->Label(-background=>$sheet{default_win_background},-textvariable=>\$frame)->pack(-side=>'left');
    
    if($frame < 3) {
	$ae_frame_forward_radio->invoke();
    } elsif ($frame < 6) {
	$ae_frame_reverse_radio->invoke();
    } elsif ($frame == 6) {
	$ae_frame_none_radio->invoke();
    }

    if($$annotation{type}->[$annotation_nr] ne 'manual') {
	$ae_start_entry->configure(-state=>'disabled');
	$ae_stop_entry->configure(-state=>'disabled');
	$ae_frame_forward_radio->configure(-state=>'disabled');
	$ae_frame_reverse_radio->configure(-state=>'disabled');
	$ae_frame_none_radio->configure(-state=>'disabled');
    }

    # name
    my $ae_name_label=$ae_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation name")->pack(-pady=>6,-side=>'top',-anchor=>'w');
    my $ae_name_entry=$ae_main_frame->Entry(-background=>$sheet{default_win_background},-textvariable=>\$name)->pack(-pady=>6,-side=>'top',-anchor=>'w',-fill=>'both',-expand=>'yes');
    my $ae_comment_heading=$ae_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotator comment")->pack(-pady=>6,-side=>'top',-anchor=>'w');

    my $ae_comment_label = "";
    defined($$annotation{comment}->[$annotation_nr]) && ($ae_comment_label = $$annotation{comment}->[$annotation_nr]);
    my $ae_comment=$ae_main_frame->Label(-relief=>'sunken',-background=>$sheet{default_win_background},-text=>"$ae_comment_label",-justify=>'left', -wraplength=>700)->pack(-side=>'top',-anchor=>'w',-fill=>'x');
    my $ae_notes_heading=$ae_main_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation notes")->pack(-pady=>6,-side=>'top',-anchor=>'w');
#  my $ae_notes_entry=$ae_main_frame->Entry(-background=>$sheet{default_win_background},-textvariable=>\$new_note)->pack(-pady=>6,-side=>'top',-anchor=>'w',-fill=>'both',-expand=>'yes');
    my $ae_notes_text=$ae_main_frame->Text(-background=>$sheet{default_win_background})->pack(-pady=>6,-side=>'top',-anchor=>'w',-fill=>'x');
    $ae_notes_text->insert('end',$new_note);

    my $ae_color_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-pady=>6,-side=>'top',-anchor=>'w');
    my $ae_color_heading=$ae_color_frame->Label(-background=>$sheet{default_win_background},-text=>"Annotation color")->pack(-side=>'left');

    my $ae_color_button;
    $ae_color_button = $ae_color_frame->Button(-activebackground=>$color, -background=>$color,-text=>$color, -command=>sub {
	$color = $ae_win->chooseColor(-title => "Choose an annotation color", -initialcolor => $color);
	$ae_color_button->configure(-activebackground => $color, -background => $color,-text => $color);
    })->pack(-side=>'left');

    my $ae_action_frame=$ae_main_frame->Frame(-background=>$sheet{default_win_background})->pack(-fill=>'x',-side=>'top',-anchor=>'w',-expand=>'yes');
    my $ae_apply_button;
    my $ae_ok_button=$ae_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$sheet{default_win_background},-text=>"Ok",-command=>sub {
	$ae_apply_button->invoke();
	$ae_win->destroy;
    })->pack(-side=>'left');

    $ae_apply_button = $ae_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$sheet{default_win_background},-text=>"Apply",-command=>sub {
      # Check annotation nr against annotation uid before saving - subsequent deletes have rendered annotation_nr obsolete...
	my $nr = annotation_what_nr_for_uid($annotation, $annotation_uid);
	if($nr == -1) {
	    $$main->messageBox(-icon=>'error',-type=>'OK',-message=>"Cannot find annotation to edit.\n\nHas it perhaps been removed after the opening of this window?");
	} else {
	    $DEVEL && ($annotation_nr != $nr) && print "DEVEL: Annotation_nr changed since edit window opened.. Obtained current nr $nr for uid $annotation_uid.\n";
	    $annotation_nr = $nr;

	    $$annotation{note}->[$annotation_nr] = $ae_notes_text->get('1.0','end');
	    $$annotation{name}->[$annotation_nr] = $name;      
	    $$annotation{color}->[$annotation_nr] = $color;

	    if($$annotation{type}->[$annotation_nr] eq 'manual') {
		# sanity checks
		if($start>$stop && ($reverse or $none)) {
		    $_=$stop;
		    $stop=$start;
		    $start=$_;
		} elsif ($start>$stop && $forward) {
		    $ae_win->messageBox(-icon=>'error',-type=>'OK',-message=>"Start is further upstream than stop, and orientation is set to forward.");
		    return 0;
		}
		if((!$start)||(!$stop)) {
		    $ae_win->messageBox(-icon=>'error',-type=>'OK',-message=>"In order to add this annotation, you must specify start and stop positions. (1st pos is 1, last is $seq_length");
		    return 0;
		}
		
		# update frame...
		if($forward) {
		    $frame=($start%3)?($start%3)-1:2;
		} elsif ($reverse) {
		    $frame=5-(($seq_length-$stop)%3);
		} elsif ($none) {
		    $frame=6;
		}
		
		$$annotation{start}->[$annotation_nr]=$start;
		$$annotation{stop}->[$annotation_nr]=$stop;
		$$annotation{frame}->[$annotation_nr]=$frame;
	    }

	    # refresh
	    redraw_annotations($canvas,$annotation,$sheet,$seq);
	    $$sheet{status}="dirty";
	}
    })->pack(-side=>'left');
    my $ae_cancel_button=$ae_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$sheet{default_win_background},-text=>"Cancel",-command=>sub {
	$ae_win->destroy;
    })->pack(-side=>'right');
}

# View menu

sub view_zoomOut {
  my $main=shift;
  my $canvas=shift;
  my $annotation=shift;
  my $sheet=shift;
  my $plotstate=shift;
  my $seqName=shift;
  my $seq=shift;

  if($$sheet{zoom} eq 'sequence') {
    $$sheet{zoom}='normal';
  } elsif($$sheet{zoom} eq 'normal') {
      $$sheet{zoom}='overview';
  } elsif($$sheet{zoom} eq 'overview') {
      $$sheet{zoom}='birdseye';
  }  

  # redraw seq...
  draw_seq($canvas,$$seqName,$seq,$sheet);

  redraw_annotations($canvas,$annotation,$sheet,$seq);

  if(($$sheet{zoom} ne 'sequence') && $$plotstate{axis}) {
      # redraw graph... *SLOW*
      $DEBUG && print "DEBUG: deleting graph....\n";
      $$canvas->delete(-tags=>"graph");
      $$plotstate{axis}=0;
      replot($canvas,$seq,$plotstate,$sheet);
#      calculate_graph($main,$canvas,$seq,$plotstate,$sheet);
  }
}

sub view_zoomIn {
  my $main=shift;
  my $canvas=shift;
  my $annotation=shift;
  my $sheet=shift;
  my $plotstate=shift;
  my $seqName=shift;
  my $seq=shift;

  if($$sheet{zoom} eq 'normal') {
      $$sheet{zoom}='sequence';
      if($$plotstate{axis}) {
	  # remove graph and shrink canvas.
	  $DEBUG && print "DEBUG: deleting graph....\n";
	  $$canvas->delete(-tags=>"graph");
	  $$canvas->configure(-height=>$$sheet{canvas_seq_height});
      }
  } elsif($$sheet{zoom} eq 'overview' ) {
      $$sheet{zoom}='normal';
  } elsif($$sheet{zoom} eq 'birdseye') {
      $$sheet{zoom}='overview';
  }

  draw_seq($canvas,$$seqName,$seq,$sheet);
  redraw_annotations($canvas,$annotation,$sheet,$seq);
  if($$sheet{zoom} ne 'sequence' && $$plotstate{axis}) {
      $$canvas->delete(-tags=>"graph");
      $DEBUG && print "DEBUG: reintroducing graph...\n";
      $$plotstate{axis}=0; # Tell calculate_graph to redraw axis...
      replot($canvas,$seq,$plotstate,$sheet);
#      calculate_graph($main,$canvas,$seq,$plotstate,$sheet);
  }
}

sub view_zoom {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $plotstate=shift;
    my $seqName=shift;
    my $seq=shift;

    my $zoom_birdseye = 0;
    my $zoom_overview = 0;
    my $zoom_normal = 0;
    my $zoom_sequence = 0;
    my $zoom_factor = 0;
    my $zoom_factor_value = 20;

    eval("\$zoom_$$sheet{zoom} = 1");
    if($zoom_factor) {
	$zoom_factor_value = $$sheet{zoomfactor};
    }

    my $zoom_factor_entry;

    # popup zoom selection window..
    # zoom dialog..
    my $zoom_win=$$main->Toplevel;
    $zoom_win->title("A GUI - Zoom settings");
    # $zoom_win->geometry('+250+0');
    $zoom_win->configure(-background=>'linen');
    
    my $zoom_frame=$zoom_win->Frame(-background=>$$sheet{default_win_background})->pack;

    my $zoom_label=$zoom_frame->Label(-background=>$$sheet{default_win_background},-text=>"Zoom level")->pack(-anchor=>"w",-side=>"top");
    my $zoom_toggle_frame=$zoom_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w');

    my $zoom_birdseye_toggle=$zoom_toggle_frame->Radiobutton(-command=> sub {
	if($zoom_birdseye) {
	    $zoom_overview=0;
	    $zoom_normal=0;
	    $zoom_sequence=0;
	    $zoom_factor=0;
	    $zoom_factor_entry->configure(-state=>'disabled');	  
	}      
    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$zoom_birdseye, -value=>'1',-text=>"Birdseye (40)")->pack(-side=>"top",-anchor=>"w");
    my $zoom_overview_toggle=$zoom_toggle_frame->Radiobutton(-command=> sub {
	if($zoom_overview) {
	    $zoom_birdseye=0;
	    $zoom_normal=0;
	    $zoom_sequence=0;
	    $zoom_factor=0;
	    $zoom_factor_entry->configure(-state=>'disabled');
	}

    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$zoom_overview,-value=>'1',-text=>"Overview (5)")->pack(-side=>"top",-anchor=>"w");
    my $zoom_normal_toggle=$zoom_toggle_frame->Radiobutton(-command=> sub {
	if($zoom_normal) {
	    $zoom_birdseye=0;
	    $zoom_overview=0;
	    $zoom_sequence=0;
	    $zoom_factor=0;
	    $zoom_factor_entry->configure(-state=>'disabled');
	}

    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$zoom_normal,-value=>'1',-text=>"Normal (1)")->pack(-side=>"top",-anchor=>"w");
    my $zoom_sequence_toggle=$zoom_toggle_frame->Radiobutton(-command=> sub {
	if($zoom_sequence) {
	    $zoom_birdseye=0;
	    $zoom_overview=0;
	    $zoom_normal=0;
	    $zoom_factor=0;
	    $zoom_factor_entry->configure(-state=>'disabled');
	}
    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$zoom_sequence,-value=>'1',-text=>"Sequence")->pack(-side=>"top",-anchor=>"w");

    my $zoom_factor_frame=$zoom_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>"w");  
    my $zoom_factor_toggle=$zoom_factor_frame->Radiobutton(-command=> sub {
	if($zoom_factor) {
	    $zoom_birdseye=0;
	    $zoom_overview=0;
	    $zoom_normal=0;
	    $zoom_sequence=0;
	    $zoom_factor_entry->configure(-state=>'normal');
	}
    }, -activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$zoom_factor,-value=>'1',-text=>"Zoom factor: ")->pack(-anchor=>'n', -side=>"left");
    $zoom_factor_entry=$zoom_factor_frame->Entry(-state=>'disabled', -background=>$$sheet{default_win_background},-textvariable=>\$zoom_factor_value)->pack(-anchor=>'n', -side=>"left");

    
    
#  my $zoom_info_frame=$zoom_frame->Frame(-background=>$sheet{default_win_background})->pack(-side=>"top");

    my $zoom_action_frame=$zoom_frame->Frame(-background=>$$sheet{default_win_background})->pack(-anchor=>"w",-side=>"top",-fill=>"x",-pady=>5,-ipady=>3);
    my $zoom_button_apply;
    my $zoom_button_ok=$zoom_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Ok",-command=>sub {
	$zoom_button_apply->invoke();
	$zoom_win->destroy;
    })->pack(-side=>"left");
    
    $zoom_button_apply=$zoom_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Apply",-command=>sub {
	if($$sheet{'zoom'} ne 'sequence' && $$plotstate{axis}) {
	    # remove graph and shrink canvas.
	    $DEBUG && print "DEBUG: deleting graph....\n";
	    $$canvas->delete(-tags=>"graph");
	    $$canvas->configure(-height=>$$sheet{canvas_seq_height});
	}

	if($zoom_birdseye) {
	    $$sheet{zoom}='birdseye';
	} elsif($zoom_overview) {
	    $$sheet{zoom}='overview';
	} elsif($zoom_normal) {
	    $$sheet{zoom}='normal';
	} elsif($zoom_sequence) {
	    $$sheet{zoom}='sequence';	    
	} elsif($zoom_factor) {
	    $$sheet{zoom}='factor';
	    $$sheet{zoomfactor} = $zoom_factor_value;
	} else {
	    $WARNING && print "WARNING: No zoom level choosen - this should never happen!\n";
	}

	my @old_scroll_pos = $$canvas->xview(); # save scroll center position
	my $old_scroll_center = $old_scroll_pos[0]+($old_scroll_pos[1]-$old_scroll_pos[0])/2;

	view_refresh($main,$canvas,$annotation,$seqName,$seq,$plotstate,$sheet);

	my @current_scroll_pos = $$canvas->xview(); # calculate new scroll width
	my $current_scroll_width = ($current_scroll_pos[1] - $current_scroll_pos[0]);
	
	$$canvas->xview(moveto=>$old_scroll_center-($current_scroll_width/2)); # scroll the canvas - it controls the scrollbar in turn!
    })->pack(-side=>"left");
    my $zoom_button_cancel=$zoom_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub {$zoom_win->destroy;})->pack(-side=>"right");
}

sub view_levels {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $seqName=shift;
    my $seq=shift;
    my $plotstatus = shift;
    my $sheet=shift;

    my $level_win = $$main->Toplevel();
    $level_win->title("A GUI - display levels");
    $level_win->configure(-background=>'linen');
    
    my $level_main_frame = $level_win->Frame(-background=>$$sheet{default_win_background})->pack;
    
    # toggles

    my $level_check_frame = $level_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top', -anchor=>'w', -pady=>5);
    
   # sheet{"display_${key}_level"} = 1;

    foreach my $sheetkey (keys %$sheet) {
	if($sheetkey =~ m/^display_\w+_level/) {
	    my ($display_sheetkey) = $sheetkey=~m/^display_(\w+)_level$/;
	    if(!defined($annotatorName{$display_sheetkey})) {
		$WARNING && print "WARNING: Problem adding ".$annotatorName{$display_sheetkey}." on key $display_sheetkey from $sheetkey - perhaps one entry in a preferences file is incorrectly spelled or formatted?\n";
	    } else {
		my $this_annotator_frame = $level_check_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top', -anchor=>'w', -fill=>'x');
		$this_annotator_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$sheet{$sheetkey},-text=>"$annotatorName{$display_sheetkey}")->pack(-side=>"left");
		$annotatorBioframed{$display_sheetkey} || $this_annotator_frame->Entry(-background=>$$sheet{default_win_background}, -textvariable=>\$annotatorHeight{$display_sheetkey}, -width=>2)->pack(-side=>'right');
	    }
	}
    }
    # actions

    my $level_action_frame=$level_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w',-fill=>"x",-pady=>5,-ipady=>3);

    # ok

    my $level_button_apply;

    my $level_button_ok=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},
						    -text=>"Ok",-command=>sub {
							$level_button_apply->invoke();
							$level_win->destroy();
						    })->pack(-side=>'left');

    $level_button_apply=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},
						    -text=>"Apply",-command=>sub {
							level_layout_all($annotation);
							view_refresh($main,$canvas,$annotation,$seqName,$seq,$plotstatus,$sheet);
						    })->pack(-side=>'left');

    my $level_button_cancel=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},
							-text=>"Cancel",-command=>sub {
							    $level_win->destroy();
							})->pack(-side=>'right');
}

sub select_levels {
    my $main=shift;
    my $canvas=shift;
    my $annotator_selection = shift;
    my $sheet = shift;

    my $level_win = $$main->Toplevel();
    $level_win->title("A GUI - select annotation types");
    $level_win->configure(-background=>'linen');
    
    my $level_main_frame = $level_win->Frame(-background=>$$sheet{default_win_background})->pack;
    
    # toggles

    my $level_check_frame = $level_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top', -anchor=>'w', -pady=>5);
    
   # sheet{"display_${key}_level"} = 1;

    foreach my $annotator (keys %annotatorName) {    
	my $this_annotator_frame = $level_check_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top', -anchor=>'w', -fill=>'x');
	if( !defined($$annotator_selection{$annotator}) ) {
	    $WARNING && print STDERR "WARNING: an undefined annotator selection field ($annotator) was encountered. Adding to dialog as unselected.\n";
	    $$annotator_selection{$annotator} = 0;
	}
	$this_annotator_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$annotator_selection{$annotator},-text=>"$annotatorName{$annotator}")->pack(-side=>"left");
    }
    # actions

    my $level_action_frame=$level_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w',-fill=>"x",-pady=>5,-ipady=>3);

    # ok

 #   my $level_button_apply;

    my $level_button_ok=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Ok",-command=>sub {
	# $level_button_apply->invoke();
	$level_win->destroy();
    })->pack(-side=>'left');

#    $level_button_apply=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Apply",-command=>sub {    })->pack(-side=>'left');

    my $level_button_cancel=$level_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub {
	$level_win->destroy();
    })->pack(-side=>'right');

    #obtain local grab? wait for destroy?

    $level_win->grab();
    $level_win->waitWindow();
}

sub view_center_on {
    my $main = shift;
    my $canvas = shift;
    my $annotation = shift;
    my $seq = shift;
    my $sheet = shift;

    my $goto_win = $$main->Toplevel();
    $goto_win->title("A GUI - Center on");
    $goto_win->configure(-background=>'linen');

    # predeclare
    my $goto_name_toggled = 0;
    my $goto_list_toggled = 0;
    my $goto_position_toggled = 0;
    my $goto_name_entry;
    my $goto_position_entry;
    
    my $goto_frame=$goto_win->Frame(-background=>$$sheet{default_win_background})->pack;

    # top

    my $goto_toggle_frame=$goto_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w');

    # position

    my $goto_toggle_position_frame=$goto_toggle_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top', -anchor=>'w');    
    my $goto_position_toggle=$goto_toggle_position_frame->Radiobutton(-command=> sub {
	if($goto_position_toggled) {
	    $goto_name_toggled = 0;
	    $goto_list_toggled = 0;
	    $goto_name_entry->configure(-state=>'disabled');
	    $goto_position_entry->configure(-state=>'normal');
	}
    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$goto_position_toggled, -value=>'1',-text=>"Goto position (1 - ".length($$seq).")")->pack(-side=>"left",-anchor=>"n");   
    

    my $goto_position;    
    $goto_position_entry=$goto_toggle_position_frame->Entry(-state=>'disabled', -background=>$$sheet{default_win_background},-textvariable=>\$goto_position)->pack(-anchor=>'n', -side=>"left");

    # name

    my $goto_name_toggle_frame = $goto_toggle_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w');
    my $goto_name_toggle=$goto_name_toggle_frame->Radiobutton(-command=> sub {
	if($goto_name_toggled) {
	    $goto_list_toggled = 0;
	    $goto_position_toggled = 0;
	    $goto_name_entry->configure(-state=>'normal');
	    $goto_position_entry->configure(-state=>'disabled');
	}
    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$goto_name_toggled, -value=>'1',-text=>"Goto annotation named ")->pack(-side=>"left",-anchor=>"n");
    my $goto_name;
    $goto_name_entry = $goto_name_toggle_frame->Entry(-state=>'disabled', -background=>$$sheet{default_win_background},-textvariable=>\$goto_name)->pack(-anchor=>'n', -side=>"left");
    
    # name list

    my $goto_list_toggle_frame = $goto_toggle_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>'top',-anchor=>'w');
    my $goto_list_toggle=$goto_list_toggle_frame->Radiobutton(-command=> sub {
	if($goto_list_toggled) {
	    $goto_name_toggled = 0;
	    $goto_position_toggled = 0;
	    $goto_name_entry->configure(-state=>'disabled');
	    $goto_position_entry->configure(-state=>'disabled');
	}
    },-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$goto_list_toggled, -value=>'1',-text=>"Goto ")->pack(-side=>"left",-anchor=>"n");
   
    my $goto_listbox=$goto_list_toggle_frame->Listbox(-relief => 'sunken', -height=>'20', -setgrid=>'true', -selectmode=> 'single',-background=>$$sheet{default_win_background})->pack(qw/-side left -expand yes -fill both/);

    for(my $i = 0; $i < $$annotation{nr}; $i++) {
	if(defined($$annotation{name}->[$i]) && ($$annotation{name}->[$i] ne "")) {
	    $goto_listbox->insert(0,$$annotation{name}->[$i]);	    
	}
    }

    # scrollbar for namelist

    my $goto_listbox_sb=$goto_list_toggle_frame->Scrollbar(-command => ['yview', $goto_listbox])->pack(-side=>'right',-fill=>'y');
  $goto_listbox->configure(-yscrollcommand => ['set', $goto_listbox_sb]);

    $goto_position_toggle->invoke();
    
    # action buttons

    my $goto_action_frame=$goto_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -anchor=>'w',-fill=>"x",-pady=>5,-ipady=>3);

    # ok
    my $goto_button_apply;

    my $goto_button_ok=$goto_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Ok",-command=>sub {
	$goto_button_apply->invoke();
	$goto_win->destroy();
})->pack(-side=>'left');

    # apply
    $goto_button_apply = $goto_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Apply",-command=>sub {
	my @current_scroll_pos = $$canvas->xview();
	my $current_scroll_width = $current_scroll_pos[1] - $current_scroll_pos[0];
	my $seq_pos = 0;
	my $nr=0;
	my $list_name = "";

	$DEBUG && print "DEBUG: Pressed apply...\n";
	if($goto_name_toggled == 1) {
	    $nr = annotation_what_nr_for_name($annotation, $goto_name);
	    if($nr >= 0) {
		$seq_pos = $$annotation{start}->[$nr] + ($$annotation{stop}->[$nr] - $$annotation{start}->[$nr] + 1)/2;
	    } else {
		$WARNING && print "WARNING: no such annotation named $goto_name found.\n";
		$seq_pos = ($current_scroll_pos[0] + $current_scroll_width/2) * length($$seq); # Goto $current_scrollwidth/2..
	    }
	    $DEBUG && print "DEBUG: Goto named $goto_name at $seq_pos...\n";
	} elsif ($goto_list_toggled == 1) {
	    $list_name = $goto_listbox->get('active');
	    $nr = annotation_what_nr_for_name($annotation,$list_name);
	    $seq_pos = $$annotation{start}->[$nr] + ($$annotation{stop}->[$nr] - $$annotation{start}->[$nr] + 1)/2;
	    $DEBUG && print "DEBUG: Goto listed $list_name at $seq_pos...\n";
	} elsif ($goto_position_toggled == 1) {
	    if ($goto_position < 0) {
		$goto_position = 0;
	    } elsif ( $goto_position > length($$seq) ) {
		$goto_position = length($$seq);
	    }
	    $seq_pos = $goto_position;
	    $DEBUG && print "DEBUG: Goto position $seq_pos...\n";
	}
	my $new_scroll_pos = $seq_pos / length($$seq) - $current_scroll_width/2 ;
	$DEBUG && print "DEBUG: Goto position $new_scroll_pos since scroll width is $current_scroll_width, seq pos $seq_pos and sequence length is ", length($$seq),"...\n";
	$$canvas->xview(moveto=>$new_scroll_pos);
    })->pack(-side=>'left');

    # cancel

    my $goto_button_cancel=$goto_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub {$goto_win->destroy;})->pack(-side=>"right");
}

# Annotation menu

sub annotation_merge {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift; #ref
    my $seqName=shift; # ref
    my $seq=shift;# ref
    my $sheet=shift; #ref

    my %merge_status;
    $merge_status{shortest_olap}=50; # Default value

    my $merge_win=$$main->Toplevel;
    $merge_win->title("A GUI - Merge annotations");
    
    my $mw_main_frame=$merge_win->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'both');
    
    # Naming: Tryposomatid, Yeast or Clonebased
    my $mw_naming_frame=$mw_main_frame->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'both',-side=>'top');
    my $mw_naming_label=$mw_naming_frame->Label(-text=>'Automatic naming')->pack(-anchor=>'w',-side=>'top');
    my $mw_naming_toggle_trypanosomatid=$mw_naming_frame->Checkbutton(-state=>'disabled',-text=>'Apply trypanosomatid homology naming',-variable=>\$merge_status{tryp_name})->pack(-side=>'top',-anchor=>'w');
    my $mw_naming_toggle_yeast=$mw_naming_frame->Checkbutton(-state=>'disabled',-text=>'Apply yeast homology naming',-variable=>\$merge_status{yeast_name})->pack(-side=>'top',-anchor=>'w');
    my $mw_naming_toggle_clone=$mw_naming_frame->Checkbutton(-text=>'Apply clone based naming',-variable=>\$merge_status{clone_name})->pack(-side=>'top',-anchor=>'w');

    my $mw_annot_type_frame=$mw_main_frame->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'both',-side=>'top');
    my $mw_annot_label=$mw_annot_type_frame->Label(-text=>'Merge annotations')->pack(-anchor=>'w',-side=>'top');
    my $mw_toggle_orf=$mw_annot_type_frame->Checkbutton(-text=>'ORF',-variable=>\$merge_status{ORF})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_glimmer=$mw_annot_type_frame->Checkbutton(-text=>'Glimmer2',-variable=>\$merge_status{glimmer2})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_testcode=$mw_annot_type_frame->Checkbutton(-text=>'Testcode',-variable=>\$merge_status{testcode})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_estorf=$mw_annot_type_frame->Checkbutton(-text=>'ESTORF',-variable=>\$merge_status{ESTORF})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_manual=$mw_annot_type_frame->Checkbutton(-text=>'Manual',-variable=>\$merge_status{manual})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_blast=$mw_annot_type_frame->Checkbutton(-text=>'Blast',-variable=>\$merge_status{blast})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_genbankcds=$mw_annot_type_frame->Checkbutton(-text=>'Genbank CDS',-variable=>\$merge_status{genbankcds})->pack(-side=>'top',-anchor=>'w');

    my $mw_method_frame=$mw_main_frame->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'both',-side=>'top');
    my $mw_method_label=$mw_method_frame->Label(-text=>'Merge method')->pack(-anchor=>'w',-side=>'top');
    my $mw_toggle_or=$mw_method_frame->Checkbutton(-text=>'or (any of the selected above)',-variable=>\$merge_status{or})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_and=$mw_method_frame->Checkbutton(-text=>'and exactly (exact matches of the above)',-variable=>\$merge_status{and})->pack(-side=>'top',-anchor=>'w');
    my $mw_toggle_comp_frame=$mw_method_frame->Frame()->pack(-expand=>'yes',-side=>'top',-fill=>'x');
    my $mw_toggle_comp=$mw_toggle_comp_frame->Checkbutton(-text=>'and overlap pair (any of the above) with overlap longer than ',-variable=>\$merge_status{comp})->pack(-side=>'left');
    my $mw_toggle_comp_entry=$mw_toggle_comp_frame->Entry(-textvariable=>\$merge_status{shortest_olap})->pack(-side=>'left');
    
    my $mw_param_frame=$mw_main_frame->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'both',-side=>'top');
    my $mw_param_label=$mw_param_frame->Label(-text=>"Options")->pack(-anchor=>'w',-side=>'top');
    # If start and stop are equal, but orientation differs someone fucked up bigtime, so that is _not_ to be considered an option.
    my $mw_toggle_frames=$mw_param_frame->Checkbutton(-text=>"Consider frame orientation for overlaps",-variable=>\$merge_status{frame})->pack(-anchor=>'w',-side=>'top');
    my $mw_toggle_framex=$mw_param_frame->Checkbutton(-text=>"Consider exact frame for overlaps",-variable=>\$merge_status{framex})->pack(-anchor=>'w',-side=>'top');

    my $mw_action_frame=$mw_main_frame->Frame()->pack(-padx=>6, -pady=>6,-expand=>'yes',-fill=>'x',-side=>'bottom');
    my $mw_apply;
    my $mw_ok=$mw_action_frame->Button(-text=>'Ok',-command=>sub { 
	$mw_apply->invoke();
	$merge_win->destroy();
    })->pack(-side=>'left');
    $mw_apply=$mw_action_frame->Button(-text=>'Apply',-command=>sub {
	merge($main,$canvas,$annotation,\%merge_status,$sheet,$seqName);
	$$sheet{status}="dirty"; # assume something changed
	redraw_annotations($canvas,$annotation,$sheet,$seq);
    })->pack(-side=>'left');
    my $mw_cancel=$mw_action_frame->Button(-text=>'Cancel',-command=>sub { 
	$merge_win->destroy();
    })->pack(-side=>'right');  
}

sub annotation_name {
    my $annotation=shift;
    my $annotation_nr=shift;
    my $seqName=shift;
    my $name_style=shift;

    # Get clone name, or use seqName if not found...

    my $name_base="";
    my @used_names=();

    if($name_style eq "clone") {
	if(open(CLONENAME,"<$cloneNameDb") || !(messageBox(-icon=>'error',-type=>'OK',-message=>"Could not open clone name database file $cloneNameDb"))) {
	    while(<CLONENAME>) {
		# read lines for "official" clone-name and check for a valid synonym.
		
	    }
	    close(CLONENAME);
	} else { 
	    $DEBUG && print "Could not open clone name database file $cloneNameDb. Reverting to seqName.\n";
	    $name_base=$$seqName;
	}

	# Check for first unused clonename... (oh, and lock the namedb first, by the way?)
	# Add a preferences tab for changing the location of these files?
	if (!open(NAMEDB,"+>$nameDb")) {
	    $DEBUG && messageBox(-icon=>'error',-type=>'OK',-message=>"Could not open CDS name database file $nameDb (for read and write)");
	    return;
	}
	
	while(<NAMEDB>) {
	    # Construct a used name (number) list from the namedb
	}
	close(NAMEDB);
	
    } else {
	# Other kind of naming.. 
	
    }
    
}

sub annotation_remove {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;

    my %rm_field;
    
    # Introduce accelerator key // pop-up menu for "Remove currently selected annotation.."?? 
    my $remove_win=$$main->Toplevel;
    $remove_win->title("A GUI - Remove annotation");
    # $remove_win->geometry('+0+0');
    $remove_win->configure(-background=>$$sheet{default_win_background},-width=>300);

    $rm_field{id}="Id: ";
    $rm_field{name}="Name: ";
    $rm_field{type}="Type: ";
    $rm_field{start}="Start: ";
    $rm_field{stop}="Stop: ";
    $rm_field{frame}="Frame: ";
    
    my $rm_frame=$remove_win->Frame(-background=>$$sheet{default_win_background})->pack(-padx=>'6', -pady=>'6',-expand=>"yes", -fill=>"both");
    my $rm_title=$rm_frame->Label(-background=>$$sheet{default_win_background},-text=>"Selected annotation")->pack(-expand=>'yes',-fill=>'x',-anchor=>'w',-side=>'top');

    #  my $rm_id_frame=$rm_frame->Frame(-background=>$$sheet{default_win_background})->pack(-expand=>"yes", -fill=>'x',-side=>'top',-anchor=>'w');  
    my $rm_id_label=$rm_frame->Label(-textvariable=>\$rm_field{id},-background=>$$sheet{default_win_background})->pack(-expand=>"yes", -fill=>'x', -anchor=>'w',-side=>'top',-anchor=>'w');
    my $rm_type_label=$rm_frame->Label(-textvariable=>\$rm_field{type},-background=>$$sheet{default_win_background})->pack(-expand=>"yes", -anchor=>'w',-fill=>'x',-side=>'top',-anchor=>'w');
    my $rm_locus_frame=$rm_frame->Frame(-background=>$$sheet{default_win_background})->pack(-expand=>"yes", -fill=>"x",-side=>"top",-anchor=>'w');
    my $rm_start_label=$rm_locus_frame->Label(-textvariable=>\$rm_field{start},-background=>$$sheet{default_win_background})->pack(-side=>'left');
    my $rm_stop_label=$rm_locus_frame->Label(-textvariable=>\$rm_field{stop},-background=>$$sheet{default_win_background})->pack(-side=>'left');
    my $rm_frame_label=$rm_locus_frame->Label(-textvariable=>\$rm_field{frame},-background=>$$sheet{default_win_background})->pack(-side=>'left');
    
    my $rm_action_frame = $rm_frame->Frame(-background=>$$sheet{default_win_background})->pack(-expand=>"yes", -fill=>"x",-side=>"top");
    my $rm_remove_button= $rm_action_frame->Button(-text=>"Remove",-background=>$$sheet{default_win_background},-activebackground=>$$sheet{dialog_button_active_background},-command=> sub {
	my ($annotation_id)=$rm_field{id}=~m/Id:\s{1}([\w\d]+)/;
	if ($annotation_id) {
	    my @annotation_nr;
	    $annotation_nr[0]=annotation_what_nr($annotation,$annotation_id);
	    $DEBUG && print "DEBUG: call remove for $annotation_nr[0]\n";
	    remove_annotations($annotation,\@annotation_nr,$sheet);

	    redraw_annotations($canvas,$annotation,$sheet,\$seq); # using seq as global..
	}
    })->pack(-side=>'left');
    my $rm_update_button=$rm_action_frame->Button(-text=>"Refresh",-background=>$$sheet{default_win_background},-activebackground=>$$sheet{dialog_button_active_background},-command=> sub { fill_in_remove_fields($canvas,\%rm_field,$annotation); })->pack(-side=>'left');
    my $rm_cancel_button=$rm_action_frame->Button(-text=>"Cancel",-background=>$$sheet{default_win_background},-activebackground=>$$sheet{dialog_button_active_background},-command=> sub { $remove_win->destroy; })->pack(-side=>'right');

    fill_in_remove_fields($canvas,\%rm_field,$annotation,$sheet);
    check_selection($remove_win,\%rm_field,$canvas,$annotation,$sheet);
}

sub fill_in_remove_fields {
    my $canvas=shift;
    my $rm_field=shift;
    my $annotation=shift;
    my $sheet = shift;

    my ($annotation_id, $annotation_type) = annotation_id_type_selected($annotation, $sheet);  
    
    if($annotation_id eq "error") {
	$$rm_field{id}="Id: ";
	$$rm_field{name}="Name: ";
	$$rm_field{type}="Type: ";
	$$rm_field{start}="Start: ";
	$$rm_field{stop}="Stop: ";
	$$rm_field{frame}="Frame: ";

	return 0;
    }

    my $annotation_nr = annotation_what_nr($annotation,$annotation_id);
    
    $$rm_field{id}="Id: $annotation_id";
    if($$annotation{name}->[$annotation_nr]) {
	$$rm_field{name}="Name: $$annotation{name}->[$annotation_nr]";
    }
    $$rm_field{type}="Type: $annotation_type";
    $$rm_field{start}="Start: $$annotation{start}->[$annotation_nr]";
    $$rm_field{stop}="Stop: $$annotation{stop}->[$annotation_nr]";
    if($annotation_type eq "EST") {
	$$rm_field{frame}="Frame: any";
    } else {
	$$rm_field{frame}="Frame: $$annotation{frame}->[$annotation_nr]";
    }
    return 0;
}

sub check_selection {
    # Idea: check selection ever-so-often.
    # It would be better/cheaper if we could check for the existance of a remove window in 
    # left-mouse-click... or send an event indicating selection update to catch..
    my $remove_win=shift;
    my $rm_field=shift;
    my $canvas=shift;
    my $annotation = shift;
    my $sheet = shift;

    fill_in_remove_fields($canvas,$rm_field,$annotation,$sheet);
    
    if (Exists($remove_win)) {
	$remove_win->after(400,sub { check_selection($remove_win,$rm_field,$canvas, $annotation,$sheet) });
	return 1;
    }   

    return 0;
}

sub annotation_removeAll {
    my $main=shift;
    my $canvas=shift;
    my $annotation=shift;
    my $seq=shift; # ref
    my $sheet=shift;

    my %removestate; # {glimmerall, orfestc1, orfestc2, orfestshadowed, ...}

    my $remove_all_win=$$main->Toplevel;
    $remove_all_win->title("A GUI - Remove all");
    # $remove_all_win->geometry('+0+25');
    $remove_all_win->configure(-width=>300);

    my $ra_frame=$remove_all_win->Frame()->pack(-expand=>"yes", -fill=>"both");

    my $ra_notebook=$ra_frame->NoteBook(-ipadx=>6, -ipady=>6);
    my %ra_page;
    foreach my $annotator (keys %annotatorName) {
	$ra_page{$annotator} = $ra_notebook->add($annotator, -label => $annotatorBrev{$annotator}[0], -underline=>$annotatorBrev{$annotator}[1]);

	# ubiquitous all & size removal -- exceptions allowed, e g ST-size is ridiculous
	$ra_page{$annotator}->Checkbutton(-variable=>\$removestate{$annotator."all"},-text=>"All",-underline=>0)->pack(-side=>"top",-anchor=>"w");

	$nosizeremove = 0;
	foreach my $pref (@{ $annotatorPrefs{$annotator} }) {
	    $DEVEL && $DEBUG && print "DEBUG: Found pref $pref when checking prefs for $annotator.\n";
	    if ($pref eq "nosizeremove") {
		$nosizeremove = 1;
	    }
	}
	unless ($nosizeremove) {
	    my $ra_size_frame=$ra_page{$annotator}->Frame()->pack(-side=>'top', -anchor=>'w');
	    $ra_size_frame->Label(-text=>$annotatorName{$annotator}." shorter than")->pack(-side=>'left');
	    $ra_size_frame->Entry(-textvariable=>\$removestate{$annotator."size"})->pack(-side=>'left');
	}
    }

    # type specific removall options
    $ra_page{glimmer2}->Label(-text=>"Glimmer comment matching any")->pack(-side=>'top', -anchor=>'w');
    my $ra_glimmer_toggle_shadowed=$ra_page{glimmer2}->Checkbutton(-variable=>\$removestate{glimmershadowed},-text=>"ShadowedBy")->pack(-side=>"top",-anchor=>"w");
    my $ra_glimmer_toggle_shorterthan=$ra_page{glimmer2}->Checkbutton(-variable=>\$removestate{glimmershorterthan},-text=>"ShorterThan")->pack(-side=>"top",-anchor=>"w");
    my $ra_glimmer_toggle_olap=$ra_page{glimmer2}->Checkbutton(-variable=>\$removestate{glimmerolapwith},-text=>"OlapWith")->pack(-side=>"top",-anchor=>"w");

    my $ra_testcode_p_frame=$ra_page{testcode}->Frame()->pack(-side=>'top', -anchor=>'w');
    $ra_testcode_p_frame->Label(-text=>"Testcode ORF with average lower than")->pack(-side=>'left');
    my $ra_testcode_p_entry=$ra_testcode_p_frame->Entry(-textvariable=>\$removestate{testcodep})->pack(-side=>'left');

    my $ra_estorf_toggle_c1=$ra_page{ESTORF}->Checkbutton(-variable=>\$removestate{estorfc1},-text=>"Class 1")->pack(-side=>"top",-anchor=>"w");
    my $ra_estorf_toggle_c2=$ra_page{ESTORF}->Checkbutton(-variable=>\$removestate{estorfc2},-text=>"Class 2")->pack(-side=>"top",-anchor=>"w");
    my $ra_estorf_toggle_shadowed=$ra_page{ESTORF}->Checkbutton(-variable=>\$removestate{estorfredundant},-text=>"Invisible")->pack(-side=>"top",-anchor=>"w");

    my $ra_est_p_frame=$ra_page{EST}->Frame()->pack(-side=>'top', -anchor=>'w');
    $ra_est_p_frame->Label(-text=>"EST hit P higher than")->pack(-side=>'left');  
    my $ra_est_p_entry=$ra_est_p_frame->Entry(-textvariable=>\$removestate{estp})->pack(-side=>'left');

    my $ra_blast_p_frame=$ra_page{blast}->Frame()->pack(-side=>'top', -anchor=>'w');
    $ra_blast_p_frame->Label(-text=>"Blast hit P higher than")->pack(-side=>'left');  
    my $ra_blast_p_entry=$ra_blast_p_frame->Entry(-textvariable=>\$removestate{blastp})->pack(-side=>'left');

    my $ra_orf_toggle_strand_heuristic=$ra_page{ORF}->Checkbutton(-variable=>\$removestate{orfheur},-text=>"According to strandedness heuristic")->pack(-side=>"top",-anchor=>"w");

    $ra_notebook->pack(-expand=>'yes',
		       -fill=>'both', 
		       -padx=>5, 
		       -pady=>5, 
		       -side=>'top');

    my $ra_action_frame=$ra_frame->Frame()->pack(-expand=>"yes", -fill=>"x",-side=>"top",-anchor=>"sw");
    my $ra_button_apply;
    my $ra_button_ok=$ra_action_frame->Button(-text=>"Ok",-command=>sub {
	$ra_button_apply->invoke();
	$remove_all_win->destroy;
    })->pack(-side=>"left");
    $ra_button_apply = $ra_action_frame->Button(-text=>"Apply",-command=>sub {
	# Known bug: sets dirty even if nothing removed.
	# Solve eg by accepting return status from remove_all.
	remove_all($annotation,\%removestate,$seq,$sheet);
	$$sheet{status}="dirty"; # Sheet is now dirty.
	redraw_annotations($canvas,$annotation,$sheet,$seq);
    })->pack(-side=>"left");
    my $ra_button_cancel=$ra_action_frame->Button(-text=>"Cancel",-command=>sub {$remove_all_win->destroy;})->pack(-side=>"right");
}

sub remove_all  {
    my $annotation=shift;
    my $removestate=shift;
    my $seq=shift; # ref
    my $sheet = shift;

    my @remove_us=();

    # Cover for the ubiquitous remove-all and remove-size
    my $removals = 0;
    foreach my $annotator (keys %annotatorName) {
	if($$removestate{$annotator."all"}) {
	    # Find nr of annotation type annotator
	    for(my $i=0;$i<$$annotation{nr};$i++)  {
		if($$annotation{type}->[$i] eq $annotator) {
		    push @remove_us,$i;
		    $removals = 1;
		}
	    }
	    $DEBUG && print "DEBUG: $annotator-all removes queued for remove...\n";

	    $$removestate{$annotator."all"}=0;
	}	
    }
    # Run this once or several times? Where is Knuth when you need him...
    if($removals) {
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
    }
   
    $removals = 0;
    foreach my $annotator (keys %annotatorName) {
	if($$removestate{$annotator."size"}) {
	    # Find nr's of this type of annotations shorter than..
	    for(my $i=0;$i<$$annotation{nr};$i++)  {
		if($$annotation{type}->[$i] eq $annotator) {
		    if($$annotation{stop}->[$i]-$$annotation{start}->[$i]<$$removestate{$annotator."size"}) {
			push @remove_us,$i;
			$removals = 1;
		    }
		}
	    }
	    $$removestate{$annotator."size"}="";
	    $DEBUG && print "DEBUG: $annotator-size removes queued for remove...\n";
	}
    }
    # Run this once or several times? Where is Knuth when you need him...
    if($removals) {
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
    }

    if($$removestate{glimmershadowed}) {
	# Find nr's of glimmer annotations with comment containing a ShadowedBy
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "glimmer2") {
		if($$annotation{comment}->[$i] =~ m/ShadowedBy/) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Glimmershadow remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{glimmershadowed}="";
    }

	if($$removestate{glimmershorterthan}) {
	# Find nr's of glimmer annotations with comment containing a ShorterThan
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "glimmer2") {
		if($$annotation{comment}->[$i] =~ m/ShorterThan/) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Glimmershorterthan remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);   
	@remove_us=();
	$$removestate{glimmershorterthan}="";
    }

	if($$removestate{glimmerolapwith}) {
	# Find nr's of glimmer annotations with comment containing OlapWith
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "glimmer2") {
		if($$annotation{comment}->[$i] =~ m/OlapWith/) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Glimmerolapwith remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{glimmerolapwith}="";
    }

    if($$removestate{estorfredundant}) {
	# Check for and remove identical ORFs, as well as invisible ORFs.
	# Should include an option to do this, to shed some more light on what the comparison routine actually 
	# does.. Perhaps its nice to have on the remove menu?
	
	my @orfnrs;
	my @orfStart;
	my @orfEnd;
	my @orfClass;
	my $orfClass; # Might be a little confusing...

	for(my $j=0;$j<$$annotation{nr};$j++)  {
	    if($$annotation{type}->[$j] eq "ESTORF") {
		# State variables
		my $shadowDetected=0;
		my $shadowResolved=0;

		$orfStart=$$annotation{start}->[$j];
		$orfStop=$$annotation{stop}->[$j];
		($orfClass)=$$annotation{comment}->[$j]=~m/Hit class (\d+)/;
		$DEBUG && print "DEBUG: looking at a class $orfClass hit $orfStart-$orfStop ... ";
#	$orfPhase=$$annotation{frame}->[$j];

		# Note that we have only one seq, so the original compare-code is greatly simplified.
		if(@orfnrs) {
		    for(my $i=0;$i<@orfnrs;$i++) {
			# Strategy: keep longest ORF if there are several shadowed ones sharing the same stop.
			if( $orfStart == $orfBegin[$i] ) {
			    $shadowDetected=1;
# C'mon, things should never have the same start and different stop. This is an identical, or an error.
#	      if(abs($orfEnd[$i]-$orfBegin[$i])<abs($orfStop-$orfStart)) {
#		# We need only adjust the unchecked end...
#		$orfEnd[$i]=$orfStop;
# 	      }
			}	  
			# For identicals, this is not an elsif... (Identicals simply cause no action and shadowFlag set.)
			if( $orfStop == $orfEnd[$i] ) {
			    if($shadowDetected) {		
				# This is an identical (or an error). Check class, and keep only the best one.   
				if($orfClass<$orfClass[$i]) {
				    $DEBUG && print "... better than old one ($orfnrs[$i]) with same coords.\n"; # DEBUG
				    # Delete the old identical with lower class..
				    push @remove_us,$orfnrs[$i];
				    # And keep this!
				    $orfClass[$i]=$orfClass;
				    $orfnrs[$i]=$j;
				    # We don't want this added as a "novel ORF", but neither do we wish to
				    # have it removed. Enable shadowResolved to bypass "novelAddition".
				    $shadowDetected=0;
				    $shadowResolved=1;
				}
			    } else {
				#Ok, worse class, identical length so dump this.
				$shadowDetected=1;	
			    }

			    # Abs needed for the reverse reading frames.
			    if(abs($orfEnd[$i]-$orfBegin[$i])<abs($orfStop-$orfStart)) {
				$DEBUG && print "... longer than old one ($orfnrs[$i]) with same stop.\n"; # DEBUG
				# Delete the old shorter one..
				push @remove_us,$orfnrs[$i];
				# ..and update the remaining ORF entry.
				$orfBegin[$i]=$orfStart;
				$orfClass[$i]=$orfClass;
				$orfnrs[$i]=$j;
				$shadowDetected=0;
				$shadowResolved=1;
			    } 
			}          
		    }			# If $orfnrs == 0 we end up here immediately, without shadowDetected...
		    
		    if(!$shadowDetected && !$shadowResolved) {
			# $DEBUG && print "DEBUG:$subjectName touches ORF $orfStart - $orfStop $orfPhase\n";
			# Ok, this is a shadow-free, novel ORF	  
			$orfBegin[@orfnrs]=$orfStart;
			$orfEnd[@orfnrs]=$orfStop;
			$orfClass[@orfnrs]=$orfClass;
			$orfnrs[@orfnrs]=$j;
			$DEBUG && print "... it is novel, so it is added.\n"; # DEBUG
		    } elsif($shadowDetected){	    
			# This ORF is already annontated, or is shorter than another ORF 
			# with the same stop codon.	    
			push @remove_us,$j;
			$DEBUG && print "... dumped due to shadow.\n"; # DEBUG
		    }
		    $shadowResolved=0;
		} else {  
		    # Create entry for the first ORF.. 
		    $orfBegin[@orfnrs]=$orfStart;
		    $orfEnd[@orfnrs]=$orfStop;
		    $orfClass[@orfnrs]=$orfClass;	  
		    $orfnrs[@orfnrs]=$j;
		    $DEBUG && print "... its the first ORF to be added.\n"; # DEBUG
		}
	    }
	}
	$DEBUG && print "DEBUG: Estorfinvis remove calling remove...\n";
	# Run this once or several times? Where is Knuth when you need him...
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{estorfredundant}=0;
    }

    if($$removestate{estorfc1}) {
	# Find nr's of ESTORF annotations with comment containing Class 1
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "ESTORF") {
		if($$annotation{comment}->[$i] =~ m/Hit class 1/) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Estorfc1 remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);   
	@remove_us=();
	$$removestate{estorfc1}="";
    }

    if($$removestate{estorfc2}) {
	# Find nr's of ESTORF annotations with comment containing Class 2
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "ESTORF") {
		if($$annotation{comment}->[$i] =~ m/Hit class 2/) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Estorfc2 remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{estorfc2}="";
    }

    if($$removestate{estorfshadowed}) {
	# Find nr's of glimmer annotations shorter than..
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "ESTORF") {
		#if ($shadowed) {
		#  push @remove_us,$i;	
		#}
	    }    
	}
	$DEBUG && print "DEBUG: Estorfshadowed remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{estorfshadowed}="";
    }  

    if($$removestate{estp}) {
	# Find nr's of EST annotations with P higher than..
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "EST") {
		($pval) = $$annotation{comment}->[$i] =~ /,P=([\d+-eE.]+),/;
		if( $pval > $$removestate{estp}) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Estsize remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{estp}="";
    }

    if($$removestate{testcodep}) {
	# Find nr's of TESTCODE annotations with "P" lower than..
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "testcode") {
		$$annotation{comment}->[$i] =~ /ORF Testcode score ([\de+-.]+)/;
		$pval = $1;

		if( $pval < $$removestate{testcodep}) {
		    push @remove_us,$i;
		}
	    }
	}
	$DEBUG && print "DEBUG: Testcodesize remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);
	@remove_us=();
	$$removestate{testcodep}="";
    }

    if($$removestate{blastp}) {
	# Find nr's of blast annotations with P higher than..
	for(my $i=0;$i<$$annotation{nr};$i++)  {
	    if($$annotation{type}->[$i] eq "blast") {
		($pval) = $$annotation{comment}->[$i] =~ /,P=([\d+-eE.]+),/;
		if( $pval > $$removestate{blastp}) {
		    push @remove_us,$i;	
		}
	    }    
	}  
	$DEBUG && print "DEBUG: Blastsize remove calling remove...\n";
	remove_annotations($annotation,\@remove_us,$sheet);    
	@remove_us=();
	$$removestate{blastp}="";
    }

    if($$removestate{orfheur}) {
	heuristic_likely_orf($annotation,$seq,$sheet); 

	@remove_us=();
	$$removestate{orfheur}=0;
    }
}

sub remove_annotations {
    my $annotation=shift;
    my $remove_us=shift;
    my $sheet = shift;

    my $newnr;  
    my $oldnr=$$annotation{nr};

    my $removedElements;
    
    # Hey, this assumes remove_us sorted. Better do that; you never know when you'll forget this..
    @$remove_us = sort {$a <=> $b} @$remove_us;

    # Do the actual moving... Dont know what Knuth would do, but this is O($field*$oldnr). 
    # Could be a lot worse, and not so extremely much better (with this datatype, ofcourse). 
    $DEBUG && print "DEBUG: attempt remove of nr: @$remove_us,\n";

    $newnr=0;
    $removedElements=0;
    for(my $i=0;$i<$oldnr;$i++) {
	if(($removedElements<@$remove_us) && ($i==$$remove_us[$removedElements])) {
	    # If the first relation is false, we have now removed all elements. 
	    # The only thing left is to do is to fill up the new vector. Besides, 
	    # the second relation would be a "use if uninit" in this case...

	    # Dont insert, dont increase newnr.
	    if ($$annotation{id}->[$i] eq $$sheet{selected}) {
		$DEBUG && print "DEBUG: Annotation $$annotation{id}->[$i] is the selected annotation ($$sheet{selected}). Deselecting on removal..\n";
		$$sheet{selected} = '';
	    }
	    $removedElements++;
	} else {
	    foreach $field (keys %$annotation) {
		next if $field eq "nr";
		next if $field eq "unr";
		# $DEVEL && print "DEVEL temp: annotation nr 1 now has $field = $$annotation{$field}->[1].\n";
		if($field eq "id") {
		    # Adjust Id-tags on the remaining annotations
		    my ($old_annotation_id_wo_nr)=$$annotation{id}->[$i]=~/(\w+)_\d+/;
		    $DEBUG && print "DEBUG: Annotation $$annotation{id}->[$i] remapped to id ",$old_annotation_id_wo_nr,"_",$newnr,"\n";
		    if ($$annotation{id}->[$i] eq $$sheet{selected}) {
			$DEBUG && print "DEBUG: Annotation $$annotation{id}->[$i] is the selected annotation ($$sheet{selected}). Changing selection...\n";
			$$annotation{id}->[$newnr]=$old_annotation_id_wo_nr . "_" . $newnr; # re-id sequences
			$$sheet{selected} = $$annotation{id}->[$newnr];
		    } else {
			$$annotation{id}->[$newnr]=$old_annotation_id_wo_nr . "_" . $newnr; # re-id sequences
		    }
		} else {
		    $$annotation{$field}->[$newnr]=$$annotation{$field}->[$i];
		    # Design question: should we keep some unchanging tag for cross-reference here, or
		    # is that simply not needed?
		    # Problems occur when new annotaions are added once again. Then the annotation{nr} gives the 
		    # current _actual_number off annotations. The problem could be avoided by setting annotation{nr} to 
		    # point at the highest used annotation nr. If we are to inroduce annotation names, the point in
		    # keeping constant ids dissapears - ids can once again be program internal.
		    
		    # Ok, keeping a uid to be able to track'em. How I wish this was object oriented instead..
		}
	    }
	    $newnr++;
	}
    }

    # Number of annotations + 1
    $$annotation{nr}=$newnr;
    
    foreach $field (keys %$annotation) {
	next if $field eq "nr";
	next if $field eq "unr";
	# $DEBUG && print "DEBUG: deleting the remaining old annotations for field $field nr $newnr up to ".($oldnr-1).".\n";
	foreach my $t ($newnr .. ($oldnr -1)) {
	    # $DEBUG && print "Counting $t..\n"; 
	   $$annotation{$field}->[$t] = undef;
	}
	# $DEVEL && print "DEVEL temp: annotation nr 1 now has $field = $$annotation{$field}->[1].\n";
    }
    update_annotation_nr_cache($annotation);
    # oops! the canvas will not be up to date at this stage! cannot trust withtag selected.
    my $selected_nr = annotation_what_selected_nr($annotation, $sheet);
    
    $$sheet{selected} = ($selected_nr > 0) ? $$annotation{id}->[$selected_nr] : '';
}

sub update_annotation_nr_cache {
    my $annotation = shift;
    # reset id cache
    %annotation_nr_cache = ();
    for(my $i=0;$i<$$annotation{nr};$i++)  {
	$annotation_nr_cache{$$annotation{id}->[$i]} = $i;
    }
}

sub annotation_what_selected_nr {
    my $annotation = shift;
    my $sheet = shift;
    
    return (($$sheet{selected} ne '') ? $annotation_nr_cache{$$sheet{selected}}:-1);
}

sub redraw_annotations {
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $seq=shift; # ref

    $DEBUG && ((ref $seq) || die "DEBUG FAIL: redraw_annotations got seq $seq instead of reference..\n");
    $$canvas->delete(-tags=>"annotation");

    # $DEBUG && print "DEBUG: Redraw $$annotation{nr} annotations\n";
    return 0 if (! defined($$annotation{nr}));
    
    for(my $i=0;$i<$$annotation{nr};$i++) {
	$DEVEL && $DEBUG && print "DEVEL DEBUG: redrawing annotation $i out of ".($$annotation{nr}-1).".\n";
        # eq "EST") {
	if($annotatorDirected{$$annotation{type}->[$i]} == 2) {
	    if ($$annotation{frame}->[$i] == $AnnotationBoxFrame) { # $frame eq Any
		draw_annotation_box($canvas,$annotation,$i,$sheet);
	    } else {
		draw_annotation_arrow($canvas,$annotation,$i,$sheet,$seq);
	    }
	} elsif ($annotatorDirected{$$annotation{type}->[$i]} == 1) {
	    draw_annotation_arrow($canvas,$annotation,$i,$sheet,$seq);
	} elsif ($$annotation{type}->[$i] eq "ST") { # directed == 0 or undef
	    draw_annotation_st($canvas,$annotation,$i,$sheet);
	} else {
	    $WARNING && print "WARNING: unknown annotation type ".$$annotation{type}->[$i]." for annotation nr $i.\n";
	}
    }
    $main->update(); # havent got a main ref!
}

sub redraw_selected_annotation {
    my $canvas=shift;
    my $annotation=shift;
    my $sheet=shift;
    my $seq=shift; # ref 
    my $last_selected_annotation_id = shift;   
    my $last_selected_annotation_nr;

    # Find the currently selected annotation before deleting
    my $selected_annotation_nr = annotation_what_selected_nr($annotation,$sheet);
    my $selected_annotation_id = $$annotation{id}->[$selected_annotation_nr];
   
    if(!defined($last_selected_annotation_id) or $last_selected_annotation_id eq "") {
	# no previous selection?
	$last_selected_annotation_id = $selected_annotation_id;
    }
    $last_selected_annotation_nr = annotation_what_nr($annotation, $last_selected_annotation_id);
    
    $$canvas->delete(-tags=>"$selected_annotation_id");
    $$canvas->delete(-tags=>"$last_selected_annotation_id");

    foreach my $i ($last_selected_annotation_nr, $selected_annotation_nr) {
	if($annotatorDirected{$$annotation{type}->[$i]} == 2) {
	    if ($$annotation{frame}->[$i] == $AnnotationBoxFrame) { # $frame eq Any
		draw_annotation_box($canvas,$annotation,$i,$sheet);
	    } else {
		draw_annotation_arrow($canvas,$annotation,$i,$sheet,$seq);
	    }
	} elsif ($annotatorDirected{$$annotation{type}->[$i]} == 1) {
	    draw_annotation_arrow($canvas,$annotation,$i,$sheet,$seq);
	} elsif ($$annotation{type}->[$i] eq "ST") { # directed == 0 or undef
	    draw_annotation_st($canvas,$annotation,$i,$sheet);
	}
    }
}

sub view_graph {
    my $main=shift;
    my $canvas=shift;
    my $seq=shift;
    
    my $plotstate=shift;
    my $sheet=shift;

    my $replot=1;

    my @xAxis;
    my @gcFraction;
    my @atFraction;
    my @agFraction;
    my @ctFraction;
    my @HWin;

#  if(!$$plotstate{windowSize}) {
#    $$plotstate{windowSize}=100; # Default window size
#  }

    # graph dialog..
    my $graph_win=$$main->Toplevel;
    $graph_win->title("A GUI - Graph settings");
    # $graph_win->geometry('+0+25');
    $graph_win->configure(-background=>'linen');
    
    my $graph_frame=$graph_win->Frame(-background=>$$sheet{default_win_background})->pack;

    my $graph_toggle_main_frame=$graph_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top", -fill=>"x");
    my $graph_toggle_left_frame=$graph_toggle_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"left");
    my $graph_toggle_right_frame=$graph_toggle_main_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"left");
    
    my $graph_GC_toggle=$graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{GC},-text=>"GC")->pack(-side=>"top",-anchor=>"w");
    my $graph_AT_toggle=$graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{AT},-text=>"AT")->pack(-side=>"top",-anchor=>"w");
    
    my $graph_GCs_toggle=$graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{GCs},-text=>"GC Skew")->pack(-side=>"top",-anchor=>"w");
    my $graph_ATs_toggle=$graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{ATs},-text=>"AT Skew")->pack(-side=>"top",-anchor=>"w");

#  my $graph_P_frame=$graph_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top");
    my $graph_AG_toggle=$graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{AG},-text=>"Pu (AG)")->pack(-anchor=>'w',-side=>"top");
    my $graph_CT_toggle=$graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{CT},-text=>"Py (CT)")->pack(-anchor=>"w",-side=>"top");
    
#  my $graph_info_frame=$graph_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top");
    my $graph_H_toggle=$graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{H},-text=>"H/2")->pack(-anchor=>"w",-side=>"top");
    my $graph_AMI_toggle=$graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{AMI},-text=>"AMI",-state=>"disabled")->pack(-anchor=>"w",-side=>"top");   

    my $graph_BS_toggle = $graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{BS},-text=>"BS")->pack(-anchor=>"w",-side=>"top");
    my $graph_PT_toggle = $graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{PT},-text=>"PT")->pack(-anchor=>"w",-side=>"top");

    my $graph_PD_toggle = $graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{PD},-text=>"PD")->pack(-anchor=>"w",-side=>"top");
    my $graph_B_toggle = $graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{B},-text=>"B")->pack(-anchor=>"w",-side=>"top");

    my $graph_PP_toggle = $graph_toggle_left_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{PP},-text=>"PP")->pack(-anchor=>"w",-side=>"top");
#    my $graph_B_toggle = $graph_toggle_right_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{B},-text=>"B")->pack(-anchor=>"w",-side=>"top");

    my $graph_wins_frame=$graph_frame->Frame(-background=>$$sheet{default_win_background})->pack(-side=>"top");
    my $graph_wins_desc=$graph_wins_frame->Label(-background=>$$sheet{default_win_background},-text=>"Window length")->pack(-anchor=>"w",-side=>"left");
    my $graph_wins_entry=$graph_wins_frame->Entry(-background=>$$sheet{default_win_background},-textvariable=>\$$plotstate{windowSize})->pack(-side=>"left");

    my $graph_slide_toggle=$graph_frame->Checkbutton(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-variable=>\$$plotstate{slide},-text=>"Sliding window")->pack(-anchor=>"w",-side=>"top");

    my $graph_action_frame=$graph_frame->Frame(-background=>$$sheet{default_win_background})->pack(-anchor=>"w",-side=>"top",-fill=>"x",-pady=>5,-ipady=>3);

    my $graph_button_apply;    
    my $graph_button_ok=$graph_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Ok",-command=>sub {
	$graph_button_apply->invoke();
	$graph_win->destroy;
    })->pack(-side=>"left");
    $graph_button_apply = $graph_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Apply",-command=>sub {
	calculate_graph($main,$canvas,$seq,$plotstate,$sheet);
	$$sheet{status}="dirty"; # Sheet is now dirty.
    })->pack(-side=>"left");
    my $graph_button_cancel=$graph_action_frame->Button(-activebackground=>$$sheet{dialog_button_active_background},-background=>$$sheet{default_win_background},-text=>"Cancel",-command=>sub {$graph_win->destroy;})->pack(-side=>"right");
}

sub heuristic_likely_orf {
    my $annotation = shift;
    my $seq = shift;
    my $sheet = shift;
    # also "secretly" uses $main and $status_text to inform the user

    my $seq_len=length($$seq);

    my @orfnr=();
    my @orflen=();
    my @orfgc=();
    my @strandpq=();

    my $norm_orf_len = 800;
    my $w_gc = 3;
    my $w_len = 10;
    my $w_strand = 5;
    my $w_strand_influence = 2200;

    my $maxIters=8;

    # get the numbers of ORF type annotations
    for ($nr=0; $nr < $$annotation{nr}; $nr++) {
	if($$annotation{type}->[$nr] eq 'ORF') {
	    push @orfnr, $nr;
	}
    }

    # sort the annotations with respect to stop position to speed up the strandedness calculation
    @orfnr = sort { $$annotation{stop}->[$a] <=> $$annotation{stop}->[$b] } @orfnr;

    # initialise
    foreach $nr (@orfnr) {
	$len = $$annotation{stop}->[$nr]-$$annotation{start}->[$nr]+1;
	push @orflen, $len;

	if($$annotation{frame}->[$nr] < 3) {
	    $_ = substr($$seq, $$annotation{start}->[$nr], $$annotation{stop}->[$nr]-$$annotation{start}->[$nr]+1);
	} else {
	    # revcomp
	    $_ = substr($$seq, $$annotation{start}->[$nr], $$annotation{stop}->[$nr]-$$annotation{start}->[$nr]+1 );
	    @_ = split(/ */,$_);
	    $_ = join('',reverse(@_));
	    tr/acgtACGT/tgcaTGCA/;
	}
	( @gc ) = /([gcGC]{1})/g;
	push @orfgc, @gc/$len;
	push @strandpq, 0;
    }

    # a first round of strandedness updating with strandpq == 1
   
    @strandedness = split(/ */, "0" x $seq_len);
    @local_cover = @strandedness;

    for(my $j = 0; $j < @orfnr; $j++) {	
	my $start = $$annotation{start}->[$orfnr[$j]];
	if($start < 0) {
	    $start = 0;
	}
	my $stop = $$annotation{stop}->[$orfnr[$j]];
	if($stop > $seq_len - 1) {
	    $stop = $seq_len - 1;
	}

	if($$annotation{frame}->[$orfnr[$j]] < 3) {
	    for($k = $start; $k < $stop; $k++) {
		$strandedness[$k] += 1; # length normalize? well, thats done in $strandpq calc, so I think we can skip it, right?
		$local_cover[$k]++;
	    }
	} else {
	    for($k = $start; $k < $stop; $k++) {
		$strandedness[$k] -= 1;
		$local_cover[$k]++;
	    }
	}
    }

    for($i = 0; $i < $seq_len; $i++) {
	if($local_cover[$i] != 0) {
	    $strandedness[$i] /= $local_cover[$i];
	}
    }


    # keep user posted on progress
    $status_text = "Removing ORF type annotations according to heuristic - computing strandedness..";
    $main->update(); # havent got a main ref!

    # iterate until strandedness doesn't change
    my @old_strandedness=();
    my $changes = 1;
    my $iterations = 0;
    while(($changes == 1) && ($iterations<$maxIters)) {
	$changes = 0;

	# update strandpq and likely
	for(my $j=0; $j < @orfnr; $j++) {
	    $strandpq[$j]=0;
	    for(my $k = $$annotation{start}->[$orfnr[$j]]; $k <= $$annotation{stop}->[$orfnr[$j]]; $k++) {
		$strandpq[$j] += $strandedness[$k];
	    }

	    if ($$annotation{frame}->[$orfnr[$j]] < 3 && $strandpq[$j] > 0) {
		$strandpq[$j] = $strandpq[$j] / $orflen[$j]; # aligned forward
	    } elsif ($$annotation{frame}->[$orfnr[$j]] >= 3 && $strandpq[$j] < 0) {
		$strandpq[$j] = - $strandpq[$j] / $orflen[$j]; # aligned reverse
	    } elsif ($$annotation{frame}->[$orfnr[$j]] < 3 && $strandpq[$j] < 0) {
		$strandpq[$j] = - $strandpq[$j] / $orflen[$j]; # opposing forward-oriented
	    } elsif ($$annotation{frame}->[$orfnr[$j]] >= 3 && $strandpq[$j] > 0) {
		$strandpq[$j] =  $strandpq[$j] / $orflen[$j]; # opposing reverse-oriented
	    }
	    
	    $likely[$orfnr[$j]] = ($w_gc * $orfgc[$j] + $w_len * $orflen[$j] / $norm_orf_len + $w_strand * $strandpq[$j])/($w_gc+$w_len+$w_strand); # we could insert data from ESTORF, Glimmer2, testcode, BLAST here if we like..
#	    $likely[$orfnr[$j]] = ($w_gc * $orfgc[$j] + $w_len * $orflen[$j] / $norm_orf_len) / ($w_gc+$w_len);
	    $DEBUG && print "DEBUG: ORF $$annotation{id}->[$orfnr[$j]] is $orflen[$j] long and contains $orfgc[$j]\% GC. The strand quotient (annotation frame is $$annotation{frame}->[$orfnr[$j]]) for the ORF is $strandpq[$j] so the heuristic measure is $likely[$orfnr[$j]].\n";
	}

	# update strandedness

 	@old_strandedness = @strandedness;

  	@strandedness = split(/ */, "0" x $seq_len);
	@local_cover = @strandedness;

  	for(my $j = 0; $j < @orfnr; $j++) {

	    $influence_region = POSIX::floor ( $likely[$orfnr[$j]] * $w_strand_influence);
	    my $start = $$annotation{start}->[$orfnr[$j]] - $influence_region;
	    if($start < 0) {
		$start = 0;
	    }
	    my $stop = $$annotation{stop}->[$orfnr[$j]] + $influence_region;
	    if($stop > $seq_len - 1 ) {
		$stop = $seq_len - 1;
	    }

  	    if($$annotation{frame}->[$orfnr[$j]] < 3) {
		for(my $k = $start; $k < $stop; $k++) {
		    $strandedness[$k] += $likely[$orfnr[$j]]; # length normalize? well, thats done in $strandpq calc, so I think we can skip it, right?
		    $local_cover[$k]++;
		}
	    } else {
		for(my $k = $start; $k < $stop; $k++) {
		    $strandedness[$k] -= $likely[$orfnr[$j]];
		    $local_cover[$k]++;
		}
	    }
  	} 

	for(my $i = 0; $i < $seq_len; $i++) {
	    if($local_cover[$i] != 0) {
		$strandedness[$i] /= $local_cover[$i];
	    }
	}

  	if (join('', @old_strandedness) ne join('',@strandedness)) {
  	    $changes = 1;
  	    $DEBUG && print "DEBUG: Update of strandedness ocurred..\n";
  	}

	$iterations++;

	# keep user posted on progress
	$DEVEL && print "DEVEL: Progress to iteration $iterations; changes is $changes.\n";

	$status_text = "Removing ORF type annotations according to heuristic - iteration $iterations of $maxIters..";
	$main->update(); # havent got a main ref!
    }

    # DEBUG - add results to notes..
    #for($j=0; $j < @orfnr; $j++) {
#	$DEVEL && $$annotation{note}->[$orfnr[$j]] .= "ORF removal heuristic final likely $likely[$orfnr[$j]]; length $orflen[$j], gc $orfgc[$j] strandpq is $strandpq[$j].";
 #   }

    # now go through the ORFs, starting with the highest ranking, and resolve any overlaps
    @orfscoresorted = reverse ( sort { $likely[$a] <=> $likely[$b] } @orfnr );
    
    my @all_overlapping = ();
    my @nonoverlapping = @orfnr;

    my $orfsdone = 0;
    my $j = 0;
    while ($orfsdone < @orfscoresorted) {
	my @overlapping = (); # overlaps for this orf, $orfscoresorted[$j]
	# check if current scoresorted is left non-overlapping
	$_= join(',',@nonoverlapping);
	$_ = ",$_,";
	if(!(m/,$orfscoresorted[$j],/)) {
	    $j++;
	    next;
	}

	foreach my $k (@nonoverlapping) {
	    if( (($$annotation{start}->[$orfscoresorted[$j]] > $$annotation{start}->[$k]) && ($$annotation{start}->[$orfscoresorted[$j]] < $$annotation{stop}->[$k]) )||
		(($$annotation{stop}->[$orfscoresorted[$j]] > $$annotation{start}->[$k]) && ($$annotation{stop}->[$orfscoresorted[$j]] < $$annotation{stop}->[$k])) ||
		(($$annotation{start}->[$orfscoresorted[$j]] < $$annotation{start}->[$k]) && ($$annotation{stop}->[$orfscoresorted[$j]] > $$annotation{stop}->[$k])) ) {
		# Overlap over start, stop or overlap completely contained in ORF.
		push @overlapping, $k;
		$orfsdone++;
	    }
	}

	$DEBUG && print "DEBUG: The following annotations overlapped $$annotation{id}->[$orfscoresorted[$j]] (likely $likely[$orfscoresorted[$j]]) and are thus scheduled for removal: \n";

	# remove from nonoverlapping list
	my $nonolap= join(',',@nonoverlapping);
	$nonolap = ",$nonolap,";
	foreach my $olap (@overlapping) {
	    $DEBUG && print "DEBUG: Orfnr $olap: $$annotation{id}->[$olap] ($likely[$olap])\n";
	    # $$annotation{note}->[$olap] .= "\nOverlapped by $$annotation{id}->[$orfscoresorted[$j]] (likely $likely[$orfscoresorted[$j]]).";
	    $nonolap=~s/,$olap,/,/;
	    # $DEBUG && print "DEBUG: nonoverlapping left: $_\n";
            # s/,,/,/;
	}
	($nonolap) = $nonolap =~/^,(.+),$/;
	
	@nonoverlapping = split(/,/,$nonolap);
					    
	push @all_overlapping, @overlapping;
	$orfsdone++; # count done orfs up with the current orf  
	$j++; # counter for looping over high-scoring orfs
    }

    $DEBUG && print "DEBUG: Remaining, nonoverlapped or highest-scoring, annotations are @nonoverlapping.\n";

    # this is slightly controversial: clean out annotations if strandedness does not agree.. 
    for($j=0; $j < @nonoverlapping; $j++) {
	$start = $$annotation{start}->[$nonoverlapping[$j]];
	$stop = $$annotation{stop}->[$nonoverlapping[$j]];
	$frame = $$annotation{frame}->[$nonoverlapping[$j]];
	
	$strandedness_sum = sum ( @strandedness[$start .. $stop] );
	if(($strandedness_sum < 0 && $frame < 3) || ($strandedness_sum >= 0 && $frame >= 3)) {
	    push @all_overlapping , $nonoverlapping[$j];
	    $DEVEL && print "DEVEL: Annotation $$annotation{name}->[$nonoverlapping[$j]] has frame $frame when strandedness over its extent ($start - $stop) is $strandpq[$j] - scheduled for removal!\n";
	}
    }

    # remove overlapping annotations..
    remove_annotations($annotation, \@all_overlapping,$sheet);

    $nr_of_removed_orfs = @all_overlapping; 
    $status_text = "Done removing $nr_of_removed_orfs ORF type annotations according to heuristic."; 
}

sub calculate_graph {
    my $main=shift;
    my $canvas=shift;
    my $seqref=shift;
    my $plotstate=shift;
    my $sheet=shift;

    my $seq = $$seqref;
    $seq=~tr/ATGC/atgc/;
    my @alphabet=('a','t','g','c'); 
    
    my %histSeq=();
    my %p=();
    my $info=0;
    
    my @sequence=split(/ */,$seq);
    
    foreach(@sequence) {
	$histSeq{$_}++;
    }
    
    print "Overall frequencies:\n";
    
    foreach my $char (@alphabet) {
	$p{$char}=$histSeq{$char}/@sequence;
	
	print "$char $p{$char}\n";
	
	if($p{$char}!=0) {
	    $info-=( $p{$char}*log2($p{$char}) );
	}
    }
    
    print "Overall H=$info\n\n";
    
    # Output data is saved indexed on window position.
    my @HWin=();
    my %pWin=();
    
    my @xAxis=();
    
    my @atFraction=();
    my @gcFraction=();
    my @agFraction=();
    my @ctFraction=();
    
    my @gcSkew = ();
    my @atSkew = ();
    my @cumATskew = ();    
    my @cumGCskew = ();
    
    my $windowNr=0;
    
    if($$plotstate{slide}) {
	#
	# Run sliding windows
	#

	my %histWin=();
	
	@sequence=split(/ */,$seq);
	# Slide window
	my $i;
	for($i=0;$i<length($seq)-$$plotstate{windowSize}-1;$i++) {
	    
	    # Fill up the window initially...
	    if($i==0) {
		my @seqWin=split(/ */,substr($seq, 0, $$plotstate{windowSize}));
#		$DEBUG && print "DEBUG: Window $i has actual length ",scalar(@seqWin),".\n"; 
		foreach(@seqWin) {
		    $histWin{$_}++;
		}
	    } else {
		# Update window
		$histWin{$sequence[$i-1+$$plotstate{windowSize}]}++;
		$histWin{$sequence[$i-1]}--;
       	    }
	    
	    # Calculate window stats
	    $HWin[$i]=0;
	    foreach my $char (@alphabet) {
		if(defined($histWin{$char}) && ($histWin{$char} != 0 )) {
		    $pWin{$char}[$i]=$histWin{$char}/$$plotstate{windowSize};
		    # $DEBUG && print "DEBUG: Window $i char $char freq $histWin{$char}.\n";
		    $HWin[$i]-=( $pWin{$char}[$i]*log2($pWin{$char}[$i]) );
		} else {
		    # No such thing as log(0).
		    $pWin{$char}[$i]=0;
		    $DEVEL && print "DEVEL: 0 prob for $char in window $i!\n";
		}
	    }

	    $gcFraction[$i]=$pWin{'g'}[$i]+$pWin{'c'}[$i];
	    $atFraction[$i]=$pWin{'a'}[$i]+$pWin{'t'}[$i];
	    $agFraction[$i]=$pWin{'a'}[$i]+$pWin{'g'}[$i];
	    $ctFraction[$i]=$pWin{'t'}[$i]+$pWin{'c'}[$i];
	    $HWin[$i]/=2;	# NOTE! This is here simply to make the plot look nice...

	    if($$plotstate{GCs} || $$plotstate{ATs}) {
		if(($pWin{'g'}[$i]+$pWin{'c'}[$i]) != 0) {
		    $gcSkew[$i]=($pWin{'g'}[$i]-$pWin{'c'}[$i])/($pWin{'g'}[$i]+$pWin{'c'}[$i]);
		} else {
		    $WARNING && print "WARNING: P(G)+P(C) = 0 in window $i - setting skew to zero.\n";
		    $gcSkew[$i]=0;
		}
		
		if($i == 0) {
		    $cumGCskew[$i] = $gcSkew[$i];
		} else {
		    $cumGCskew[$i] = $cumGCskew[$i-1] + $gcSkew[$i];
		}
		
		if(($pWin{'a'}[$i]+$pWin{'t'}[$i]) != 0) {
		    $atSkew[$i]=($pWin{'a'}[$i]-$pWin{'t'}[$i])/($pWin{'a'}[$i]+$pWin{'t'}[$i]);
		} else {
		    $WARNING && print "WARNING: P(A)+P(T) = 0 in window $i - setting skew to zero.\n";
		    $atSkew[$i]=0;
		}
		
		if($i == 0) {
		    $cumATskew[$i] = $atSkew[$i];
		} else {			
		    $cumATskew[$i] = $cumATskew[$i-1] + $atSkew[$i];
		}
	    }
	    
	    $xAxis[$i]=$i+$$plotstate{windowSize}/2; # Where to put the dot is disputable, but why not centered around the window...
	}
	$windowNr=$i; # Cut & paste code... Not that nice; i might even do a S&R.. =)
	
    } else {
	
	#
	# Non-overlapping windows
	#

	for(my $i=0;$i<length($seq);$i+=$$plotstate{windowSize}) {
	    my @seqWin=split(/ */,substr($seq, $i, $$plotstate{windowSize}));
	    my $actualWinSize = @seqWin;
		
	    my %histWin=();
		
	    # Determine char frequencies in window
	    foreach(@seqWin) {
		$histWin{$_}++; # Might be a little low on nx and such..
	    }
	    
	    $HWin[$windowNr]=0;
	    foreach my $char (@alphabet) {
		# No such thing as log(0).
		    if(defined($histWin{$char}) && ($histWin{$char} != 0 )) {
			$pWin{$char}[$windowNr]=$histWin{$char}/$actualWinSize;
			$HWin[$windowNr]-=( $pWin{$char}[$windowNr]*log2($pWin{$char}[$windowNr]) );
		    } else {
			$pWin{$char}[$windowNr]=0;
		    }
		}
		
	    $atFraction[$windowNr]=$pWin{'a'}[$windowNr]+$pWin{'t'}[$windowNr];
	    $gcFraction[$windowNr]=$pWin{'g'}[$windowNr]+$pWin{'c'}[$windowNr];
	    $agFraction[$windowNr]=$pWin{'a'}[$windowNr]+$pWin{'g'}[$windowNr];
	    $ctFraction[$windowNr]=$pWin{'t'}[$windowNr]+$pWin{'c'}[$windowNr];
		
	    if($$plotstate{GCs} || $$plotstate{ATs}) {
		if(($pWin{'g'}[$windowNr]+$pWin{'c'}[$windowNr]) != 0) {
		    $gcSkew[$windowNr]=($pWin{'g'}[$windowNr]-$pWin{'c'}[$windowNr])/($pWin{'g'}[$windowNr]+$pWin{'c'}[$windowNr]);
		} else {
		    $WARNING && print "WARNING: P(G)+P(C) = 0 in window $windowNr - setting skew to zero.\n";
		    $gcSkew[$windowNr]=0;
		}
		
		if($windowNr == 0) {
		    $cumGCskew[$windowNr] = $gcSkew[$windowNr];
		} else {
			$cumGCskew[$windowNr] = $cumGCskew[$windowNr-1] + $gcSkew[$windowNr];
		    }
		
		if(($pWin{'a'}[$windowNr]+$pWin{'t'}[$windowNr]) != 0) {
		    $atSkew[$windowNr]=($pWin{'a'}[$windowNr]-$pWin{'t'}[$windowNr])/($pWin{'a'}[$windowNr]+$pWin{'t'}[$windowNr]);
		} else {
		    $WARNING && print "WARNING: P(A)+P(T) = 0 in window $windowNr - setting skew to zero.\n";
		    $atSkew[$windowNr]=0;
		}
		
		if($windowNr == 0) {
		    $cumATskew[$windowNr] = $atSkew[$windowNr];
		} else {			
		    $cumATskew[$windowNr] = $cumATskew[$windowNr-1] + $atSkew[$windowNr];
		}
	    }
	    
	    $HWin[$windowNr]/=2;	# NOTE! This normalisation is here simply to make the plot look nice...      
	    $xAxis[$windowNr]=$i+$actualWinSize/2; # + windowsize/2 for a q&d centered window...
	    $windowNr++;
	}        
    }
    
    if($$plotstate{GC}) {
	@{$$plotstate{gcFraction}} = @gcFraction;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    if($$plotstate{AT}) {
	@{$$plotstate{atFraction}} = @atFraction;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    if($$plotstate{AG}) {
	@{$$plotstate{agFraction}} = @agFraction;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    if($$plotstate{CT}) {
	@{$$plotstate{ctFraction}} = @ctFraction;
	@{$$plotstate{xAxis}} = @xAxis;
    }
    
    if($$plotstate{GCs}) {
	my $max_cum_GCs = max(@cumGCskew);
	my $min_cum_GCs = min(@cumGCskew);
	my $norm_const;	

	if(abs($min_cum_GCs) > abs($max_cum_GCs)) {
	    $norm_const = abs($min_cum_GCs) * 2;
	} else {
	    $norm_const = abs($max_cum_GCs) * 2;
	}

	my @normCumGCskew = map { $_ / $norm_const + 0.5} @cumGCskew;

	@{$$plotstate{normCumGCskew}} = @normCumGCskew;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    if($$plotstate{ATs}) {
	my $max_cum_ATs = max(@cumATskew);
	my $min_cum_ATs = min(@cumATskew);
	my $norm_const;

	if(abs($min_cum_ATs) > abs($max_cum_ATs)) {
	    $norm_const = abs($min_cum_ATs) * 2;
	} else {
	    $norm_const = abs($max_cum_ATs) * 2;
	}

	my @normCumATskew = map { $_ / $norm_const + 0.5} @cumATskew;
	
	@{$$plotstate{normCumATskew}} = @normCumATskew;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    if($$plotstate{H} ) {  
	@{$$plotstate{HWin}} = @HWin;
	@{$$plotstate{xAxis}} = @xAxis;
    }

    my %dinuc;
    my @x;
    my $points;

    if( $$plotstate{BS} || $$plotstate{PT} || $$plotstate{PD} ) {
	dinucleotide_descriptor($seqref, \%dinuc, $$plotstate{windowSize},\@x);
	$points = @x;
    }

    if($$plotstate{BS} ) {  
	@{$$plotstate{dinucBS}} = @{$dinuc{BS}};
	@{$$plotstate{xAxis}} = @x;
    }
    
    if($$plotstate{PT} ) {
	@{$$plotstate{dinucPT}} = @{$dinuc{PT}}; 
	@{$$plotstate{xAxis}} = @x;
    }
    
    if($$plotstate{PD} ) {
	@{$$plotstate{dinucPD}} = @{$dinuc{PD}};
	@{$$plotstate{xAxis}} = @x;
    }

    my %trinuc;

    if( $$plotstate{B} || $$plotstate{PP} ) {
	trinucleotide_descriptor($seqref, \%trinuc, $$plotstate{windowSize}, \@x); 
	$points = @x;
    }

    if( $$plotstate{B} ) {
	@{$$plotstate{trinucB}} = @{$trinuc{B}};
	@{$$plotstate{xAxis}} = @x;
    }

    if($$plotstate{PP} ) {
	@{$$plotstate{trinucPP}} = @{$trinuc{PP}};
	@{$$plotstate{xAxis}} = @x;
    }
    	
    # Update legend? Keep legend in a separate window??  

    replot($canvas, $seqref, $plotstate, $sheet);
}

sub replot {   
    my $canvas = shift;
    my $seq = shift;
    my $plotstate = shift;
    my $sheet = shift;

    $WARNING && (ref ($plotstate) || print "WARNING: plotstate is not a reference in replot!\n");

    my $windowNr = (defined @{$$plotstate{xAxis}}) ? @{$$plotstate{xAxis}} : 0;

    # Create axis if needed.
    if($$plotstate{axis} == 0) {
	$$canvas->delete(-tags=>"plotAxis");
	plot_axis($canvas,1,length($$seq),$sheet);
	$$plotstate{axis}=1;
    }

    # Plot selected curves, and delete unselected old curves.
    
    # Delete the old curve of this type if it exist before plotting a new 
    # with e g a potentially different window size.

    if($$plotstate{GC}) { 
	$$canvas->delete(-tags=>"gcfrac");
	my $plotname="gcfrac";
	$DEBUG && print "DEBUG: attemplting to run plot_curve from replot of GC content. Xaxis is ".scalar(@{$$plotstate{xAxis}})."and the actual gc-data is ".scalar(@{$$plotstate{gcFraction}}).".\n";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{gcFraction}}, $plotname,'red1',$sheet);
    } else {
#  $DEBUG && print "DEBUG: deleting gcplot upon request..\n"
	$$canvas->delete(-tags=>"gcfrac");
    }

    if($$plotstate{AT} ) {
	$$canvas->delete(-tags=>"atfrac");
	my $plotname="atfrac";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{atFraction}}, $plotname,'red3',$sheet);
    } else {
	$$canvas->delete(-tags=>"atfrac");
    }
    if($$plotstate{AG} ) {
	$$canvas->delete(-tags=>"agfrac");
	my $plotname="agfrac";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{agFraction}}, $plotname,'green1',$sheet);
    } else {
	$$canvas->delete(-tags=>"agfrac");
    }
    if($$plotstate{CT} ) {
	$$canvas->delete(-tags=>"ctfrac");
	my $plotname="ctfrac";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{ctFraction}}, $plotname,'green3',$sheet);
    } else {
	$$canvas->delete(-tags=>"ctfrac");
    }
    
    if($$plotstate{GCs}) {
	$$canvas->delete(-tags=>"gcskew");
	my $plotname="gcfrac";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{normCumGCskew}}, $plotname,'magenta',$sheet);
    } else {
	$$canvas->delete(-tags=>"gcskew");
    }
    if($$plotstate{ATs}) {
	$$canvas->delete(-tags=>"atskew");
	my $plotname="atfrac";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{normCumATskew}}, $plotname,'cyan',$sheet);
    } else {
	$$canvas->delete(-tags=>"atskew");
    }

    if($$plotstate{H} ) {  
	$$canvas->delete(-tags=>"hwin");    
	my $plotname="hwin";
	plot_curve($canvas,$windowNr,\@{$$plotstate{xAxis}},\@{$$plotstate{HWin}}, $plotname,'blue',$sheet);    
    } else {
	$$canvas->delete(-tags=>"hwin");
    }

    if($$plotstate{BS} ) {  
	$$canvas->delete(-tags=>"bsprop");
	my $plotname="bsprop";
	plot_curve($canvas, $windowNr,\@{$$plotstate{xAxis}}, \@{$$plotstate{dinucBS}}, $plotname,'blue2',$sheet);
    } else {
	$$canvas->delete(-tags=>"bsprop");
    }
    
    if($$plotstate{PT} ) {
	$$canvas->delete(-tags=>"ptprop");
	my $plotname="ptprop";
	plot_curve($canvas, $windowNr, \@{$$plotstate{xAxis}}, \@{$$plotstate{dinucPT}}, $plotname,'blue4',$sheet);
    } else {
	$$canvas->delete(-tags=>"ptprop");
    }
    
    if($$plotstate{PD} ) {
	$$canvas->delete(-tags=>"pdprop");
	my $plotname="pdprop";
	plot_curve($canvas, $windowNr, \@{$$plotstate{xAxis}}, \@{$$plotstate{dinucPD}}, $plotname,'orange1',$sheet);
    } else {
	$$canvas->delete(-tags=>"pdprop");
    }

    if( $$plotstate{B} ) {
	$$canvas->delete(-tags=>"bprop");
	my $plotname="bprop";
	# $DEVEL && print "DEVEL: Got B @{$trinuc{B}}.\n";
	plot_curve($canvas, $windowNr, \@{$$plotstate{xAxis}}, \@{$$plotstate{trinucB}}, $plotname,'orange3',$sheet);
    } else {
	$$canvas->delete(-tags=>"bprop");
    }

    if($$plotstate{PP} ) {
	$$canvas->delete(-tags=>"ppprop");
	my $plotname="ppprop";
	plot_curve($canvas, $windowNr, \@{$$plotstate{xAxis}}, \@{$$plotstate{trinucPP}}, $plotname,'green2',$sheet);
    } else {
	$$canvas->delete(-tags=>"ppprop");
    }
        
    if(!$$plotstate{GC} and !$$plotstate{AT} and !$$plotstate{AG} and !$$plotstate{CT} and !$$plotstate{GCs} and !$$plotstate{ATs} and !$$plotstate{H} and !$$plotstate{BS} and !$$plotstate{PT} and !$$plotstate{PD} and !$$plotstate{B} and !$$plotstate{PP} and $$plotstate{axis}) {
	# Nothing to show - remove axis, and shrink canvas.
	$$canvas->delete(-tags=>"plotAxis");
	$$plotstate{axis}=0;
	$$canvas->configure(-height=>$$sheet{canvas_seq_height});
    }

}

sub view_legend {
    my $main = shift;
    my $annotation = shift;
    my $plotstatus = shift;
    my $sheet = shift;

    my $legend=$$main->Toplevel;
    $legend->title("A GUI - Legend");
    $legend->configure(-width=>'300',-height=>'340',-background=>$$sheet{default_win_background});

    my $legend_main_frame=$legend->Frame(-background=>$$sheet{default_win_background})->pack(-expand=>'yes',-fill=>'both');
#    my $legend_label=$legend_main_frame->Label(-text=>'Annotation',-background=>$$sheet{default_win_background})->pack(-side=>'top',-anchor=>'w');
    my $legend_canvas=$legend_main_frame->Canvas(-background=>$$sheet{default_win_background})->pack(-expand => 'yes', -fill => 'both', -side=>'top');
    
    # Check for actually present annotations?

    my ($x1,$x2,$x3)=(10,55,70);
    my $y=15;

    my $heading1=$legend_canvas->createText($x1,$y,-text=>"Annotations",-anchor=>'w');
    $y+=15;
    my $arrow1 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'red1', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label1 = $legend_canvas->createText($x3,$y,-text=>"Glimmer 2 putative",-anchor=>'w');
    $y+=15;
    my $arrow2 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'red3', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label2 = $legend_canvas->createText($x3,$y,-text=>"Glimmer 2 uncertain putative",-anchor=>'w');
    $y+=15;
    my $arrow3 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'blue1', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label3 = $legend_canvas->createText($x3,$y,-text=>"cDNA tagged ORF - high score",-anchor=>'w');
    $y+=15;
    my $arrow4 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'blue3', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label4 = $legend_canvas->createText($x3,$y,-text=>"cDNA tagged ORF - intermediate score",-anchor=>'w');
    $y+=15;
    my $box1 = $legend_canvas->createLine($x1,$y,$x2,$y, -width=>10,-fill=>'yellow1');
    my $label5 = $legend_canvas->createText($x3,$y,-text=>"cDNA blast hit - high score",-anchor=>'w');
    $y+=15;
    my $box2 = $legend_canvas->createLine($x1,$y,$x2,$y, -width=>10,-fill=>'yellow3');
    my $label6 = $legend_canvas->createText($x3,$y,-text=>"cDNA blast hit - intermediate score",-anchor=>'w');
    $y+=15;
    my $arrow5 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'gray', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label7 = $legend_canvas->createText($x3,$y,-text=>"Open Reading Frame",-anchor=>'w');
    $y+=15;
    my $line1 = $legend_canvas->createLine($x1, $y, $x1+3, $y, -width=>10,-fill=>'green');
    my $line2 = $legend_canvas->createLine($x2-3, $y, $x2, $y, -width=>10,-fill=>'black');
    my $label8 = $legend_canvas->createText($x3,$y,-text=>"Start and stop codons",-anchor=>'w');
    $y+=15;
    my $arrow6 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'magenta', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label9 = $legend_canvas->createText($x3,$y,-text=>"Manually added with orientation",-anchor=>'w');
    $y+=15;
    my $box3 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'magenta');
    my $label10 = $legend_canvas->createText($x3,$y,-text=>"Manually added without orientation",-anchor=>'w');
    $y+=15;
    my $arrow7 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'purple', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label11 = $legend_canvas->createText($x3,$y,-text=>"Merged annotation",-anchor=>'w');
    $y+=15;
    my $arrow8 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'green3', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label12 = $legend_canvas->createText($x3,$y,-text=>"Poly pyrimidine stretch",-anchor=>'w');
    $y+=15;
    my $arrow9 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'DarkOrange1', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label13 = $legend_canvas->createText($x3,$y,-text=>"Blast hit",-anchor=>'w');
    $y+=15;
    my $arrow10 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'SkyBlue3', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label13b = $legend_canvas->createText($x3,$y,-text=>"Testcode positive ORF",-anchor=>'w');
    $y+=15;
    my $arrow11 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'dark turquoise', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label13c = $legend_canvas->createText($x3,$y,-text=>"Genbank CDS",-anchor=>'w');
    $y+=15;
    my $arrow12 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'aquamarine3', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label13d = $legend_canvas->createText($x3,$y,-text=>"General Feature Format",-anchor=>'w');
    $y+=15;
    my $arrow13 = $legend_canvas->createLine($x1,$y,$x2,$y , -width=>10,-fill=>'Hot Pink', -arrow=>'last', -arrowshape=>[5,5,0]);
    my $label13e = $legend_canvas->createText($x3,$y,-text=>"Regexp match",-anchor=>'w');

    # Check for actually active plots?
    $y+=30;
    my $heading2=$legend_canvas->createText($x1,$y,-text=>"Graph",-anchor=>'w');
    $y+=15;
    my $graph1 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'red1');
    my $label14 = $legend_canvas->createText($x3,$y,-text=>"GC content",-anchor=>'w');
    $y+=15;
    my $graph2 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'red3');
    my $label15 = $legend_canvas->createText($x3,$y,-text=>"AT content",-anchor=>'w');
    $y+=15;
    my $graph3 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'green1');
    my $label16 = $legend_canvas->createText($x3,$y,-text=>"AG (Purine) content",-anchor=>'w');
    $y+=15;
    my $graph4 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'green3');
    my $label17 = $legend_canvas->createText($x3,$y,-text=>"CT (Pyrimidine) content ",-anchor=>'w');
    $y+=15;
    my $graph5 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'blue');
    my $label18 = $legend_canvas->createText($x3,$y,-text=>"H/2 Information content (Source entropy)",-anchor=>'w');
    $y+=15;
    my $graph6 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'magenta');
    my $label19 = $legend_canvas->createText($x3,$y,-text=>"GC skew (cumulative)",-anchor=>'w');
    $y+=15;
    my $graph7 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'cyan');
    my $label20 = $legend_canvas->createText($x3,$y,-text=>"AT skew (cumulative)",-anchor=>'w');
    $y+=15;
    my $graph8 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'blue2');
    my $label21 = $legend_canvas->createText($x3,$y,-text=>"Base stacking (dinucleotide measure)",-anchor=>'w');
    $y+=15;
    my $graph9 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'orange1');
    my $label22 = $legend_canvas->createText($x3,$y,-text=>"Protein deformability (dinucleotide measure)",-anchor=>'w');
    $y+=15;
    my $graph10 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'blue4');
    my $label23 = $legend_canvas->createText($x3,$y,-text=>"Propeller twist (dinucleotide measure)",-anchor=>'w');
    $y+=15;
    my $graph11 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'green2');
    my $label24 = $legend_canvas->createText($x3,$y,-text=>"Position preference (trinucleotide measure)",-anchor=>'w');
    $y+=15;
    my $graph12 = $legend_canvas->createLine($x1,$y+5,$x2,$y-5,-fill=>'orange3');
    my $label25 = $legend_canvas->createText($x3,$y,-text=>"Bendability (trinucleotide measure)",-anchor=>'w');
}

sub view_stats {
    my $annotation=shift;
    my $seqName=shift;
    my $seq=shift;

    my $nr=$$annotation{nr};

    my %freq;

    for(my $i=0;$i<$nr;$i++) {       
	$freq{$$annotation{type}->[$i]}++; # Count number of annotations for each type		
    }

    print "Currently, the display contains :\n";
    foreach my $type (keys %freq) {
	print "$annotatorName{$type} : $freq{$type}\n"; 
    }
}

sub run {
    my $main=shift; 
    my $canvas=shift;
    my $status_text=shift; #ref
    my $annotation=shift; #ref
    my $seqName=shift; #ref 
    my $seq=shift;#ref
    my $sheet=shift; #ref
    my $first_tab=shift;

    my %runstate=();
    
    $runstate{glim_shortest_orf} = 100;

    my $run_win=$$main->Toplevel;
    $run_win->title("A GUI - run");
    # no geometry specified - actually the placement feels fairly ok now?!

    my $run_main_frame=$run_win->Frame()->pack(-expand=>'y',-fill=>'both');
    
    my $run_notebook=$run_main_frame->NoteBook(-ipadx=>6, -ipady=>6)->pack();
    # Note: a page in a notebook is a Frame reference...

    # GLIMMER2 run page

    my $run_glimmer_page=$run_notebook->add("glimmer", -label => "Glimmer2", -underline=>0);

    my $r_glim_run_label=$run_glimmer_page->Label(-text=>"Run glimmer on")->pack(-anchor=>'w',-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_run_seq_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_run_seq_cb=$r_glim_run_seq_frame->Checkbutton(-text=>"Current sequence $$seqName",-variable=>\$runstate{glim_run_seq},-command=>sub{ $runstate{glim_run}=0; })->pack(-side=>'left');
    my $r_glim_run_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_run_cb=$r_glim_run_frame->Checkbutton(-text=>'File',-variable=>\$runstate{glim_run},-command=>sub{ $runstate{glim_run_seq}=0; })->pack(-side=>'left');
    my $r_glim_run_entry=$r_glim_run_frame->Entry(-textvariable=>\$runstate{glim_run_file})->pack(-side=>'left');
    my $r_glim_run_browse=$r_glim_run_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['Sequence FASTA file','*']);
	$runstate{glim_run_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Run glimmer on");
    })->pack(-side=>'right');
    
    my $r_glim_model_label=$run_glimmer_page->Label(-text=>"Use model")->pack(-anchor=>'w',-side=>'top',-expand=>'y',-fill=>'x');

    my $r_glim_train_seq_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_train_seq_cb=$r_glim_train_seq_frame->Checkbutton(-text=>"Train on current sequence $$seqName",-variable=>\$runstate{glim_train_seq},-command=>sub{ $runstate{glim_train}=0; $runstate{glim_model}=0;})->pack(-side=>'left');

    my $r_glim_train_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_train_cb=$r_glim_train_frame->Checkbutton(-text=>'Train with file',-variable=>\$runstate{glim_train},-command=>sub{ $runstate{glim_train_seq}=0; $runstate{glim_model}=0;})->pack(-side=>'left');
    my $r_glim_train_entry=$r_glim_train_frame->Entry(-textvariable=>\$runstate{glim_train_file})->pack(-side=>'left');
    my $r_glim_train_browse=$r_glim_train_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['Sequence FASTA file','*']);
	$runstate{glim_train_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Build ICM from sequence file");
    })->pack(-side=>'right');

    my $r_glim_model_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_glim_model_cb=$r_glim_model_frame->Checkbutton(-text=>'Use model file',-variable=>\$runstate{glim_model},-command=>sub{ $runstate{glim_train}=0; $runstate{glim_train_seq}=0;})->pack(-side=>'left');
    my $r_glim_model_entry=$r_glim_model_frame->Entry(-textvariable=>\$runstate{glim_model_file})->pack(-side=>'left');
    my $r_glim_model_browse=$r_glim_model_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['ICM model file','*']);
	$runstate{glim_model_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select ICM model");
    })->pack(-side=>'right');

    my $r_glim_options_label=$run_glimmer_page->Label(-text=>"Options")->pack(-anchor=>'w',-side=>'top',-expand=>'y',-fill=>'x');
    
    my $r_glim_altstart_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $r_glim_altstart_cb=$r_glim_altstart_frame->Checkbutton(-state=>'disabled',-text=>"Use alternative start codons",-variable=>\$runstate{glim_altstart})->pack(-side=>'left');
    my $r_glim_shortest_orf_frame=$run_glimmer_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $r_glim_shortest_orf_label=$r_glim_shortest_orf_frame->Label(-text=>"Shortest ORF")->pack(-side=>'left');
    my $r_glim_shortest_orf_entry=$r_glim_shortest_orf_frame->Entry(-textvariable=>\$runstate{glim_shortest_orf})->pack(-side=>'left');

    my $run_estorf_page=$run_notebook->add("ESTORF", -label => "EST ORF", -underline=>4);

    my $run_estorf_label=$run_estorf_page->Label(-text=>'Find cDNA matches to sequence, as well as cDNA tagged ORFs.')->pack(-side=>'top',-anchor=>'w');
    my $run_estorf_target_label=$run_estorf_page->Label(-text=>"Subject sequence")->pack(-side=>'top',-anchor=>'w');
    my $run_estorf_currentseq_rb=$run_estorf_page->Radiobutton(-text=>"Current sequence $$seqName",-variable=>\$runstate{estorf_currentseq}, -value=>'1', -command=>sub{ $runstate{estorf_seq} = 0;})->pack(-side=>'top',-anchor=>'w');
    my $run_estorf_seq_frame=$run_estorf_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_estorf_seq_rb=$run_estorf_seq_frame->Radiobutton(-text=>'Use sequence file',-variable=>\$runstate{estorf_seq}, -value=>'1',-command=>sub{ $runstate{estorf_currentseq} = 0;})->pack(-side=>'left');
    my $run_estorf_seq_entry=$run_estorf_seq_frame->Entry(-textvariable=>\$runstate{estorf_seq_file})->pack(-side=>'left');
    my $run_estorf_seq_browse=$run_estorf_seq_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['Sequence FASTA file','*']);
	$runstate{estorf_seq_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select FASTA sequence file");
    })->pack(-side=>'right');
    my $run_estorf_est_label=$run_estorf_page->Label(-text=>"cDNA sequences")->pack(-side=>'top',-anchor=>'w');  
    my $run_estorf_est_frame=$run_estorf_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_estorf_est_cb=$run_estorf_est_frame->Checkbutton(-text=>'Use sequence file',-variable=>\$runstate{estorf_est})->pack(-side=>'left');
    my $run_estorf_est_entry=$run_estorf_est_frame->Entry(-textvariable=>\$runstate{estorf_est_file})->pack(-side=>'left');
    my $run_estorf_est_browse=$run_estorf_est_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['EST sequence file','*']);
	$runstate{estorf_est_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select EST FASTA file");})->pack(-side=>'right');

    my $run_estorf_blast_label = $run_estorf_page->Label(-text=>"Local blast package")->pack(-side=>'top',-anchor=>'w');  
    my $run_estorf_blast_frame = $run_estorf_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_estorf_blast_rb_ncbi = $run_estorf_blast_frame->Radiobutton(-text=>'NCBI (megablast)',-variable=>\$runstate{estorf_ncbi}, -value=>'1',-command=>sub { $runstate{estorf_wu} = 0; })->pack(-side=>'left');
    my $run_estorf_blast_rb_wu = $run_estorf_blast_frame->Radiobutton(-text=>'WU-BLAST2 (blastn)',-variable=>\$runstate{estorf_wu}, -value=>'1',-command=>sub { $runstate{estorf_ncbi} = 0; })->pack(-side=>'left');

#  my $run_estorf_db_frame=$run_estorf_page->Frame()->pack(-side=>'top',-fill=>'x');
#  my $run_estorf_db_cb=$run_estorf_db_frame->Checkbutton(-text=>'Use blastable db file',-variable=>\$runstate{estorf_db})->pack(-side=>'left');
#  my $run_estorf_db_entry=$run_estorf_db_frame->Entry(-textvariable=>\$runstate{estorf_db_file})->pack(-side=>'left');
#  my $run_estorf_db_browse=$run_estorf_db_frame->Button(-text=>'Browse...',-command=>sub {
#						      my @filetypes=(['WU-BLAST 2 db file','*']);
#						      $runstate{estorf_db_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select WU-BLAST 2 db file");})->pack(-side=>'right');

    my $run_blast_page=$run_notebook->add("blast", -label => "BLAST", -underline=>0);

    my $run_blast_label=$run_blast_page->Label(-text=>'Run blast')->pack(-side=>'top',-anchor=>'w');
    my $run_blast_target_label=$run_blast_page->Label(-text=>"Subject sequence")->pack(-side=>'top',-anchor=>'w');
    my $run_blast_currentseq_rb=$run_blast_page->Radiobutton(-text=>"Current sequence $$seqName",-variable=>\$runstate{blast_currentseq}, -value=>'1', -command=>sub{ $runstate{blast_seq} = 0;})->pack(-side=>'top',-anchor=>'w');
    my $run_blast_seq_frame=$run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');

    my $run_blast_seq_rb=$run_blast_seq_frame->Radiobutton(-text=>'Use sequence file',-variable=>\$runstate{blast_seq}, -value=>'1',-command=>sub{ $runstate{blast_currentseq} = 0;})->pack(-side=>'left');
    my $run_blast_seq_entry=$run_blast_seq_frame->Entry(-textvariable=>\$runstate{blast_seq_file})->pack(-side=>'left');
    my $run_blast_seq_browse=$run_blast_seq_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['Sequence FASTA file','*']);
	$runstate{blast_seq_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select FASTA sequence file");
    })->pack(-side=>'right');

    # mutually exclusive options! add support..
    my $run_blast_db_label=$run_blast_page->Label(-text=>"Database")->pack(-side=>'top',-anchor=>'w');

    my $run_blast_db_frame=$run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_db_cb=$run_blast_db_frame->Checkbutton(-text=>'Local blast db',-variable=>\$runstate{blast_db}, -command=>sub {$runstate{blast_db_seq}=0; $runstate{blast_db_remote}=0;})->pack(-side=>'left');
    my $run_blast_db_entry=$run_blast_db_frame->Entry(-textvariable=>\$runstate{blast_db_file})->pack(-side=>'left');
    my $run_blast_db_browse=$run_blast_db_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['blast db','*']);
	$runstate{blast_db_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select prepressed db sequence FASTA file");})->pack(-side=>'right');
    
    my $run_blast_db_seq_frame=$run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_db_seq_cb=$run_blast_db_seq_frame->Checkbutton(-text=>'Use sequence FASTA file',-variable=>\$runstate{blast_db_seq},-command=>sub {$runstate{blast_db}=0; $runstate{blast_db_remote}=0;})->pack(-side=>'left');
    my $run_blast_db_seq_entry=$run_blast_db_seq_frame->Entry(-textvariable=>\$runstate{blast_db_seq_file})->pack(-side=>'left');
    my $run_blast_db_seq_browse=$run_blast_db_seq_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['FASTA file','*']);
	$runstate{blast_db_seq_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Select db sequence FASTA file");
    })->pack(-side=>'right');
    my $run_blast_db_remote_frame=$run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_db_remote_cb=$run_blast_db_remote_frame->Checkbutton(-text=>'Remote blast db',-variable=>\$runstate{blast_db_remote}, -command=>sub {$runstate{blast_db_seq}=0; $runstate{blast_db}=0;})->pack(-side=>'left');
    my $run_blast_db_remote_entry=$run_blast_db_remote_frame->Entry(-textvariable=>\$runstate{blast_remote_db_name})->pack(-side=>'left');

    my $run_blast_pro_frame=$run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_pro_label=$run_blast_db_remote_frame->Label(-text=>'')->pack(-side=>'left');

    #$run_blast_pro_frame->  # drop down menu with blast programs...

    my $run_blastprog_label = $run_blast_page->Label(-text=>"Blast type")->pack(-side=>'top',-anchor=>'w');
    my $run_blastprog_frame = $run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blastprog_rb_x = $run_blastprog_frame->Radiobutton(-text=>'blastx',-variable=>\$runstate{blastx}, -value=>'1',-command=>sub { $runstate{blastn} = 0; })->pack(-side=>'left');
    my $run_blastprog_rb_n = $run_blastprog_frame->Radiobutton(-text=>'blastn',-variable=>\$runstate{blastn}, -value=>'1',-command=>sub { $runstate{blastx} = 0; })->pack(-side=>'left');
    
    # complexity filter 
    my $run_blastprogcompl_frame = $run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blastprogcompl_cb=$run_blastprogcompl_frame->Checkbutton(-text=>'Low complexity filter',-variable=>\$runstate{blast_complexityfilter})->pack(-side=>'left');

    # p threshold or use hitlist
    my $run_blast_scorefilter_frame = $run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_scorefilter_cb = $run_blast_scorefilter_frame->Checkbutton(-text=>'P value threshold',-variable=>\$runstate{blast_threshold})->pack(-side=>'left');
    my $run_blast_scorefilter_entry=$run_blast_scorefilter_frame->Entry(-textvariable=>\$runstate{blast_pvalue})->pack(-side=>'left');
    
    my $run_blast_hitlist_frame = $run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_hitlist_cb = $run_blast_hitlist_frame->Checkbutton(-text=>'Results in huge hitlist',-variable=>\$runstate{blast_hitlist})->pack(-side=>'left');

    # local package
    my $run_blast_blast_label = $run_blast_page->Label(-text=>"Local blast package")->pack(-side=>'top',-anchor=>'w');
    my $run_blast_blast_frame = $run_blast_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_blast_blast_rb_ncbi = $run_blast_blast_frame->Radiobutton(-text=>'NCBI',-variable=>\$runstate{blast_ncbi}, -value=>'1',-command=>sub { $runstate{blast_wu} = 0; })->pack(-side=>'left');
    my $run_blast_blast_rb_wu = $run_blast_blast_frame->Radiobutton(-text=>'WU-BLAST2',-variable=>\$runstate{blast_wu}, -value=>'1',-command=>sub { $runstate{blast_ncbi} = 0; })->pack(-side=>'left');
    
    my $run_testcode_page = $run_notebook->add("Testcode", -label => "Testcode", -underline=>0);
    
    my $run_testcode_label=$run_testcode_page->Label(-text=>"Run Testcode on")->pack(-anchor=>'w',-side=>'top',-expand=>'y',-fill=>'x');

    my $r_testcode_run_testcode_frame=$run_testcode_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_testcode_run_seq_cb=$r_testcode_run_testcode_frame->Checkbutton(-text=>"Current sequence $$seqName",-variable=>\$runstate{testcode_current_seq},-command=>sub{ $runstate{testcode_run}=0; })->pack(-side=>'left');
    my $r_testcode_run_frame=$run_testcode_page->Frame()->pack(-side=>'top',-expand=>'y',-fill=>'x');
    my $r_testcode_run_cb=$r_testcode_run_frame->Checkbutton(-text=>'File',-variable=>\$runstate{testcode_run},-command=>sub{ $runstate{testcode_current_seq}=0; })->pack(-side=>'left');
    my $r_testcode_run_entry=$r_testcode_run_frame->Entry(-textvariable=>\$runstate{testcode_seq_file})->pack(-side=>'left');
    my $r_testcode_run_browse=$r_testcode_run_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['Sequence FASTA file','*']);
	$runstate{testcode_run_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Run testcode on");
    })->pack(-side=>'right');
    
    my $r_testcode_options_label=$run_testcode_page->Label(-text=>"Options")->pack(-anchor=>'w',-side=>'top',-expand=>'y',-fill=>'x');

    my $r_testcode_mwinl_frame=$run_testcode_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $r_testcode_mwinl_label=$r_testcode_mwinl_frame->Label(-text=>"Mwinl")->pack(-side=>'left');
    my $r_testcode_mwinl_entry=$r_testcode_mwinl_frame->Entry(-textvariable=>\$runstate{testcode_mwinl})->pack(-side=>'left');

    my $r_testcode_increment_frame=$run_testcode_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $r_testcode_increment_label=$r_testcode_increment_frame->Label(-text=>"Increment")->pack(-side=>'left');
    my $r_testcode_increment_entry=$r_testcode_increment_frame->Entry(-textvariable=>\$runstate{testcode_increment})->pack(-side=>'left');

#    my $r_testcode_shortest_orf_frame=$run_testcode_page->Frame()->pack(-side=>'top',-fill=>'x');
#    my $r_testcode_shortest_orf_label=$r_testcode_shortest_orf_frame->Label(-text=>"Shortest ORF")->pack(-side=>'left');
#    my $r_testcode_shortest_orf_entry=$r_testcode_shortest_orf_frame->Entry(-textvariable=>\$runstate{testcode_shortest_orf})->pack(-side=>'left');

    # splice model

    my $run_sm_page=$run_notebook->add("SM", -label => "Splicemodel", -underline=>0);
    my $run_sm_label=$run_sm_page->Label(-text=>"Splice model")->pack(-side=>'top',-anchor=>'w');

     my $run_sm_currentseq_rb=$run_sm_page->Checkbutton(-text=>"Current sequence $$seqName",-variable=>\$runstate{sm_currentseq})->pack(-side=>'top',-anchor=>'w'); # -fill=>'x' perhaps?
     
    my $run_sm_report_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_report_cb=$run_sm_report_frame->Checkbutton(-text=>'Create report file',-variable=>\$runstate{sm_report})->pack(-side=>'left');
    my $run_sm_report_entry=$run_sm_report_frame->Entry(-textvariable=>\$runstate{sm_report_file})->pack(-side=>'left');
    my $run_sm_report_browse=$run_sm_report_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['splicemodel report','*']);
	$runstate{sm_report_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Create Splicemodel report file");})->pack(-side=>'right'); 

    my $run_sm_stats_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_stats_cb=$run_sm_stats_frame->Checkbutton(-text=>'Create statistics file',-variable=>\$runstate{sm_stats})->pack(-side=>'left');
    my $run_sm_stats_entry=$run_sm_stats_frame->Entry(-textvariable=>\$runstate{sm_stats_file})->pack(-side=>'left');
    my $run_sm_stats_browse=$run_sm_stats_frame->Button(-text=>'Browse...',-command=>sub {
	my @filetypes=(['R file','*']);
	$runstate{sm_stats_file}=$run_win->getOpenFile(-filetypes=>\@filetypes,-initialfile => "",-defaultextension => '',-title=>"Create Splicemodel R statistics file");})->pack(-side=>'right');

    my $run_sm_uorf_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_uorf_cb=$run_sm_uorf_frame->Checkbutton(-text=>'Predict uORFs',-variable=>\$runstate{sm_uorf})->pack(-side=>'left');


    # short/long (radio w/ alternative splicing as default isn't half bad)
    my $run_sm_short_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_short_cb=$run_sm_short_frame->Checkbutton(-text=>'Only shortest possible UTRs',-variable=>\$runstate{sm_short})->pack(-side=>'left');

    my $run_sm_long_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_long_cb=$run_sm_long_frame->Checkbutton(-text=>'Only longest possible UTRs',-variable=>\$runstate{sm_long})->pack(-side=>'left');

    my $run_sm_polypyag_frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
    my $run_sm_polypyag_label=$run_sm_polypyag_frame->Label(-text=>'Maximum polyPy - AG distance')->pack(-side=>'left');
    my $run_sm_polypyag_entry=$run_sm_polypyag_frame->Entry(-textvariable=>\$runstate{sm_polypyag})->pack(-side=>'left');

#      my $run_sm__frame=$run_sm_page->Frame()->pack(-side=>'top',-fill=>'x');
#      my $run_sm_uorf_cb=$run_sm_uorf_frame->Checkbutton(-text=>'Predict uORFs',-variable=>\$runstate{sm_uorf})->pack(-side=>'left');

    # pEST

    my $run_est_page=$run_notebook->add("EST", -label => "EST processing", -underline=>0);
    my $run_est_label=$run_est_page->Label(-text=>"EST Processing.\nTo be implemented upon request.\nPlease invoke the separate command line programs for this purpouse.")->pack(-side=>'top',-anchor=>'w');

    my $run_action_frame=$run_main_frame->Frame()->pack(-expand=>'yes',-fill=>'both',-side=>'top');

    my $run_ok_button=$run_action_frame->Button(-text=>'Ok',-command=> sub { run_external($main,$canvas,$status_text,$run_win,\%runstate,$annotation,$seqName,$seq,$sheet);
									     $run_win->destroy();})->pack(-side=>'left',-padx=>6,-pady=>6);
    my $run_cancel_button=$run_action_frame->Button(-text=>'Cancel',-command=> sub { $run_win->destroy();})->pack(-padx=>6,-pady=>6,-side=>'left');

    if($first_tab eq 'glimmer2') {
	$run_notebook->raise("glimmer");
    } elsif($first_tab eq 'estorf') {
	$run_notebook->raise("ESTORF");
    } elsif($first_tab eq 'testcode') {
	$run_notebook->raise("Testcode");
    } elsif($first_tab eq 'pest') {
	$run_notebook->raise("EST");
    } elsif($first_tab eq 'blast') {
	$run_notebook->raise("blast");
    } elsif($first_tab eq 'sm') {
	$run_notebook->raise("sm");
    }
}

### Graphic routines

sub plot_axis {
    my $canvas=shift;
    my $begin=shift;
    my $end=shift;
    my $sheet=shift;

    my $zoomfactor = $$sheet{zoomfactor};

    #    my $canvasheight=$$canvas->cget("height");
    my $yoffset=$$sheet{canvas_seq_height};
    my $yscale=100; # y are typically frequencies 0 < x < 1.
    
    # Create axis..
    $$canvas->createLine($begin/$zoomfactor,$yscale/4+$yoffset,$end/$zoomfactor,$yscale/4+$yoffset,-tags=>(['plotAxis','graph']),-stipple=>"gray50",-fill=>"gray");
    $$canvas->createLine($begin/$zoomfactor,3*$yscale/4+$yoffset,$end/$zoomfactor,3*$yscale/4+$yoffset,-tags=>(['plotAxis','graph']),-stipple=>"gray50",-fill=>"gray");
    $$canvas->createLine($begin/$zoomfactor,$yoffset,$begin/$zoomfactor,$yscale+$yoffset,-tags=>(['plotAxis','graph']));
    $$canvas->createLine($end/$zoomfactor,$yoffset,$end/$zoomfactor,$yscale+$yoffset,-tags=>(['plotAxis','graph']));
    $$canvas->createLine($begin/$zoomfactor,$yoffset,$end/$zoomfactor,$yoffset,-tags=>(['plotAxis','graph']));
    $$canvas->createLine($begin/$zoomfactor,$yscale+$yoffset,$end/$zoomfactor,$yscale+$yoffset,-tags=>(['plotAxis','graph']));
    $$canvas->createLine($begin/$zoomfactor,$yscale/2+$yoffset,$end/$zoomfactor,$yscale/2+$yoffset,-tags=>(['plotAxis','graph']),-stipple=>"gray50",-fill=>"black");

    # Then show it! 
    my $newHeight=$$sheet{canvas_seq_height}+$yscale;

    $DEVEL && print "DEVEL: plotting axis with zoomfactor $zoomfactor, yoffset $$sheet{canvas_seq_height}, yscale $yscale, start ",$begin/$zoomfactor," stop ", $end/$zoomfactor,".\n";
    $$canvas->configure(-height=>$newHeight);
}

sub plot_curve {
    my $canvas=shift;
    my $points=shift;
#    my @x=@_[0..$points-1];
#    my @y=@_[$points..2*$points-1];
    my $x=shift;
    my $y=shift;
    my $plotname=shift;
    my $color=shift;
    my $sheet=shift;

    my $zoomfactor = $$sheet{zoomfactor};

    my $yoffset=$$sheet{canvas_seq_height}; # Globals are sometimes nice.. If this program goes object oriented, this would be accessed by a display method..
    my $yscale=100; # y are typically frequencies 0 < x < 1. (Or normalised to this.. =)
    
    # determine max(@y) and scale accordingly.
    # Errhm. Might not be such a groovy thing to scale the frequencies. They just might look a bit odd.
    # my $ymax=max(@$y);
    # print "DEBUG: The maximum of $plotname is $ymax. Scaling data accordingly..\n";
    #    for($i=0;$i<$points;$i++) {
    #      $$y[$i]=$$y[$i]/$ymax;
    #    }

    $DEBUG && print "DEBUG: plotting curve with $points points (x is ".scalar(@$x)." and y ".scalar(@$y)." named $plotname, color $color at zoomf $zoomfactor.\n";

    for(my $i=0;$i<$points-1;$i++) {
	$$canvas->createLine(${$x}[$i]/$zoomfactor,(1-${$y}[$i])*$yscale+$yoffset,${$x}[$i+1]/$zoomfactor,(1-${$y}[$i+1])*$yscale+$yoffset,-tags=>(['graph',$plotname]),-fill=>$color);
#      print "DEBUG: plot (x,y)=($$x[$i],$$y[$i])\n";
    }
}

sub set_main_view_height {
    my $sheet = shift;

    my @annotators = keys %annotatorOrder;
    @annotators = sort { $annotatorOrder{$a} <=> $annotatorOrder{$b} } @annotators;

    %annotatorDisplayLevel = (); # global
    $mainViewHeight = 1; # global
    foreach my $annotator (@annotators) {
	if($annotator eq "ST" && $$sheet{display_ORF_level} == 1) {	
	    # these occupy same levels as ORF
	} else {
	    if($$sheet{"display_${annotator}_level"} == 1) {
		$annotatorDisplayLevel{$annotator} = $mainViewHeight;
		$mainViewHeight += $annotatorHeight{$annotator};
	    }
	}
    }
    $annotatorDisplayLevel{ST}=$annotatorDisplayLevel{ORF};
    
    if($mainViewHeight > $mainViewMaxHeight) {
	$WARNING && print "WARNING: with current level display selection, main view height is $mainViewHeight, which is larger than the maximum $mainViewMaxHeight. Continuing anyway...\n";
    }

    $$sheet{canvas_seq_height} = 40 + 20 * $mainViewHeight;

    $DEVEL && print "DEVEL: main view is $mainViewHeight units high, so canvas seq height is $$sheet{canvas_seq_height}.\n";
}

sub draw_seq {
    my $canvas=shift;
    my $seqName=shift;
    my $seq=shift; 
    my $sheet=shift;

    $DEBUG && ((ref $seq) || die "Seq is not a ref in draw_seq!\n");

    my $length=length($$seq);

    # set main view height
    set_main_view_height($sheet);

    my $y=$$sheet{canvas_seq_height}/2;

    my $tickHPos;
    my $tickHNeg;
    my $numberingGap;
    my $lineColor='black';
    
    my $tickmark_spacing;
    my $seq_block_len;
    my $zoomfactor;

    my($x1,$x2);
    
    # delete old
    $$canvas->delete(tags=>'sequence');
    $$canvas->delete(tags=>'tickmark');
    $$canvas->delete(tags=>'guide');

    if($$sheet{zoom} eq "normal") {
	$zoomfactor=1;
	$$sheet{zoomfactor}=$zoomfactor;
	$tickHPos=3;
	$tickHNeg=3;
	$numberingGap=10*$mainViewHeight+$tickHPos;
	$tickmark_spacing=50;
	($x1,$x2)=(1,$length);
    } elsif ($$sheet{zoom} eq "overview") {
	$zoomfactor=5;
	$$sheet{zoomfactor}=$zoomfactor;
	$tickHPos=3;
	$tickHNeg=3;
	$numberingGap=10*$mainViewHeight+$tickHPos;
	$tickmark_spacing=250;
	($x1,$x2)=(1,$length);
    } elsif ($$sheet{zoom} eq "birdseye") {
	$zoomfactor=40;
	$$sheet{zoomfactor}=$zoomfactor;
#      $zoomfactor=120; # just for a picture... =)
	$tickHPos=3;
	$tickHNeg=3;
	$numberingGap=10*$mainViewHeight+$tickHPos;
#    $tickmark_spacing=2000;
	$tickmark_spacing=6000;
	($x1,$x2)=(1,$length);
    } elsif ($$sheet{zoom} eq "factor") {
	$zoomfactor=$$sheet{zoomfactor};
	$tickHPos=3;
	$tickHNeg=3;
	$numberingGap=10*$mainViewHeight+$tickHPos;
	$tickmark_spacing= 50 * POSIX::ceil($zoomfactor);
	($x1,$x2)=(1,$length);
    } elsif ($$sheet{zoom} eq "sequence") {
#    print "DEBUG: ",$$canvas->itemconfigure('sequence'),"\n";
	my $seqlength = $$canvas->fontMeasure('fixed',$$seq);
	my $seqascent = $$canvas->fontMetrics('fixed',-ascent);
	my $seqdescent = $$canvas->fontMetrics('fixed',-descent);
	$DEBUG && print "DEBUG: Accordning to fontMeasure, the sequence is $seqlength pixels long, ascent $seqascent and descent $seqdescent.\n";
	$zoomfactor=($length)/$seqlength;
	$$sheet{zoomfactor}=$zoomfactor;
	$numberingGap=10*$mainViewHeight+3;
	$tickmark_spacing=25;
	$seq_block_len = 5000;
	($x1,$x2)=(1,$length+1);
	$DEBUG && print "DEBUG: So, zoomfactor is $zoomfactor, y is $y, x1 is $x1, x2 is $x2, length is $length, tickmarkspacing is $tickmark_spacing, seq_block_len is $seq_block_len and numberingGap is $numberingGap.\n";

	# Draw faint guide-boxes
	for(my $i = 1; $i < $mainViewHeight; $i+=2) {
	    my $guide_y = $$sheet{canvas_seq_height}/2+$i*10+5;
	    my $box = $$canvas->createLine(1/$zoomfactor,$guide_y,($length+1)/$zoomfactor, $guide_y, -width=>10,-fill=>$$sheet{faint_guide_color}, -tags=>(['guide']));

	    # And on the reverse levels
	    $guide_y = $$sheet{canvas_seq_height}/2-$i*10-5;
	    $box = $$canvas->createLine(1/$zoomfactor,$guide_y,($length+1)/$zoomfactor, $guide_y, -width=>10,-fill=>$$sheet{faint_guide_color}, -tags=>(['guide']));
	}

	# Top & bottom rulers

	for(my $i=0;$i<$length;$i+=$tickmark_spacing) {
	    my $xseqpos = $x1+$i;
	    my $xpos = $xseqpos/$zoomfactor;
	    
	    $$canvas->createLine($xpos,$y+$numberingGap,$xpos,$y-$numberingGap,-tags=>'tickmark',-fill=>'black',-stipple=>'gray50');
	    $$canvas->createText($xpos,$y-$numberingGap,-text=>"$xseqpos",-anchor=>'s',-tags=>'tickmark');
	    $$canvas->createText($xpos,$y+$numberingGap,-text=>"$xseqpos",-anchor=>'n',-tags=>'tickmark');
	}

	$$canvas->createLine($x2/$zoomfactor,$y+$numberingGap,$x2/$zoomfactor,$y-$numberingGap,-fill=>'black',-stipple=>'gray50');
	$$canvas->createText($x2/$zoomfactor,$y-$numberingGap,-text=>"$length",-anchor=>'s',-tags=>'tickmark');
	$$canvas->createText($x2/$zoomfactor,$y+$numberingGap,-text=>"$length",-anchor=>'n',-tags=>'tickmark');

	$$canvas->configure(-scrollregion=>[1/$zoomfactor,1,($length+1)/$zoomfactor,$y]);
	# $$canvas->configure(-scrollregion=>[1-($tickmark_spacing)/$zoomfactor,1,($length+$tickmark_spacing)/$zoomfactor,$y]); 
	# Display sequence selection
	
	draw_seq_text($canvas,$seqName,$seq,$sheet);
	
	return;
    }

    # Draw faint guide-boxes 
    # $DEBUG && print "DEBUG: there are ",max(values(%annotatorLevel))+2," levels (2 + max out of @tempLevels to draw annotations on..\n"; # +2 since the last one has three levels..

    for(my $i = 1; $i < $mainViewHeight; $i+=2) {

	my $guide_y = $$sheet{canvas_seq_height}/2+$i*10+5;
	my $box = $$canvas->createLine(1/$zoomfactor,$guide_y,($length+1)/$zoomfactor, $guide_y, -width=>10,-fill=>$$sheet{faint_guide_color}, -tags=>(['guide'])); 

	# And on the reverse levels
	$guide_y = $$sheet{canvas_seq_height}/2-$i*10-5;
	$box = $$canvas->createLine(1/$zoomfactor,$guide_y,($length+1)/$zoomfactor, $guide_y, -width=>10,-fill=>$$sheet{faint_guide_color}, -tags=>(['guide']));          
    }

    $DEVEL && print "DEVEL: Placing tickmarked line on height $y\n";
    
    # Top & bottom rulers
    $$canvas->createLine($x1,$y-$numberingGap,$x2/$zoomfactor,$y-$numberingGap,-fill =>$lineColor,-tags=>'sequence');
    $$canvas->createLine($x1,$y+$numberingGap,$x2/$zoomfactor,$y+$numberingGap,-fill =>$lineColor,-tags=>'sequence');
    
    $$canvas->createLine($x1,$y,$x2/$zoomfactor,$y,-fill =>$lineColor,-tags =>(['sequence',$seqName]));

#  $$canvas->createLine($x1,$y-$tickHNeg,$x1,$y+$tickHPos,-fill =>$lineColor,-tags=>('tickmark'));
#  $$canvas->createText($x1,$y-$numberingGap,-text=>"$x1",-anchor=>'n');

    for(my $i=0;$i<$length-$tickmark_spacing;$i+=$tickmark_spacing) {
	my $xseqpos=$x1+$i;
	my $xpos=$xseqpos/$zoomfactor;

	# TEST: Extra rulers...
	$$canvas->createLine($xpos,$y-$numberingGap-$tickHNeg,$xpos,$y-$numberingGap+$tickHPos,-fill =>$lineColor,-tags=>'tickmark');
	$$canvas->createLine($xpos,$y+$numberingGap-$tickHNeg,$xpos,$y+$numberingGap+$tickHPos,-fill =>$lineColor,-tags=>'tickmark');
	$$canvas->createLine($xpos,$y+$numberingGap,$xpos,$y-$numberingGap,-tags=>'tickmark',-fill=>'black',-stipple=>'gray50');    

	$$canvas->createLine($xpos,$y-$tickHNeg,$xpos,$y+$tickHPos,-fill =>$lineColor,-tags=>'tickmark');
	$$canvas->createText($xpos,$y-$numberingGap,-text=>"$xseqpos",-anchor=>'s',-tags=>'tickmark');
	$$canvas->createText($xpos,$y+$numberingGap,-text=>"$xseqpos",-anchor=>'n',-tags=>'tickmark');
    }

    # TEST: Extra rulers...
    $$canvas->createLine($x2/$zoomfactor,$y-$numberingGap-$tickHNeg,$x2/$zoomfactor,$y-$numberingGap+$tickHPos,-fill =>$lineColor,-tags=>'tickmark');
    $$canvas->createLine($x2/$zoomfactor,$y+$numberingGap-$tickHNeg,$x2/$zoomfactor,$y+$numberingGap+$tickHPos,-fill =>$lineColor,-tags=>'tickmark');
    $$canvas->createLine($x2/$zoomfactor,$y+$numberingGap,$x2/$zoomfactor,$y-$numberingGap,-tags=>'tickmark',-fill=>'black',-stipple=>'gray50');    

    $$canvas->createLine($x2/$zoomfactor,$y-$tickHNeg,$x2/$zoomfactor,$y+$tickHPos,-fill =>$lineColor,-tags=>('tickmark'));
    $$canvas->createText($x2/$zoomfactor,$y-$numberingGap,-text=>"$length",-anchor=>'s',-tags=>'tickmark');
    $$canvas->createText($x2/$zoomfactor,$y+$numberingGap,-text=>"$length",-anchor=>'n',-tags=>'tickmark');
    
    $$canvas->configure(-scrollregion=>[1-$tickmark_spacing/$zoomfactor,1,($length+$tickmark_spacing)/$zoomfactor,$y]);
}

sub draw_seq_text {
    my $canvas=shift;
    my $seqName=shift;
    my $seq=shift; 
    my $sheet=shift;

    my $length=length($$seq);

    # count active levels..
    set_main_view_height($sheet);

    # my $y=$canvas->cget("height")/2;
    my $y=$$sheet{canvas_seq_height}/2;

    my $tickHPos;
    my $tickHNeg;
    my $numberingGap;
    my $lineColor='black';
    
    my $seq_block_len;
    my $zoomfactor;

    my($x1,$x2);
    
    # delete old
    $$canvas->delete(tags=>'sequence');
    
    if ($$sheet{zoom} eq "sequence") {
#    print "DEBUG: ",$$canvas->itemconfigure('sequence'),"\n";
	my $seqlength=$$canvas->fontMeasure('fixed',$$seq);
	my $seqascent=$$canvas->fontMetrics('fixed',-ascent);
	my $seqdescent=$$canvas->fontMetrics('fixed',-descent);
	$DEBUG && print "DEBUG: Accordning to fontMeasure, the sequence is $seqlength pixels long, ascent $seqascent and descent $seqdescent.\n";
	my $zoomfactor=$length/$seqlength;
	$$sheet{zoomfactor}=$zoomfactor;
	$numberingGap=10*$mainViewHeight+3;
	$seq_block_len = 5000;
	($x1,$x2)=(1,$length+1);
	$DEBUG && print "DEBUG: So, zoomfactor is $zoomfactor, y is $y, x1 is $x1, x2 is $x2, length is $length, seq_block_len is $seq_block_len and numberingGap is $numberingGap.\n";

	# Top & bottom rulers are regarded as sequence
	$$canvas->createLine($x1,$y-$numberingGap,$x2/$zoomfactor,$y-$numberingGap,-fill =>$lineColor,-tags=>(['sequence',$seqName]));
	$$canvas->createLine($x1,$y+$numberingGap,$x2/$zoomfactor,$y+$numberingGap,-fill =>$lineColor,-tags=>(['sequence',$seqName]));

	# Display sequence selection
	
	$DEVEL && print "DEVEL: entering sequence block display loop..\n";
	
	# seq_sel_start & seq_sel_stop are 1-based, $i is 0-based
	for(my $i=0;$i<$length;$i+=$seq_block_len) {
	    my $xseqpos = $x1+$i;
	    my $xpos = $xseqpos/$zoomfactor;	    

	    if( ( $i > ($$sheet{seq_sel_start}-1 - $seq_block_len ) and ($i < $$sheet{seq_sel_start}-1 ) ) ) { # within blocklength bases upstream of selection start?
		# deal with the leftmost, unselected part of this sequence block
		$DEVEL && print "DEVEL: deal with the leftmost, unselected part of this sequence block since seq_sel_start $$sheet{seq_sel_start} and i $i.\n";
		my $delta_l = $$sheet{seq_sel_start}-1 - $i + 1 - 1; # leftmost, unselected sequence length: i left of (<) selection start. Take away one since selection has not yet started.

		my $subseq = substr($$seq,$i,$delta_l);
		$$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i],-font=>'fixed',-text=>$subseq);

		$subseq=~tr/atgcATGC/tacgTACG/;
		$$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i, "comp"],-font=>'fixed',-text=>$subseq);
		$DEVEL && print "DEVEL: added seqence segment with first pos ",$i,".\n";

		# draw the first selection piece - to end of block or end of selection
		$xseqpos = $x1+$i+$delta_l; # update position
		$xpos = $xseqpos/$zoomfactor;

		my $sel_len = $$sheet{seq_sel_stop} - $$sheet{seq_sel_start} + 1;

		if($$sheet{seq_sel_stop} >= $i + $seq_block_len ) {
		    # selection stretches to end of block, or further
		    $DEVEL && print "DEVEL: selection stretches to end of block, or further since seq_sel_stop $$sheet{seq_sel_stop}, i $i and seq_block_len $seq_block_len.\n";
		    $subseq = substr($$seq,$i+$delta_l, $seq_block_len-1 - $delta_l +1); # can't use sel length because selection can be larger than blocksize (and stopping this left of the block end..)
		    $$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i+$delta_l)],-font=>'fixed',-fill=>'red',-text=>$subseq);

		    $subseq=~tr/atgcATGC/tacgTACG/;
		    $$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i+$delta_l), "comp"],-font=>'fixed',-fill=>'red',-text=>$subseq);
		    $DEVEL && print "DEVEL: added selected seqence segment with first pos ",($i+$delta_l),",and length ", $seq_block_len - $delta_l + 1,".\n";
		    
		} elsif($$sheet{seq_sel_stop} < $i + $seq_block_len ) { 
		    # selection ends within block
		    $DEVEL && print "DEVEL: selection ends within block, since seq_sel_stop $$sheet{seq_sel_stop}, i $i and seq_block_len is $seq_block_len.\n";
		    $subseq = substr($$seq,$i+$delta_l, $sel_len);
		    $$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i+$delta_l)],-font=>'fixed',-fill=>'red',-text=>$subseq);

		    $subseq=~tr/atgcATGC/tacgTACG/;
		    $$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i+$delta_l), "comp"],-font=>'fixed',-fill=>'red',-text=>$subseq);
		    $DEVEL && print "DEVEL: added selected seqence segment with first pos ",($i+$delta_l)," and length ",$sel_len ,".\n";

		    my $delta_r = $seq_block_len - ($$sheet{seq_sel_stop}-1 - $i) +1;

		    # update position
		    $xseqpos = $x1+$i+$delta_l+$sel_len; 
		    $xpos = $xseqpos/$zoomfactor;

		    $subseq = substr($$seq, $i + $delta_l + $sel_len ,$delta_r);
		    $$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i + $delta_l + $sel_len )],-font=>'fixed',-text=>$subseq);

		    $subseq=~tr/atgcATGC/tacgTACG/;
		    $$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".($i + $delta_l + $sel_len ), "comp"],-font=>'fixed',-text=>$subseq);
		    $DEVEL && print "DEVEL: added seqence segment with first pos ",($i + $delta_l + $sel_len)," and length ",$delta_r,".\n";
		}

	    } elsif ( $i >= $$sheet{seq_sel_start}-1 && ( ( $i + $seq_block_len ) <= $$sheet{seq_sel_stop}-1 ) )  { # entire block within selection
		my $subseq = substr($$seq,$xseqpos-1,$seq_block_len);
		$$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i],-font=>'fixed',-fill=>'red',-text=>$subseq);
		$subseq=~tr/atgcATGC/tacgTACG/;
		$$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i, "comp"],-font=>'fixed',-fill=>'red',-text=>$subseq);
		$DEVEL && print "DEVEL: added selected sequence segment with first pos ",$i,".\n";
	    } elsif ( $i >= $$sheet{seq_sel_start}-1 && $i < $$sheet{seq_sel_stop} ) { # marking ends within current block, but there are no unselected bases upstream in block
		$DEVEL && print "DEVEL: leftmost part of block selected.\n";
		# Draw rightmost selection piece
		my $delta_r = $seq_block_len - ($$sheet{seq_sel_stop}-1 - $i) +1;
		my $delta_s_r = $$sheet{seq_sel_stop}-1 - $i + 1;

		my $subseq = substr($$seq,$i,$delta_s_r);
		$$canvas->createText($xpos, $y-5, -anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i],-font=>'fixed',-fill=>'red',-text=>$subseq);

		$subseq=~tr/atgcATGC/tacgTACG/;
		$$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i, "comp"],-font=>'fixed',-fill=>'red',-text=>$subseq);
		$DEVEL && print "DEVEL: added selected seqence segment with first pos ",$i," and length ",$delta_s_r," (seq_sel_start $$sheet{seq_sel_start} to seq_sel_stop $$sheet{seq_sel_stop}).\n";
				
		# Deal with the rightmost, unselected part of this sequence block
		
		# update position
		$xseqpos = $x1+$i+$delta_s_r;
		$xpos = $xseqpos/$zoomfactor;

		$subseq = substr($$seq, ($$sheet{seq_sel_stop}-1)+1, $delta_r);
		$$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence',  $seqName, "xseqpos".($i+$delta_s_r)], -font=>'fixed', -text=>$subseq);
		
		$subseq=~tr/atgcATGC/tacgTACG/;
		$$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence',  $seqName, "xseqpos".($i+$delta_s_r), "comp"],-font=>'fixed',-text=>$subseq);
		$DEVEL && print "DEVEL: added seqence segment with first pos ",$i,".\n";

	    } else {
		my $subseq = substr($$seq,$xseqpos-1,$seq_block_len);
		$$canvas->createText($xpos,$y-5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i],-font=>'fixed',-text=>$subseq);
		$subseq=~tr/atgcATGC/tacgTACG/;
		$$canvas->createText($xpos,$y+5,-anchor=>'w',-tags=>['sequence', $seqName, "xseqpos".$i, "comp"],-font=>'fixed',-text=>$subseq);
		$DEVEL && print "DEVEL: added seqence segment with first pos ",$i,".\n";
	    }
	}
	
	$DEVEL && print "DEVEL: exiting sequence block display loop..\n";

	return;

    } elsif($$sheet{zoom} eq "normal") {
	$DEBUG && print "DEBUG: Draw seq text skipped since zoom is $$sheet{zoom}.\n.";
    } elsif ($$sheet{zoom} eq "overview") {
	$DEBUG && print "DEBUG: Draw seq text skipped since zoom is $$sheet{zoom}.\n.";
    } elsif ($$sheet{zoom} eq "birdseye") {
	$DEBUG && print "DEBUG: Draw seq text skipped since zoom is $$sheet{zoom}.\n.";
    } elsif ($$sheet{zoom} eq "factor")  {
	$DEBUG && print "DEBUG: Draw seq text skipped since zoom is $$sheet{zoom}.\n.";
    }
}

sub draw_annotation_st {
    my $canvas = shift;
    my $annotation = shift;  
    my $nr = shift;
    my $sheet = shift;

    my $start_ref = $$annotation{start}->[$nr];
    my $stop_ref = $$annotation{stop}->[$nr];

    my $annotator_type = $$annotation{type}->[$nr];
    if($$sheet{"display_${annotator_type}_level"} != 1) {
	$DEVEL && print "DEVEL: draw_annotation_st repressing display of $annotator_type levels..\n";
	return;
    }

    my $level=$$annotation{level}->[$nr];
    
    # $DEVEL && print "DEVEL DEBUG: the following annotators have display level entries in draw ann st: ",keys %annotatorDisplayLevel,"\n";
    # $DEVEL && print "DEVEL DEBUG: called draw st for type $annotator_type, drawing at base level $annotatorDisplayLevel{$annotator_type}.\n";

    if ($$annotation{frame}->[$nr] < 3) {
	$level += $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] < 6) {
	$level -= $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] == 6) {
	$level += $annotatorDisplayLevel{$annotator_type};
    }

    my $annotation_id = $$annotation{id}->[$nr];
    my $color = $$annotation{color}->[$nr];
    my $name = $$annotation{name}->[$nr];

    my $zoomfactor;

    if($$sheet{zoom} eq "normal") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "overview") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "birdseye") {
	# Start and stop codons not really useful at this zoomlevel.
	return 0;
    } elsif ($$sheet{zoom} eq "factor") {
	$zoomfactor=$$sheet{zoomfactor};
	if($zoomfactor > 30) {
	    $DEBUG && print "DEBUG: Zoomfactor > 30 - start/stop not shown!";
	    return 0;
	}
    } elsif ($$sheet{zoom} eq "sequence") {
	$zoomfactor=$$sheet{zoomfactor};
	# or get it by means of seq-length, fontname, etc..
    }

    my ($line, $y);

    foreach my $start (@$start_ref) {
	$color = 'Green1';
	if($level>0) {
	    $y = $$sheet{canvas_seq_height}/2-$level*10-5;
	    # Draw line just upstream of start.
	    $line = $$canvas->createLine($start/$zoomfactor, $y, ($start+3)/$zoomfactor, $y,-width=>10, -fill=>$color, -tags=>(['annotation',$annotator_type,$annotation_id]));
	} else {
	    $y = $$sheet{canvas_seq_height}/2-$level*10+5;
	    if ($y<0) {
		print "ERROR: can't draw annotation box at height $y (level $level)\n";
	    } 
	    # Draw line just downstream of reverse-start.
	    $line = $$canvas->createLine(($start+1)/$zoomfactor, $y, ($start-2)/$zoomfactor, $y,-width=>10,-fill=>$color,-tags=>(['annotation',$annotator_type,$annotation_id]));
	}
	# print "Draw $annotation_id ($annotator_type) as start at $start on $level (i e y=$y) in $color\n";
	
    }
    foreach my $stop (@$stop_ref) {
	$color = 'Black';
	if($level>0) {
	    $y = $$sheet{canvas_seq_height}/2-$level*10-5;
	    # Draw line just downstream of stop.
	    $line = $$canvas->createLine(($stop-2)/$zoomfactor, $y, ($stop+1)/$zoomfactor, $y,-width=>10, -fill=>$color, -tags=>(['annotation',$annotator_type,$annotation_id]));
	} else {
	    $y = $$sheet{canvas_seq_height}/2-$level*10+5;
	    if ($y<0) {
		print "ERROR: can't draw annotation box at height $y (level $level)\n";
	    } 
	    # Draw line just upstream of reverse-stop.
	    $line = $$canvas->createLine($stop/$zoomfactor, $y, ($stop+3)/$zoomfactor, $y,-width=>10,-fill=>$color,-tags=>(['annotation',$annotator_type,$annotation_id]));
	}
	# print "Draw $annotation_id ($annotator_type) as stop at $stop on $level (i e y=$y) in $color\n";   
   }
}

sub draw_annotation_arrow {
    my $canvas = shift;
    my $annotation = shift;
    my $nr = shift;
    my $sheet = shift;
    
    my $seq = shift; # ref -- tentative -- check calls.. 
    
    my $start=$$annotation{start}->[$nr];
    my $stop=$$annotation{stop}->[$nr];
    my $annotator_type=$$annotation{type}->[$nr];
    if($$sheet{"display_${annotator_type}_level"} != 1) {
	$DEVEL && print "DEVEL: draw_annotation_arrow repressing display of $annotator_type levels..\n";
	return 0;
    }
    my $level=$$annotation{level}->[$nr];

    if ($$annotation{frame}->[$nr] < 3) {
	$level += $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] < 6) {
	$level -= $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] == 6) {
	$level += $annotatorDisplayLevel{$annotator_type};
    }

    my $annotation_id = $$annotation{id}->[$nr];
    my $color = $$annotation{color}->[$nr];
    my $name = $$annotation{name}->[$nr];

    $WARNING && (ref $seq || die("ERROR: seq has stopped being a ref once again in draw_annotaion arrow... A bug has been ressurected..\n"));

    my $y; 
    my $box;

    my $zoomfactor;

    my $annotated_protein;

    my $stip = ''; # solid
    if($annotation_id eq $$sheet{selected}) {
	$stip = 'gray50';
    }

    if($$sheet{zoom} eq "normal") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "overview") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "factor") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "birdseye") {
	$zoomfactor=$$sheet{zoomfactor};
#    $zoomfactor=120; # For a bigger picture
    } elsif ($$sheet{zoom} eq "sequence") {
	$zoomfactor=$$sheet{zoomfactor};  # or get it by means of seq-length, fontname, etc..

	if($annotatorProtein{$annotator_type} == 1 || ($annotatorProtein{$annotator_type} == 2 && $annotatorDirected{$annotator_type} == 1)) { # should differentiate - excludes gff though.. ok?
	    $DEVEL && print STDERR "DEVEL: start is $start, stop $stop, seq ref $seq for annotation of type $annotator_type\n";
	    my $annotated_sequence = substr($$seq,$start-1,$stop-$start+1);
	    
	    if($level>0) {
		# Forward reading frame
		$annotated_protein=nt2aa($annotated_sequence); # We could remove M and + if we'd like..
	    } else {
		# Reverse reading frame
		$_=reverse(split(/ */,$annotated_sequence));
		tr/atgcATGC/tacgTACG/;
		# Translate, and then turn it around again to get it display-friendly...
		$annotated_protein=reverse(split(/ */,nt2aa($_))); # We could remove M and + if we'd like..
	    }
	    $annotated_protein=~s/(.{1})/ $1 /g; # Introduce spacing to match codon size..   
	}
    }
    
    if($level>0) {
	$y=$$sheet{canvas_seq_height}/2-$level*10-5;
	$box = $$canvas->createLine(($start)/$zoomfactor,$y,($stop+1)/$zoomfactor, $y, -width=>10,-fill=>$color, -stipple=>$stip, -arrow=>'last', -arrowshape=>[5,5,0],-tags=>(['annotation',$annotator_type,$annotation_id]));
	if($$sheet{zoom} eq "sequence") {
	    # Display protein translation
	    if( $annotatorProtein{$annotator_type} == 1 ) {
		$$canvas->createText(($start)/$zoomfactor,$y,-fill=>'white',-text=>"$annotated_protein",-font=>'fixed',-anchor=>'w',-tags=>(['annotation',$annotator_type,$annotation_id]));
	    }
	} elsif (defined($name) && $name ne "" ) {
	    my $namelen = $$canvas->fontMeasure('fixed',$name);
	    # $DEBUG && print "DEBUG: The name $name is $namelen pixels wide and the arrow is ".($stop - $start + 1)/$zoomfactor."pixels.\n";
	    if ($namelen > ($stop - $start + 1)/$zoomfactor ) {
		my $namelenfactor = $namelen / (($stop - $start + 1)/$zoomfactor);
		my $cut_name_len = POSIX::floor(length($name) / $namelenfactor) -1; 
		if($cut_name_len < 0) {
		    $cut_name_len = 0;
		}
		$name = substr($name, 0, $cut_name_len);
		$namelen = $$canvas->fontMeasure('fixed',$name);
		# $DEBUG && print "DEBUG: and is thus reduced to $name of length $cut_name_len which is only $namelen pixels wide.\n";
	    }
	    # $DEBUG && print "DEBUG: length of name is ", length($name),"\n";
	    $$canvas->createText( ($start + ($stop - $start)/2)/$zoomfactor, $y,-fill=>'white',-anchor=>'center',-text=>"$name",-font=>'fixed',-tags=>(['annotation',$annotator_type,$annotation_id]));
	}
    } else {
	$y=$$sheet{canvas_seq_height}/2-$level*10+5;
	if ($y<0) {
	    print "ERROR: can't draw annotation box at height $y (level $level)\n";
	}
	$box = $$canvas->createLine((1+$start-1)/$zoomfactor,$y,($stop+1)/$zoomfactor, $y, -width=>10,-fill=>$color, -stipple=>$stip, -arrow=>'first', -arrowshape=>[5,5,0],-tags=>(['annotation',$annotator_type,$annotation_id]));
	if($$sheet{zoom} eq "sequence") {
	    # Display protein translation
	    if($annotatorProtein{$annotator_type} == 1 || ($annotatorProtein{$annotator_type} == 2 && $annotatorDirected{$annotator_type} == 1)) {
		$$canvas->createText(($stop+1)/$zoomfactor ,$y,-fill=>'white',-anchor=>'e',-text=>"$annotated_protein",-font=>'fixed',-tags=>(['annotation',$annotator_type,$annotation_id]));
	    }
	} elsif (defined($name) && $name ne "" ) {
	    my $namelen = $$canvas->fontMeasure('fixed',$name);
	    # $DEBUG && print "DEBUG: The name $name is $namelen pixels wide and the arrow is ".($stop - $start + 1)/$zoomfactor."pixels.\n";
	    if ($namelen > ($stop - $start + 1)/$zoomfactor ) {
		my $namelenfactor = $namelen / (($stop - $start + 1)/$zoomfactor);
		my $cut_name_len = POSIX::floor(length($name) / $namelenfactor) -1;
		if($cut_name_len < 0) {
		    $cut_name_len = 0;
		}
		$name = substr($name, 0, $cut_name_len);
		$namelen = $$canvas->fontMeasure('fixed',$name);
		# $DEBUG && print "DEBUG: and is thus reduced to $name of length $cut_name_len which is only $namelen pixels wide.\n";
	    }

	    $$canvas->createText( ($stop - ($stop - $start) / 2 )/$zoomfactor,$y,-fill=>'white',-anchor=>'center',-text=>"$name",-font=>'fixed',-tags=>(['annotation',$annotator_type,$annotation_id])); 
	}
    }
    # $DEBUG && print "Draw $annotation_id ($annotator_type) from $start to $stop on $level (i e y=$y) in $color\n";
    
    return $box;
}

sub draw_annotation_box{
    my $canvas=shift;
    my $annotation = shift;  
    my $nr = shift;
    my $sheet = shift;

    my $start=$$annotation{start}->[$nr];
    my $stop=$$annotation{stop}->[$nr];
    my $annotator_type=$$annotation{type}->[$nr];
    if($$sheet{"display_${annotator_type}_level"} != 1) {
	$DEVEL && print "DEVEL: draw_annotation_box repressing display of $annotator_type levels..\n";
	return 0;
    }

    my $level=$$annotation{level}->[$nr];
    if ($$annotation{frame}->[$nr] < 3) {
	$level += $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] < 6) {
	$level -= $annotatorDisplayLevel{$annotator_type};
    } elsif ($$annotation{frame}->[$nr] == 6) {
	$level += $annotatorDisplayLevel{$annotator_type};
    }
    my $annotation_id = $$annotation{id}->[$nr];
    my $color = $$annotation{color}->[$nr];
    my $name = $$annotation{name}->[$nr];

    my $stip = ''; # solid
    if($annotation_id eq $$sheet{selected}) {
	$stip = 'gray50';
    }

    my $zoomfactor;

    if($$sheet{zoom} eq "normal") {
	$zoomfactor=1;
    } elsif ($$sheet{zoom} eq "birdseye") {
	$zoomfactor=40;
#    $zoomfactor=120; # for a bigger picture
    } elsif ($$sheet{zoom} eq "overview") {
	$zoomfactor=5;
    } elsif ($$sheet{zoom} eq "factor") {
	$zoomfactor=$$sheet{zoomfactor};
    } elsif ($$sheet{zoom} eq "sequence") {
	$zoomfactor=$$sheet{zoomfactor};
	# or get it by means of seq-length, fontname, etc..
    }
    
    my $y; 
    my $box;
    $y=$$sheet{canvas_seq_height}/2+$level*10+5;
    $box = $$canvas->createLine($start/$zoomfactor,$y,($stop+1)/$zoomfactor, $y, -width=>10,-fill=>$color, -stipple=>$stip, -tags=>(['annotation',$annotator_type,$annotation_id]));
# $DEBUG &&   print "Draw $annotation_id ($annotator_type) from $start to $stop on $level (i e y=$y) in $color\n"; 
    $y=$$sheet{canvas_seq_height}/2-$level*10-5;
    $box = $$canvas->createLine($start/$zoomfactor,$y,($stop+1)/$zoomfactor, $y, -width=>10,-fill=>$color, -stipple=>$stip, -tags=>(['annotation',$annotator_type,$annotation_id]));
#  $DEBUG && print "Draw box $annotation_id ($annotator_type) from $start to $stop on $level (i e y=$y) in $color\n";

    if($$sheet{zoom} ne "sequence" && defined($name) && $name ne "" ) {
	my $namelen = $$canvas->fontMeasure('fixed',$name);
	
	if ($namelen > ($stop - $start + 1)/$zoomfactor ) {
	    my $namelenfactor = $namelen / (($stop - $start + 1)/$zoomfactor);
	    my $cut_name_len = POSIX::floor(length($name) / $namelenfactor) -1;  
	    if($cut_name_len < 0) {
		$cut_name_len = 0;
	    }
	    $name = substr($name, 0, $cut_name_len);
	    $namelen = $$canvas->fontMeasure('fixed',$name);
	    # $DEBUG && print "DEBUG: and is thus reduced to $name of length $cut_name_len which is only $namelen pixels wide (box is ".($stop - $start + 1)/$zoomfactor."pixels).\n";
	}
	$y=$$sheet{canvas_seq_height}/2+$level*10+5;
	$$canvas->createText( ($start + ($stop - $start)/2)/$zoomfactor, $y,-fill=>'white',-anchor=>'center',-text=>"$name",-font=>'fixed',-tags=>(['annotation',$annotator_type,$annotation_id]));
	$y=$$sheet{canvas_seq_height}/2-$level*10-5;
	$$canvas->createText( ($start + ($stop - $start)/2)/$zoomfactor, $y,-fill=>'white',-anchor=>'center',-text=>"$name",-font=>'fixed',-tags=>(['annotation',$annotator_type,$annotation_id]));
    }

    return $box;
}

sub view_refresh {
    my $main = shift;
    my $canvas = shift;
    my $annotation = shift;
    my $seqName = shift;
    my $seq = shift;
    my $plotstatus = shift;
    my $sheet = shift;

    draw_seq($canvas,$$seqName,$seq,$sheet);
    $$canvas->configure(-height=>$sheet{canvas_seq_height});
    redraw_annotations($canvas,$annotation,$sheet,$seq);
    if($plotstatus{axis} && ($sheet{zoom} ne 'sequence')) {
	$DEBUG && print "DEBUG: reintroducing graph...\n";
	$plotstatus{axis}=0; # Tell calculate_graph to redraw axis...
	replot($canvas,$seq,$plotstatus,$sheet);
	# calculate_graph($main,$canvas,$seq,$plotstatus,$sheet);
    }
}
