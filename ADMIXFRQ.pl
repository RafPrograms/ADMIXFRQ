#!/usr/bin/perl
#
#    The ADMIXFRQ Perl script designed to specify genetic marker allele frequencies
#    in admixed samples as described in Nafikov et al. (2018).
#    Copyright (C) 2018  Rafael A. Nafikov email <nrafscience@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

use strict;
use warnings;
use diagnostics;
########### controls for main subroutine ##################################################

# $global_adm_comp = 0 or 1 (0 = use supplied global admixture proportions, 1 = compute global admixture proportions from data)

# $adm_printout = 0 or 1 (0 = do not print out or 1 = print out admixture proportions)

# $local_adm_extrapolation = 0 or 1 (0 = extrapolate local admixture estimates beyond bundaries within which they were estimated or 1 = do not extrapolate)

# $fbcw_adm_where = 0 or 1 (0 = use $fbcw_adm_ref or 1 = use $fbcw_adm_b_ref for adm prinout)

# $control_g = 0 (nothing) or 1 for Global admixture model
# $control_fbgw = 0 (nothing) or 1 for FBGW admixture model
# $control_fbcw = 0 (nothing) or 1 for FBCW admixture model
# $control_l = 0 (nothing) or 1 for Local admixture model

# $control_1 = 0 (nothing) or 1 for computing allele frequencies for geno file
# $control_2 = 0 (nothing) or 1 for computing allele frequencies for user specified base pair positions

# $control_geno_extra = 0 (nothing) or 1 for printing out geno file frequencies into separate file

# $comp_type = 0 if frequencies for geno files are needed; 1 = if frequencies not for geno files are needed (FBGW, FBCW, Local); 2 = if admixture proportions for Local model are needed; 3 = if admixture proportions for Global, FBGW, or FBCW model are needed; 4 = directory for error log file; 5 = if frequencies not for geno files are needed (Global)
##########################################################################################


	main ();

sub main {

	$| = 1; # output autoflush
	my $parameters_ref = read_par_file ();
	my %parameters = %$parameters_ref;
	script_info ($parameters_ref);


	my @chrs; # retrieving information about chromosome numbers
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'CHROMOSOMES'}) {
	for my $key2 (sort {$a <=> $b} keys %{$parameters{$key1}{'CHROMOSOMES'}}) { 
	if ($parameters{$key1}{'CHROMOSOMES'}{$key2} =~ m/\.\./) {
	my @temp1 = split (/\.\./, $parameters{$key1}{'CHROMOSOMES'}{$key2});
	my @temp2 = ($temp1[0]..$temp1[1]);
	push @chrs, @temp2;
	} else {
	my @temp1 = split (/\.\./, $parameters{$key1}{'CHROMOSOMES'}{$key2});
	push @chrs, @temp1;
	}
	}
	}
	}
	error ("Error: chromosome numbers were not specified!") if (! @chrs);

	my @families; # retrieving information about family names
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'FAMILYNAMES'}) {
	for my $key2 (sort {$a <=> $b} keys %{$parameters{$key1}{'FAMILYNAMES'}}) { 
	push @families, $parameters{$key1}{'FAMILYNAMES'}{$key2};
	}
	}
	}
	error ("Error: family names were not provided!") if (! @families);

	my @models; # retrieving information about admixture models to use
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'ADMIXTUREMODEL'}) {
	for my $key2 (sort {$a <=> $b} keys %{$parameters{$key1}{'ADMIXTUREMODEL'}}) { 
	push @models, $parameters{$key1}{'ADMIXTUREMODEL'}{$key2};
	}
	}
	}
	error ("Error: admixture models were not specified!") if (! @models);

	my $global_adm_comp = 0; # 0 = use provided global admixture proportions; 1 = computed global admixture proportions from data
	for my $key1 (sort keys %parameters) {
	for my $key2 (sort keys %{$parameters{$key1}}) {
	if ($key2 eq 'COMPUTEGLOBALADMIXTUREMODELPROPORTIONS') {
	if (uc $parameters{$key1}{'COMPUTEGLOBALADMIXTUREMODELPROPORTIONS'} eq 'YES') {
	$global_adm_comp = 1;
	}
	}
	}
	}
	my $global_adm_comp_ref = \$global_adm_comp;

	my $adm_printout = 0; # 0 = do not export admixture proportions; 1 = export admixture proportions
	for my $key1 (sort keys %parameters) {
	for my $key2 (sort keys %{$parameters{$key1}}) {
	if ($key2 eq 'ADMIXTUREPROPORTIONSPRINTOUT') {
	if (uc $parameters{$key1}{'ADMIXTUREPROPORTIONSPRINTOUT'} eq 'YES') {
	$adm_printout = 1;
	}
	}
	}
	}
	my $adm_printout_ref = \$adm_printout;

	my $local_adm_extrapolation = 0; # 0 = extrapolate local admixture estimates beyond bundaries within which they were estimated or 1 = do not extrapolate
	for my $key1 (sort keys %parameters) {
	for my $key2 (sort keys %{$parameters{$key1}}) {
	if ($key2 eq 'EXTRAPOLATIONOFLOCALADMIXTUREPROPORTIONS') {
	if (uc $parameters{$key1}{'EXTRAPOLATIONOFLOCALADMIXTUREPROPORTIONS'} eq 'NO') {
	$local_adm_extrapolation = 1;
	}
	}
	}
	}
	my $local_adm_extrapolation_ref = \$local_adm_extrapolation;

	my $control_g = 0;
	my $control_fbgw = 0;
	my $control_fbcw = 0;
	my $control_l = 0;

	for my $el (@models) {
	$el = uc $el;
	if ($el eq 'GLOBAL') {
	$control_g = 1;
	}
	elsif ($el eq 'FBGW') {
	$control_fbgw = 1;
	}
	elsif ($el eq 'FBCW') {
	$control_fbcw = 1;
	}
	elsif ($el eq 'LOCAL') {
	$control_l = 1;
	}
	}

	my $control_1 = 0; # option for geno file related computations (0 = nothing or 1 = yes)
	my $control_geno_extra = 0; # option to output geno file frequencies into a separate file (0 = nothing or 1 = yes)
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'GENOFILEFORGL_AUTO'}) {
	if (exists $parameters{$key1}{'FILEPATH'}{1}) {
	$control_1 = 1;
	} else {
	error ("Error: location of geno files was not provided! Please specify the location or remove 'GENO FILE FOR GL_AUTO' option from parameter file.");
	}
	if (uc $parameters{$key1}{'DETAILEDALLELEFREQUENCYFILE'}{1} eq "YES") {
	$control_geno_extra = 1;
	}
	}
	}

	my $control_2 = 0; # option for computing allele frequencies for user specified base pair positions (0 = nothing or 1 = yes)
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'MODIFYYOURBPALLELEFREQUENCIES'}) {
	if ($parameters{$key1}{'FILEPATH'}{1}) {
	$control_2 = 1;
	} else {
	error ("Error: location of file with base pair positions for which allele frequencies were to be specified has not been provided! Please specify the location or remove 'MODIFY YOUR BP ALLELE FREQUENCIES' option from parameter file.");
	}
	}
	}


	my $fbgw_adm_ref;
	my $global_adm_ref; 
	my $fbcw_adm_b_ref;
	my $fbcw_adm_where = 0;


	if ($control_fbgw == 1 or $control_g = 1) { # computation of Global and FBGW admixtures using local admixture estimates for the entire dataset

	my $average_type = 0;
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'AVERAGETYPEFORGLOBALANDFBGW'}) {

	if (uc $parameters{$key1}{'AVERAGETYPEFORGLOBALANDFBGW'}{1} eq 'WEIGHTED') {
	$average_type = 1;
	}

	elsif (uc $parameters{$key1}{'AVERAGETYPEFORGLOBALANDFBGW'}{1} eq 'NORMAL') {
	$average_type = 2;
	}

	}
	}

	my $pop_adms_genome_ref;
	my $average_type_ref = \$average_type;
	($pop_adms_genome_ref, my $var_identifier_genome_adm_ref) = read_file_admixtures_genome ($parameters_ref);

	($fbgw_adm_ref, $global_adm_ref) = adms_global_fbgw ($average_type_ref, $pop_adms_genome_ref, $global_adm_comp_ref) if ($global_adm_comp == 1);

	$global_adm_ref = adm_global_supplied ($parameters_ref) if ($global_adm_comp == 0);

	$fbgw_adm_ref = adms_global_fbgw ($average_type_ref, $pop_adms_genome_ref, $global_adm_comp_ref) if ($global_adm_comp == 0);

	print_fbgw_adm_file ($fbgw_adm_ref, $parameters_ref) if ($control_fbgw == 1 && $adm_printout == 1);

	print_global_adm_file ($global_adm_ref, $parameters_ref) if ($control_g == 1 && $adm_printout == 1 && $global_adm_comp == 1);

	$fbcw_adm_b_ref = adms_fbcw_b ($pop_adms_genome_ref) if ($adm_printout == 1 && $control_fbcw == 1); 
	print_fbcw_adm_file ($fbcw_adm_b_ref, $parameters_ref, 0, 0) if ($adm_printout == 1 && $control_fbcw == 1);
	$fbcw_adm_where = 1 if ($adm_printout == 1 && $control_fbcw == 1);

	} # end of the loop for the computation of the Global and FBGW admixtures



	for my $chr (@chrs) { # chr loop
	my $chr_ref = \$chr;

###
	my $pop_adms_ref;
	my ($pop_frqs_ref, $rs_to_bp_ref) = read_file_population ($parameters_ref, $chr_ref);
	($pop_adms_ref, my $var_identifier_chr_adm) = read_file_admixtures_chr ($parameters_ref, $chr_ref);

	if (uc $$var_identifier_chr_adm eq "RS_NUMBER") {
	$pop_adms_ref = get_bp_for_adm_from_ref_pop ($pop_adms_ref, $rs_to_bp_ref);
	}

	my $sample_frqs_ref = read_file_sample_allele_frqs ($parameters_ref, $chr_ref) if ($control_2 == 1);
###


	my ($bp_pos_for_frq_ref, $frq_template_ref, $subpop_frq_template_ref, $global_frqs_ref, $sample_bp_rs_ref);
	if ($control_2 == 1) { # loop for modifying variant allele frequencies at user samplied bp positions - part 1
	($bp_pos_for_frq_ref, $sample_bp_rs_ref) = read_bp_modify_frq ($parameters_ref, $chr_ref) if ($control_2 == 1);
	($frq_template_ref, $subpop_frq_template_ref) = frq_template ($parameters_ref, $sample_frqs_ref, $pop_frqs_ref, $bp_pos_for_frq_ref, $chr_ref, $sample_bp_rs_ref) if ($control_2 == 1);

	$global_frqs_ref = global_adm_frqs ($frq_template_ref, $chr_ref, $subpop_frq_template_ref, $global_adm_ref) if ($control_g == 1 && $control_2 == 1);
	print_frqs ($parameters_ref, $global_frqs_ref, $chr_ref, 0, 1) if ($control_g == 1 && $control_2 == 1);

	} # end of the loop for modifying variant ellele frequencies at user samplied bp positions - part 1


	for my $family (@families) { # family loop
	my $family_ref = \$family;

### computation and printout of FBCW admixtures
	my $fbcw_adms_ref = adms_fbcw ($pop_adms_ref) if ($control_fbcw == 1);
	remove_old_fbcw_adm_file ($parameters_ref, $chr_ref, $family_ref) if ($adm_printout == 1 && $control_fbcw == 1 && $fbcw_adm_where == 0);
	print_fbcw_adm_file ($fbcw_adms_ref, $parameters_ref, $chr_ref, $family_ref) if ($adm_printout == 1 && $control_fbcw == 1 && $fbcw_adm_where == 0);
###


	if ($control_2 == 1) { # loop for modifying variant allele frequencies at user samplied bp positions - part 2

### FBGW admixture model
	my $fbgw_frqs_ref = fbgw_frqs ($frq_template_ref, $subpop_frq_template_ref, $fbgw_adm_ref,  $family_ref, $chr_ref) if ($control_fbgw == 1);
	print_frqs ($parameters_ref, $fbgw_frqs_ref, $chr_ref, $family_ref, 2) if ($control_fbgw == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local

### FBCW admixture model
	my $fbcw_frqs_ref = fbcw_frqs ($frq_template_ref, $subpop_frq_template_ref, $fbcw_adms_ref, $family_ref, $chr_ref) if ($control_fbcw == 1);
	print_frqs ($parameters_ref, $fbcw_frqs_ref, $chr_ref, $family_ref, 3) if ($control_fbcw == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local

### Local admixture model
	my ($local_adms_ref, $subpop_local_ref) = adms_local ($parameters_ref, $pop_adms_ref, $bp_pos_for_frq_ref, $chr_ref, $family_ref, $local_adm_extrapolation_ref) if ($control_l == 1);
	print_local_adm_file ($local_adms_ref, $subpop_local_ref, $parameters_ref, $chr_ref, $family_ref) if ($adm_printout == 1);
	my $local_frqs_ref = local_frqs ($frq_template_ref, $subpop_frq_template_ref, $local_adms_ref,  $family_ref, $chr_ref, $subpop_local_ref, 0) if ($control_l == 1);
	print_frqs ($parameters_ref, $local_frqs_ref, $chr_ref, $family_ref, 4) if ($control_l == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local

	} # end of the loop for modifying variant ellele frequencies at user samplied bp positions part 2


	if (($control_1 == 1)) { # loop for geno file allele frequencies
	my ($marker_cm_pos_ref, $geno_bp_rs_ref, $geno_genotypes_ref, $geno_line_1_ref, $geno_bps_ref) = read_geno_file_info ($parameters_ref, $chr_ref, $family_ref);

	my ($frq_template_geno_ref, $subpop_frq_template_geno_ref) = frq_template_geno ($parameters_ref, $pop_frqs_ref, $rs_to_bp_ref, $geno_bp_rs_ref, $chr_ref);

### Global admixture model
	my $global_frqs_geno_ref = global_adm_frqs ($frq_template_geno_ref, $chr_ref, $subpop_frq_template_geno_ref, $global_adm_ref) if ($control_g == 1);
	print_geno_file ($geno_genotypes_ref, $geno_line_1_ref, $global_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 1) if ($control_g == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	print_geno_file_frqs ($marker_cm_pos_ref, $global_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 1) if ($control_geno_extra == 1 && $control_g == 1);

### FBGW admixture model
	my $fbgw_frqs_geno_ref = fbgw_frqs ($frq_template_geno_ref, $subpop_frq_template_geno_ref, $fbgw_adm_ref,  $family_ref, $chr_ref) if ($control_fbgw == 1);
	print_geno_file ($geno_genotypes_ref, $geno_line_1_ref, $fbgw_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 2) if ($control_fbgw == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	print_geno_file_frqs ($marker_cm_pos_ref, $fbgw_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 2) if ($control_geno_extra == 1 && $control_fbgw == 1);

### FBCW admixture model
	my $fbcw_frqs_geno_ref = fbcw_frqs ($frq_template_geno_ref, $subpop_frq_template_geno_ref, $fbcw_adms_ref, $family_ref, $chr_ref) if ($control_fbcw == 1);
	print_geno_file ($geno_genotypes_ref, $geno_line_1_ref, $fbcw_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 3) if ($control_fbcw == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	print_geno_file_frqs ($marker_cm_pos_ref, $fbcw_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 3) if ($control_geno_extra == 1 && $control_fbcw == 1);

### Local admixture model
	my ($local_adms_geno_ref, $subpop_local_geno_ref) = adms_local ($parameters_ref, $pop_adms_ref, $geno_bps_ref, $chr_ref, $family_ref, $local_adm_extrapolation_ref) if ($control_l == 1);

	my ($local_frqs_geno_ref, $local_adms_geno_pr_ref);
	$local_frqs_geno_ref = local_frqs ($frq_template_geno_ref, $subpop_frq_template_geno_ref, $local_adms_geno_ref,  $family_ref, $chr_ref, $subpop_local_geno_ref, 0) if ($control_l == 1 && $adm_printout == 0);
	($local_frqs_geno_ref, $local_adms_geno_pr_ref) = local_frqs ($frq_template_geno_ref, $subpop_frq_template_geno_ref, $local_adms_geno_ref,  $family_ref, $chr_ref, $subpop_local_geno_ref, 1) if ($control_l == 1 && $adm_printout == 1);
	print_geno_file ($geno_genotypes_ref, $geno_line_1_ref, $local_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 4) if ($control_l == 1); # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	print_geno_file_frqs ($marker_cm_pos_ref, $local_frqs_geno_ref, $parameters_ref, $chr_ref, $family_ref, 4) if ($control_geno_extra == 1 && $control_l == 1);
	print_local_adm_geno_file ($marker_cm_pos_ref, $local_adms_geno_pr_ref, $parameters_ref, $chr_ref, $family_ref) if ($adm_printout == 1 && $control_l == 1);

	} # end of the 0loop for geno file allele frequencies

	} # end of family loop

	} # end of chr loop

	my $date_time = localtime;
	print LOGFILE "	Analysis finished at: $date_time\n";
	print "	Analysis finished at: $date_time\n";
	close LOGFILE;

}


sub warning
{
	my (@warning) = @_;
	if ( scalar @warning ) {
	print "	@warning\n\n";
	print LOGFILE "	@warning\n\n"; 
	}
}


sub error
{
	my (@errors) = @_;
	if ( scalar @errors ) {
	print "	@errors\n\n";
	print LOGFILE "	@errors\n\n"; 
	}
	my $date_time = localtime;
	print "	Program was terminated at: $date_time\n\n";
	print LOGFILE "	Program was terminated at: $date_time\n\n";
	exit 1;
}


sub error_2
{
	my (@errors) = @_;
	if ( scalar @errors ) { warn "@errors\n"; }
	my $date_time = localtime;
	print "	Analysis started at: $date_time\n\n";
	print
	"  \@----------------------------------------------------------\@\n",
	 "  |    ADMIXFRQ!  |     v1.00      |   14/May/2018           |\n",
	 "  |----------------------------------------------------------|\n",
	 "  |Copyright (C) 2018  Rafael A. Nafikov  GNU GPLv3          |\n",
	 "  |This program comes with ABSOLUTELY NO WARRANTY;           |\n",
	 "  |for details see <https://www.gnu.org/licenses/>.          |\n",
	 "  |This is free software, and you are welcome to redistribute|\n",
	 "  |it under certain conditions; for details see              |\n",
	 "  |<https://www.gnu.org/licenses/>.                          |\n",
	 "  |----------------------------------------------------------|\n",
	 "  | For documentation, citation & bug-reporting instructions:|\n",
	 "  |        https://github.com/RafPrograms/ADMIXFRQ           |\n",
	"  \@----------------------------------------------------------\@\n",
	 "\n";
	$date_time = localtime;
	print "	Program was terminated at: $date_time\n\n";
	exit 1;
}



### section for input of data files
sub read_par_file {
	error_2 ("\nError: parameter file was not specified!\n") if (!defined $ARGV[0]);

	my %options = ();
	while (defined(my $arg = shift(@ARGV))) {
	open(PARAM, $arg) or error ("Error: cannot open parameter file $arg!") && error_2 ("Error: cannot open parameter file $arg!");
	while (defined (my $line1 = <PARAM>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if (! defined $line1);
	next if ($line1 =~ m/^#/);
	my @temp0 = split (/#/, $line1);
	goto LABEL1 if (! defined $temp0[0]);

	if ($temp0[0] =~ m/:/) {
	my @temp1 = split (/:/, $temp0[0]);
	$temp1[0] = uc $temp1[0];
	$temp1[0] =~ s/\s+//g;
	$temp1[1] =~ s/^\s+//g if ($temp1[1]);
	$temp1[1] =~ s/\s+$//g if ($temp1[1]);
	$options{$.}{$temp1[0]} = 1 if (defined $temp1[1] && ($temp1[0] ne 'ADMIXTUREPROPORTIONSPRINTOUT' && $temp1[0] ne 'EXTRAPOLATIONOFLOCALADMIXTUREPROPORTIONS' && $temp1[0] ne 'COMPUTEGLOBALADMIXTUREMODELPROPORTIONS'));
	$options{$.}{$temp1[0]} = $temp1[1] if (defined $temp1[1] && ($temp1[0] eq 'ADMIXTUREPROPORTIONSPRINTOUT' || $temp1[0] eq 'EXTRAPOLATIONOFLOCALADMIXTUREPROPORTIONS' || $temp1[0] eq 'COMPUTEGLOBALADMIXTUREMODELPROPORTIONS'));
	my @temp2 = split (/\|/, $temp1[1]);
	for my $el (@temp2)	{
	my @temp3 = split (/=/, $el);
	goto LABEL1 if (!exists $temp3[1]);
	$temp3[1] =~ s/^\s+//g;
	$temp3[1] =~ s/\s+$//g;
	if ($temp3[1] =~ m/^\[/ && $temp3[1] =~ m/\]$/) {
	$temp3[1] =~ s/^\[//g;
	$temp3[1] =~ s/\]$//g;
	my @temp4 = split (/\s+/, $temp3[1]);
	my $count = 0;
	for my $el2 (@temp4) {
	$count++;
	$temp3[0] =~ s/^\s+//g;
	$temp3[0] =~ s/\s+$//g;
	$temp3[0] = uc $temp3[0];
	$temp3[0] =~ s/\s+//g;
	$options{$.}{$temp3[0]}{$count} = $el2;
	}
	} else {
	my @temp4 = split (/\s+/, $temp3[1]);
	my $count = 0;
	for my $el2 (@temp4) {
	$count++;
	$temp3[0] =~ s/^\s+//g;
	$temp3[0] =~ s/\s+$//g;
	$temp3[0] = uc $temp3[0];
	$temp3[0] =~ s/\s+//g;
	$options{$.}{$temp3[0]}{$count} = $el2;
	}
	}
	}
	} else {

	my @temp2 = split (/\|/, $temp0[0]);
	for my $el (@temp2)	{
	my @temp3 = split (/=/, $el);
	goto LABEL1 if (!exists $temp3[1]);
	$temp3[1] =~ s/^\s+//g;
	$temp3[1] =~ s/\s+$//g;
	if ($temp3[1] =~ m/^\[/ && $temp3[1] =~ m/\]$/) {
	$temp3[1] =~ s/^\[//g;
	$temp3[1] =~ s/\]$//g;
	my @temp4 = split (/\s+/, $temp3[1]);
	my $count = 0;
	for my $el2 (@temp4) {
	$count++;
	$temp3[0] =~ s/^\s+//g;
	$temp3[0] =~ s/\s+$//g;
	$temp3[0] = uc $temp3[0];
	$temp3[0] =~ s/\s+//g;
	$options{$.}{$temp3[0]}{$count} = $el2;
	}
	} else {
	my @temp4 = split (/\s+/, $temp3[1]);
	my $count = 0;
	for my $el2 (@temp4) {
	$count++;
	$temp3[0] =~ s/^\s+//g;
	$temp3[0] =~ s/\s+$//g;
	$temp3[0] = uc $temp3[0];
	$temp3[0] =~ s/\s+//g;
	$options{$.}{$temp3[0]}{$count} = $el2;
	}
	}
	}


	}

	LABEL1:
	}
	}

	return \%options;
}


sub read_file_population {

	my ($parameters_ref, $chr_ref) = @_;
	my $chr = $$chr_ref;
	my %parameters = %$parameters_ref;

	my %rs_to_bp = ();
	my %ref_frqs = ();
	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'REFERENCEPOPULATIONNAME'}) {
	my $subpop_name = $parameters{$key1}{'REFERENCEPOPULATIONNAME'}{1};
	error ("Error: population name for 'REFERENCE POPULATION NAME' option was not provided!") if (!defined $subpop_name);
	my $file_path = $parameters{$key1}{'FILEPATH'}{1};

	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	error ("Error: file path for 'REFERENCE POPULATION NAME' option was not provided!") if (!defined $file_path);
	my $header = $parameters{$key1}{'HEADER'}{1};
	my $bp_position = $parameters{$key1}{'BP_POSITION'}{1};
	my $rs_number = $parameters{$key1}{'RS_NUMBER'}{1};
	my $ref_allele = $parameters{$key1}{'REF_ALLELE'}{1};
	my $alt_allele = $parameters{$key1}{'ALT_ALLELE'}{1};
	my $ref_allele_frq = $parameters{$key1}{'REF_ALLELE_FRQ'}{1};

	open(FILE1, $file_path) or error ("Error: cannot open file $file_path containing reference population allele frequencies!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line1);

	if (defined $temp1[$rs_number-1] && defined $temp1[$ref_allele-1] && defined $temp1[$ref_allele_frq-1] && defined $temp1[$alt_allele-1]) {
	$ref_frqs{$chr}{$temp1[$bp_position-1]}{$subpop_name}{'RS_NUMBER'} = $temp1[$rs_number-1];
	$ref_frqs{$chr}{$temp1[$bp_position-1]}{$subpop_name}{'REF_ALLELE'} = $temp1[$ref_allele-1];
	$ref_frqs{$chr}{$temp1[$bp_position-1]}{$subpop_name}{'REF_ALLELE_FRQ'} = $temp1[$ref_allele_frq-1];
	$ref_frqs{$chr}{$temp1[$bp_position-1]}{$subpop_name}{'ALT_ALLELE'} = $temp1[$alt_allele-1];

	$rs_to_bp{$chr}{$temp1[$rs_number-1]} = $temp1[$bp_position-1] if (!defined $rs_to_bp{$chr}{$temp1[$rs_number-1]});
	}
	}
	close FILE1;
	}
	}

	return \%ref_frqs, \%rs_to_bp;
}


sub read_file_admixtures_chr {

	my ($parameters_ref, $chr_ref) = @_;
	my $chr = $$chr_ref;
	my %parameters = %$parameters_ref;

	my %adms = ();
	my $var_identifier;

	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'LOCALADMIXTUREPERFAMILYBYPOPULATION'}) {
	my $subpop_name = $parameters{$key1}{'LOCALADMIXTUREPERFAMILYBYPOPULATION'}{1};
	my $file_path = $parameters{$key1}{'FILEPATH'}{1};
	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	my $header = $parameters{$key1}{'HEADER'}{1};
	my $var_id = $parameters{$key1}{'VARIANTIDENTIFIER'}{1};
	$var_identifier = $var_id if (!defined $var_identifier);
	my $var_id_col = $parameters{$key1}{'VARIANTIDENTIFIERCOLUMN'}{1};
	my @families;
	my @col_nums;
	for my $key2 (sort keys %{$parameters{$key1}{'FAMILIES'}}) { # consecutive numbers from 1 to the total number of entries
	push @families, $parameters{$key1}{'FAMILIES'}{$key2};
	push @col_nums, $parameters{$key1}{'COLUMNNUMBERS'}{$key2};
	}

	open(FILE1, $file_path) or error ("Error: cannot open file $file_path containing local admixture proportions for reference population!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line1);

	for (my $el = 0; $el <= $#families; $el++) {
	$adms{$chr}{$temp1[$var_id_col-1]}{$subpop_name}{$families[$el]} = $temp1[$col_nums[$el]-1];
	}

	}
	close FILE1;

	}
	}

	# check if admixture proportions across all populations add up to 1
	for my $key1 (sort {$a<=>$b} keys %adms) { # chr
	for my $key2 (sort {$a cmp $b} keys %{$adms{$key1}}) { # var identifier which is either rs number or bp position
	my %temp1 = ();
	for my $key3 (sort {$a cmp $b} keys %{$adms{$key1}{$key2}}) { # reference population name
	for my $key4 (sort {$a cmp $b} keys %{$adms{$key1}{$key2}{$key3}}) { # family
	push @{$temp1{$key4}}, $adms{$key1}{$key2}{$key3}{$key4};
	}
	}

	for my $key5 (sort {$a cmp $b} keys %temp1) { # family again
	my @temp2 = @{$temp1{$key5}};
	my $sum;
	foreach (@temp2) {
	$sum += $_;
	}
	if ($sum <= 1.01 && $sum >= 0.99) {
	# do nothing
	} else {
	error ("Sum of admixture proportions across specified reference populations for varriant $key2 on chromosome $key1 in family $key5 does not equal 1.\n");
	}
	}

	}
	}

	return \%adms, \$var_identifier;
}


sub read_file_admixtures_genome {

	my ($parameters_ref) = $_[0];
	my %parameters = %$parameters_ref;

	my @chrs = (1..22);
	my %adms = ();
	my $var_identifier;

	for my $chr (@chrs) {
	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'LOCALADMIXTUREPERFAMILYBYPOPULATION'}) {
	my $subpop_name = $parameters{$key1}{'LOCALADMIXTUREPERFAMILYBYPOPULATION'}{1};
	my $file_path = $parameters{$key1}{'FILEPATH'}{1};
	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	my $header = $parameters{$key1}{'HEADER'}{1};
	my $var_id = $parameters{$key1}{'VARIANTIDENTIFIER'}{1};
	$var_identifier = $var_id if (!defined $var_identifier);
	my $var_id_col = $parameters{$key1}{'VARIANTIDENTIFIERCOLUMN'}{1};

	my @families;
	my @col_nums;
	for my $key2 (sort keys %{$parameters{$key1}{'FAMILIES'}}) { # consecutive numbers from 1 to the total number of entries
	push @families, $parameters{$key1}{'FAMILIES'}{$key2};
	push @col_nums, $parameters{$key1}{'COLUMNNUMBERS'}{$key2};
	}

	open(FILE1, $file_path) or error ("Error: cannot open file $file_path containing local admixture proportions for reference population!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line1);

	for (my $el = 0; $el <= $#families; $el++) {
	$adms{$chr}{$temp1[$var_id_col-1]}{$subpop_name}{$families[$el]} = $temp1[$col_nums[$el]-1];
	}

	}
	close FILE1;
	}
	}
	}


	# check if admixture proportions across all populations add up to 1
	for my $key1 (sort {$a<=>$b} keys %adms) { # chr
	for my $key2 (sort {$a cmp $b} keys %{$adms{$key1}}) { # var identifier which is either rs number or bp position
	my %temp1 = ();
	for my $key3 (sort {$a cmp $b} keys %{$adms{$key1}{$key2}}) { # reference population name
	for my $key4 (sort {$a cmp $b} keys %{$adms{$key1}{$key2}{$key3}}) { # family
	push @{$temp1{$key4}}, $adms{$key1}{$key2}{$key3}{$key4};
	}
	}

	for my $key5 (sort {$a cmp $b} keys %temp1) { # family again
	my @temp2 = @{$temp1{$key5}};
	my $sum;
	foreach (@temp2) {
	$sum += $_;
	}
	if ($sum <= 1.01 && $sum >= 0.99) {
	# do nothing
	} else {
	error ("Sum of admixture proportions across specified reference populations for varriant $key2 on chromosome $key1 in family $key5 does not equal 1.\n");
	}
	}

	}
	}

	return \%adms, \$var_identifier;
}


sub get_bp_for_adm_from_ref_pop {
	my ($adms_rs_ref, $rs_to_bp_ref) = @_;
	my %adms_rs = %$adms_rs_ref;
	my %rs_to_bp = %$rs_to_bp_ref;

	my %adms = ();

	for my $key1 (sort {$a <=> $b} keys %adms_rs) { # chromosome
	for my $key2 (sort {$a cmp $b} keys %{$adms_rs{$key1}}) { # rs_number
	if (defined $rs_to_bp{$key1}{$key2}) {
	for my $key3 (sort {$a cmp $b} keys %{$adms_rs{$key1}{$key2}}) { # population
	for my $key4 (sort {$a cmp $b} keys %{$adms_rs{$key1}{$key2}{$key3}}) { # families
	$adms{$key1}{$rs_to_bp{$key1}{$key2}}{$key3}{$key4} = $adms_rs{$key1}{$key2}{$key3}{$key4};
	}
	}
	} else {
	print LOGFILE2 "rs number $key2 for variant on chromosome $key1 was not matched to the same rs number in reference population dataset because either this rs number had changed during genome build update or does not exist.\n";
	}
	}
	}

	return \%adms;
}


sub read_geno_file_info {

	my ($parameters_ref, $chr_ref, $family_ref) = @_;
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;

	my %marker_cm_pos = ();
	my %marker_info = ();
	my %genotypes = ();
	my $geno_line_1;
	my @geno_bps;

	my $file_path;
	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	for my $key2 (sort keys %{$parameters{$key1}}) { # options from specific line of parameter file
	if ($key2 eq 'GENOFILEFORGL_AUTO') {
	$file_path = $parameters{$key1}{'FILEPATH'}{1};
	}
	}
	}

	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	$file_path =~ s/\+/$family/g;

	open(FILE1, $file_path) or error ("Error: cannot open geno file $file_path!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;

	next if (!$line1);
	$geno_line_1 = $line1 if ($line1 =~ m/map/);
	my @temp1 = split (/\D+\s+/, $line1) if ($. == 1);
	my @temp2 = split (/\s+/, $temp1[1]) if ($. == 1);
	for my $el (@temp2) {
	$marker_cm_pos{$el} = 1;
	}

	if ($line1 !~ m/marker/ && $. != 1) {
	my @temp3 = split (/\s+/, $line1);
	$genotypes{$.} = [@temp3];
	}


	}
	close FILE1;


	my $file_path_2;
	for my $key1 (sort keys %parameters) {
	for my $key2 (sort keys %{$parameters{$key1}}) {
	if ($key2 eq 'INFORMATIONFORGENOFILE') {
	$file_path_2 = $parameters{$key1}{'FILEPATH'}{1};
	error ("Error: location of file with extra information for geno file is not provided!") if (!defined $file_path_2);
	$file_path_2 =~ s/^"//g;
	$file_path_2 =~ s/"$//g;
	$file_path_2 =~ s/^'//g;
	$file_path_2 =~ s/'$//g;
	$file_path_2 =~ s/\*/$chr/g;
	$file_path_2 =~ s/\+/$family/g;

	my $header = $parameters{$key1}{'HEADER'}{1};
	my $bp_position = $parameters{$key1}{'BP_POSITION'}{1};
	my $cm_position = $parameters{$key1}{'CM_POSITION'}{1};
	my $rs_number = $parameters{$key1}{'RS_NUMBER'}{1};


	open(FILE2, $file_path_2) or error ("Error: cannot open file $file_path containing additional information for geno file!");
	while (defined (my $line2 = <FILE2>)) {
	chomp $line2;
	$line2 =~ s/^\s+//g;
	$line2 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line2);
	if (defined $temp1[$cm_position-1] && defined $temp1[$rs_number-1]) {
	$marker_info{$temp1[$bp_position-1]}{'CM_POSITION'} = $temp1[$cm_position-1];
	$marker_info{$temp1[$bp_position-1]}{'RS_NUMBER'} = $temp1[$rs_number-1];
	}
	}
	close FILE2;

	}
	}
	} # 2nd major loop end


	my %geno_rs_bp = ();
	for my $key1 (sort {$a <=> $b} keys %marker_info) { # bp position
	if (defined $marker_cm_pos{$marker_info{$key1}{'CM_POSITION'}}) {
	$geno_rs_bp{$key1} = $marker_info{$key1}{'RS_NUMBER'};
	push @geno_bps, $key1;
	}
	}

	if (scalar keys %marker_cm_pos != scalar @geno_bps || scalar @geno_bps != scalar keys %geno_rs_bp) {
	error ("Geno file $file_path and/or its corresponding information file are incomplete. Please check your files for completeness.");
	}
 
	return \%marker_cm_pos, \%geno_rs_bp, \%genotypes, \$geno_line_1, \@geno_bps;
}


sub read_file_sample_allele_frqs {

	my ($parameters_ref, $chr_ref) = @_;
	my $chr = $$chr_ref;
	my %parameters = %$parameters_ref;

	my %sample_frqs = ();

	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'YOURSAMPLEALLELEFREQUENCIES'}) {
	my $file_path = $parameters{$key1}{'FILEPATH'}{1};
	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	my $header = $parameters{$key1}{'HEADER'}{1};
	my $bp_position = $parameters{$key1}{'BP_POSITION'}{1};
	my $al1 = $parameters{$key1}{'REF_ALLELE'}{1};
	my $al1_frq = $parameters{$key1}{'REF_ALLELE_FRQ'}{1};
	my $al2 = $parameters{$key1}{'ALT_ALLELE'}{1};
	open(FILE1, $file_path) or error ("Error: cannot open file $file_path containing your sample specific allele frequencies!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line1);
	if (defined $temp1[$al1-1] && defined $temp1[$al1_frq-1] && defined $temp1[$al2-1]) {
	$sample_frqs{$chr}{$temp1[$bp_position-1]}{'al1'} = $temp1[$al1-1];
	$sample_frqs{$chr}{$temp1[$bp_position-1]}{'al1_frq'} = $temp1[$al1_frq-1];
	$sample_frqs{$chr}{$temp1[$bp_position-1]}{'al2'} = $temp1[$al2-1];
	$sample_frqs{$chr}{$temp1[$bp_position-1]}{'al2_frq'} = sprintf ("%.6f",(1 - $temp1[$al1_frq-1]));
	}
	}
	close FILE1;

	}
	}

return \%sample_frqs;
}


sub read_bp_modify_frq {

	my ($parameters_ref, $chr_ref) = @_;
	my $chr = $$chr_ref;
	my %parameters = %$parameters_ref;

	my @bp_pos = ();
	my %sample_bp_rs = ();
	my $rs_num_check = 0;

	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'MODIFYYOURBPALLELEFREQUENCIES'}) {
	my $file_path = $parameters{$key1}{'FILEPATH'}{1};
	$file_path =~ s/^"//g;
	$file_path =~ s/"$//g;
	$file_path =~ s/^'//g;
	$file_path =~ s/'$//g;
	$file_path =~ s/\*/$chr/g;
	my $header = $parameters{$key1}{'HEADER'}{1};
	my $bp_position = $parameters{$key1}{'BP_POSITION'}{1};
	my $rs_number = $parameters{$key1}{'RS_NUMBER'}{1} if (defined $parameters{$key1}{'RS_NUMBER'}{1});
	$rs_num_check = 1 if (defined $parameters{$key1}{'RS_NUMBER'}{1});

	open(FILE1, $file_path) or error ("Error: cannot open file $file_path containing base pair positions for which allele frequencies were requested to be modified!");
	while (defined (my $line1 = <FILE1>)) {
	chomp $line1;
	$line1 =~ s/^\s+//g;
	$line1 =~ s/\s+$//g;
	next if ($header eq 'T' && $. == 1);
	my @temp1 = split (/\s+/, $line1);
	push @bp_pos, $temp1[$bp_position-1];
	$sample_bp_rs{$chr}{$temp1[$bp_position-1]} = $temp1[$rs_number-1] if ($rs_num_check == 1);
	}
	close FILE1;
	}
	}

return \@bp_pos, 7 if ($rs_num_check == 0);
return \@bp_pos, \%sample_bp_rs if ($rs_num_check == 1);
}
### end of section for input of data files



### section containing subroutines for data print out into files
sub make_dir {
	my ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref) = @_;
# $comp_type = 0 if frequencies for geno files are needed; 1 = if frequencies not for geno files are needed (FBGW, FBCW, Local); 2 = if admixture proportions for Local model are needed; 3 = if admixture proportions for Global, FBGW, or FBCW model are needed; 4 = directory for error log file; 5 = if frequencies not for geno files are needed (Global); 6 = printout of local admixture proportions for each variant in geno file
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref if ($chr_ref != 0);
	my $family = $$family_ref if ($family_ref != 0);
	my $ancestry_model = $$ancestry_model_ref if ($ancestry_model_ref != 0);
	my $comp_type = $$comp_type_ref;

	my $dir_out;
	for my $key1 (sort keys %parameters) { # line numbers in parameter file
	if (defined $parameters{$key1}{'OUTPUTDIRECTORY'}) {
	$dir_out = $parameters{$key1}{'OUTPUTDIRECTORY'}{1};
	$dir_out =~ s/^"//g;
	$dir_out =~ s/"$//g;
	$dir_out =~ s/^'//g;
	$dir_out =~ s/'$//g;
	}
	}

	if (!defined $dir_out) {error_2 ("\n	Error: Output directory is not specified!\n	Log file was not created!\n")};

	my $dir0;
	my $dir1;
	my $dir2;
	my $dir3;
	my $dir4;
	my $dir5;
	my $dir6;

	if ($comp_type == 0) {
	$dir0 = "$dir_out/gl_auto/$ancestry_model/$family/chr$chr";
	if (! -e $dir0) {
	`mkdir -p $dir0`;
	}
	return \$dir0;
	}

	elsif ($comp_type == 1) {
	$dir1 = "$dir_out/allele_frequencies/$ancestry_model/$family/chr$chr";
	if (! -e $dir1) {
	`mkdir -p $dir1`;
	}
	return \$dir1;
	}

	elsif ($comp_type == 2) {
	$dir2 = "$dir_out/admixture_proportions/$ancestry_model/$family/chr$chr";
	if (! -e $dir2) {
	`mkdir -p $dir2`;
	}
	return \$dir2;
	}

	elsif ($comp_type == 3) {
	$dir3 = "$dir_out/admixture_proportions/$ancestry_model";
	if (! -e $dir3) {
	`mkdir -p $dir3`;
	}
	return \$dir3;
	}

	elsif ($comp_type == 4) {
	$dir4 = "$dir_out";
	if (! -e $dir4) {
	`mkdir -p $dir4`;
	}
	return \$dir4;
	}

	elsif ($comp_type == 5) {
	$dir5 = "$dir_out/allele_frequencies/$ancestry_model/chr$chr";
	if (! -e $dir5) {
	`mkdir -p $dir5`;
	}
	return \$dir5;
	}

	elsif ($comp_type == 6) {
	$dir6 = "$dir_out/gl_auto_local_adm/$family/chr$chr";
	if (! -e $dir6) {
	`mkdir -p $dir6`;
	}
	return \$dir6;
	}

}


sub print_geno_file {

	my ($geno_genotypes_ref, $geno_line_1_ref, $geno_frqs_ref, $parameters_ref, $chr_ref, $family_ref, $ancestry_m) = @_;

	my %geno_genotypes = %$geno_genotypes_ref;
	my $geno_line_1 = $$geno_line_1_ref;
	my %geno_frqs = %$geno_frqs_ref;
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;

	my $ancestry_model; # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	if ($ancestry_m == 1) {
	$ancestry_model = 'Global';
	}
	elsif ($ancestry_m == 2) {
	$ancestry_model = 'FBGW';
	}
	elsif ($ancestry_m == 3) {
	$ancestry_model = 'FBCW';
	}
	elsif ($ancestry_m == 4) {
	$ancestry_model = 'Local';
	}
	my $ancestry_model_ref = \$ancestry_model;

	my $comp_type = 0;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/chr$chr.geno");
	print FILEOUT1 "$geno_line_1\n\n";

	my @geno_frqs_out;
	for my $key1 (sort {$a <=> $b} keys %{$geno_frqs{$chr}}) { # bp
	push @geno_frqs_out, "$geno_frqs{$chr}{$key1}{1}[1]\n";
	}
	my $number_m = scalar @geno_frqs_out;

	for (my $el = 0; $el <= $#geno_frqs_out; $el++) {
	my $temp1 = $el + 1;
	my $temp2 = sprintf ("%.6f", $geno_frqs_out[$el]);
	my $temp3 = 1 - $temp2;
	$temp3 = sprintf ("%.6f", $temp3);
	print FILEOUT1 "set markers $temp1 allele freq $temp2 $temp3\n";
	}
	print FILEOUT1 "\n";
	print FILEOUT1 "set markers $number_m data\n\n";

	for my $key1 (sort {$a <=> $b} keys %geno_genotypes) {
	print FILEOUT1 "@{$geno_genotypes{$key1}}\n";
	}

close FILEOUT1;

}


sub print_geno_file_frqs {

	my ($marker_cm_pos_ref, $geno_frqs_ref, $parameters_ref, $chr_ref, $family_ref, $ancestry_m) = @_;
	my %marker_cm_pos = %$marker_cm_pos_ref;
	my %geno_frqs = %$geno_frqs_ref;
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;

	my $ancestry_model; # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	if ($ancestry_m == 1) {
	$ancestry_model = 'Global';
	}
	elsif ($ancestry_m == 2) {
	$ancestry_model = 'FBGW';
	}
	elsif ($ancestry_m == 3) {
	$ancestry_model = 'FBCW';
	}
	elsif ($ancestry_m == 4) {
	$ancestry_model = 'Local';
	}
	my $ancestry_model_ref = \$ancestry_model;

	my $comp_type = 0;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/chr$chr.frq");
	print FILEOUT1 "bp_position cM_position ref_allele alt_allele ref_allele_frq alt_allele_frq\n";

	my @cm_pos;
	for my $key1 (sort {$a <=> $b} keys %marker_cm_pos) { # cM positions
	push @cm_pos, $key1;
	}

	my $count = 0;
	for my $key1 (sort {$a <=> $b} keys %{$geno_frqs{$chr}}) { # bp position
	print FILEOUT1 " $key1 $cm_pos[$count] $geno_frqs{$chr}{$key1}{1}[0] $geno_frqs{$chr}{$key1}{2}[0] $geno_frqs{$chr}{$key1}{1}[1] $geno_frqs{$chr}{$key1}{2}[1]\n";
	$count++;
	}

}


sub print_local_adm_geno_file {

	my ($marker_cm_pos_ref, $local_adms_geno_ref, $parameters_ref, $chr_ref, $family_ref) = @_;
	my %marker_cm_pos = %$marker_cm_pos_ref;
	my %local_adms_geno = %$local_adms_geno_ref;
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;

	my $ancestry_model = 'Local';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 6;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/chr$chr.local.adm.txt");
	print FILEOUT1 "bp_position cM_position";

	my @cm_pos;
	for my $key1 (sort {$a <=> $b} keys %marker_cm_pos) {
	push @cm_pos, $key1;
	}

	for my $key1 (sort {$a <=> $b} keys %{$local_adms_geno{$chr}}) { # bp position
	for my $key2 (sort {$a cmp $b} keys %{$local_adms_geno{$chr}{$key1}}) { # population
	print FILEOUT1 " local_adm_$key2";  
	}
	last;
	}
	print FILEOUT1 "\n";

	my $count = 0;
	for my $key1 (sort {$a <=> $b} keys %{$local_adms_geno{$chr}}) { # bp position
	print FILEOUT1 " $key1 $cm_pos[$count]";
	$count++;
	for my $key2 (sort {$a cmp $b} keys %{$local_adms_geno{$chr}{$key1}}) { # population
	if ($local_adms_geno{$chr}{$key1}{$key2} == 0 || $local_adms_geno{$chr}{$key1}{$key2} == 1) {
	my $temp0 = sprintf ("%.0f", $local_adms_geno{$chr}{$key1}{$key2});
	print FILEOUT1 " $temp0";
	} else {
	my $temp1 = sprintf ("%.6f", $local_adms_geno{$chr}{$key1}{$key2});
	print FILEOUT1 " $temp1";
	}
	}
		print FILEOUT1 "\n";
	}

}


sub print_local_adm_file {

	my ($local_adms_ref, $subpop_local_ref, $parameters_ref, $chr_ref, $family_ref) = @_;
	my %local_adms = %$local_adms_ref;
	my %subpop_local = %$subpop_local_ref;
	my %parameters = %$parameters_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;

	my $ancestry_model = 'Local';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 2;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/chr$chr.$family.txt");
	print FILEOUT1 "bp_position";

	for my $key (sort { $subpop_local{$a} <=> $subpop_local{$b} } keys %subpop_local) {
	print FILEOUT1 " local_adm_$key";  
	}
	print FILEOUT1 "\n";


	for my $key ( sort {$a <=> $b} keys %local_adms) {
	my @temp1 = split (/\s+/, "@{$local_adms{$key}}");
	my @temp2;
	for my $el (@temp1) {
	if ($el == 0 || $el == 1) {
	push @temp2, $el;
	} else {
	push @temp2, sprintf ("%.6f", $el);
	}
	}
	print FILEOUT1 "$key @temp2\n";
	}
	close FILEOUT1;

}


sub remove_old_fbcw_adm_file {

	my ($parameters_ref, $chr_ref, $family_ref) = @_;
	my %parameters = %$parameters_ref;

	my $ancestry_model = 'FBCW';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 3;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;
	open (FILEOUT1, ">", "$dir/fbcw.admixture.txt");
	close FILEOUT1;

}
 

sub print_fbcw_adm_file {

	my ($fbcw_adm_ref, $parameters_ref, $chr_ref, $family_ref) = @_;
	my %fbcw_adm = %$fbcw_adm_ref;
	my %parameters = %$parameters_ref;

	my $ancestry_model = 'FBCW';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 3;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	my @subpop;
	for my $key1 (sort {$a <=> $b} keys %fbcw_adm) { # chr
	for my $key2 (sort {$a cmp $b} keys %{$fbcw_adm{$key1}}) {# family
	for my $key3 (sort {$a cmp $b} keys %{$fbcw_adm{$key1}{$key2}}) {# population
	push @subpop, $key3;
	}
	last;
	}
	last;
	}

	open (FILEOUT1, ">", "$dir/fbcw.admixture.txt") if ($chr_ref == 0);
	open (FILEOUT1, ">>", "$dir/fbcw.admixture.txt") if ($chr_ref != 0);
	print FILEOUT1 "chr family";
	for my $el (@subpop) {
	print FILEOUT1 " FBCW_adm_$el";
	}
	print FILEOUT1 "\n";

	for my $key1 (sort {$a <=> $b} keys %fbcw_adm) { # chr
	for my $key2 (sort {$a cmp $b} keys %{$fbcw_adm{$key1}}) {# family
	print FILEOUT1 "$key1 $key2";
	for my $key3 (sort {$a cmp $b} keys %{$fbcw_adm{$key1}{$key2}}) {# population
	if ($fbcw_adm{$key1}{$key2}{$key3} == 0 || $fbcw_adm{$key1}{$key2}{$key3} == 1) {
	print FILEOUT1 " $fbcw_adm{$key1}{$key2}{$key3}"; 
	} else {
	my $temp1 = sprintf ("%.6f", $fbcw_adm{$key1}{$key2}{$key3});
	print FILEOUT1 " $temp1";
	}
	}
	print FILEOUT1 "\n";
	}
	}
	close FILEOUT1;

}


sub print_fbgw_adm_file {

	my ($fbgw_adm_ref, $parameters_ref) = @_;
	my %fbgw_adm = %$fbgw_adm_ref;
	my $ancestry_model = 'FBGW';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 3;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, 0, 0, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	my @subpop;
	for my $key1 ( sort {$a cmp $b} keys %fbgw_adm) {
	for my $key2 ( sort {$a cmp $b} keys %{ $fbgw_adm{$key1}} ) {
	push @subpop, $key2;
	}
	last;
	}

	open (FILEOUT1, ">", "$dir/fbgw.admixture.txt");
	print FILEOUT1 "family";
	for my $el (@subpop) {
	print FILEOUT1 " FBGW_adm_$el";
	}
	print FILEOUT1 "\n";


	for my $key1 ( sort {$a cmp $b} keys %fbgw_adm) {
	print FILEOUT1 "$key1";
	for my $key2 ( sort {$a cmp $b} keys %{ $fbgw_adm{$key1}} ) {

	if ($fbgw_adm{$key1}{$key2} == 0 || $fbgw_adm{$key1}{$key2} == 1) {
	print FILEOUT1 " $fbgw_adm{$key1}{$key2}";
	} else {
	my $temp1 = sprintf ("%.6f", $fbgw_adm{$key1}{$key2});
	print FILEOUT1 " $temp1";
	}
	}
	print FILEOUT1 "\n";
	}
	close FILEOUT1;

}


sub print_global_adm_file {

	my ($global_ref, $parameters_ref) = @_;
	my %global = %$global_ref;
	my %parameters_ref = %$parameters_ref;

	my $ancestry_model = 'Global';
	my $ancestry_model_ref = \$ancestry_model;
	my $comp_type = 3;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, 0, 0, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/global.admixture.txt");
	print FILEOUT1 "population global_adm\n";
	for my $key1 ( sort {$a cmp $b} keys %global) {
	if ($global{$key1} == 0 || $global{$key1} == 1) {
	print FILEOUT1 "$key1 $global{$key1}\n";
	} else {
	my $temp1 = sprintf ("%.6f", $global{$key1});
	print FILEOUT1 "$key1 $temp1\n";
	}
	}
	close FILEOUT1;

}


sub print_frqs {

	my ($parameters_ref, $frqs_ref, $chr_ref, $family_ref, $ancestry_m) = @_;
	my %parameters = %$parameters_ref;
	my %frqs = %$frqs_ref;
	my $chr = $$chr_ref;

	my $ancestry_model; # 1 = Global; 2 = FBGW; 3 = FBCW; 4 = Local
	if ($ancestry_m == 1) {
	$ancestry_model = 'Global';
	}
	elsif ($ancestry_m == 2) {
	$ancestry_model = 'FBGW';
	}
	elsif ($ancestry_m == 3) {
	$ancestry_model = 'FBCW';
	}
	elsif ($ancestry_m == 4) {
	$ancestry_model = 'Local';
	}
	my $ancestry_model_ref = \$ancestry_model;

	my $comp_type;
	$comp_type = 1;
	$comp_type = 5 if ($ancestry_m == 1);
	my $comp_type_ref = \$comp_type;

	my $dir_ref = make_dir ($parameters_ref, $chr_ref, $family_ref, $ancestry_model_ref, $comp_type_ref);
	my $dir = $$dir_ref;

	open (FILEOUT1, ">", "$dir/global.frq.txt") if ($ancestry_m == 1);
	print FILEOUT1 "bp_position ref_allele alt_allele ref_allele_frq alt_allele_frq\n" if ($ancestry_m == 1);

	open (FILEOUT1, ">", "$dir/fbgw.frq.txt") if ($ancestry_m == 2);
	print FILEOUT1 "bp_position ref_allele alt_allele ref_allele_frq alt_allele_frq\n" if ($ancestry_m == 2);

	open (FILEOUT1, ">", "$dir/fbcw.frq.txt") if ($ancestry_m == 3);
	print FILEOUT1 "bp_position ref_allele alt_allele ref_allele_frq alt_allele_frq\n" if ($ancestry_m == 3);

	open (FILEOUT1, ">", "$dir/local.frq.txt") if ($ancestry_m == 4);
	print FILEOUT1 "bp_position ref_allele alt_allele ref_allele_frq alt_allele_frq\n" if ($ancestry_m == 4);

	for my $key1 (sort {$a <=> $b} keys %{ $frqs{$chr}}) { # bp position
	print FILEOUT1 "$key1 $frqs{$chr}{$key1}{1}[0] $frqs{$chr}{$key1}{2}[0] $frqs{$chr}{$key1}{1}[1] $frqs{$chr}{$key1}{2}[1]\n";
	}
	close FILEOUT1;

}
### end of section containing subroutines for data print out into files



### section for computation of admixture proportions
sub adms_fbcw { # compute FBCW admixture proportions only for allele frequency computations

	my $adm_ref = $_[0];
	my %adm = %$adm_ref;
	my %fbcw = ();
	my %fbcw_out = ();

	for my $key1 (sort keys %adm) { # chr
	for my $key2 (sort keys %{ $adm{$key1} }) { # bp position
	for my $key3 (sort keys %{ $adm{$key1}{$key2} }) { # population
	for my $key4 (sort keys %{ $adm{$key1}{$key2}{$key3} }) { # family
	push @{$fbcw{$key1}{$key4}{$key3}}, $adm{$key1}{$key2}{$key3}{$key4};
	}
	}
	}
	}


	for my $key1 (sort keys %fbcw) {
	for my $key2 (sort keys %{ $fbcw{$key1} } ) {
	for my $key3 (sort keys %{ $fbcw{$key1}{$key2} } ) { # $key1 = chr, $key2 = family, $key3 = population
	my @temp1 = @{$fbcw{$key1}{$key2}{$key3}};
	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$fbcw_out{$key1}{$key2}{$key3} = $result;
	}
	}
	}

return \%fbcw_out;
}


sub adm_global_supplied {

	my $parameters_ref = $_[0];
	my %parameters = %$parameters_ref;

	my %global_adm = ();
	for my $key1 (sort keys %parameters) {
	for my $key2 (sort keys %{$parameters{$key1}}) {
	if ($key2 eq 'GLOBALADMIXTUREMODELPROPORTIONS') {
	for my $key3 (sort keys %{$parameters{$key1}}) {
	$global_adm{$key3} = $parameters{$key1}{$key3}{1} if ($key3 ne 'GLOBALADMIXTUREMODELPROPORTIONS');
	}
	}
	}
	}

	return \%global_adm;
}


sub adms_global_fbgw { # computes Global and FBGW admixture proportions for print out into a file and for allele frequency computations

	my ($average_type_ref, $adm_ref, $global_adm_comp_ref) = @_;
	my $average_type = $$average_type_ref; # 1 = weighted, 2 = normal
	my %adm = %$adm_ref;
	my $global_adm_comp = $$global_adm_comp_ref; # 0 = use specified Global admixture; 1 = compute Global admixture from the data

	my %fbgw = (); # weighted
	my %fbgw_out = (); # weighted
	my %fbcw = ();
	my %fbcw_out = ();
	my %fbgw_n = (); # normal
	my %fbgw_n_out = (); # normal
	my %global = ();
	my %global_out = ();
	my %global_n = ();
	my %global_n_out = ();


	for my $key1 (sort keys %adm) { # chr
	for my $key2 (sort keys %{ $adm{$key1} }) { # bp position
	for my $key3 (sort keys %{ $adm{$key1}{$key2} }) { # population
	for my $key4 (sort keys %{ $adm{$key1}{$key2}{$key3} }) { # family
	push @{$fbgw{$key4}{$key3}}, $adm{$key1}{$key2}{$key3}{$key4} if ($average_type == 1);
	push @{$fbcw{$key1}{$key4}{$key3}}, $adm{$key1}{$key2}{$key3}{$key4} if ($average_type == 2);
	}
	}
	}
	}


	if ($average_type == 2) {

	for my $key1 (sort keys %fbcw) { # chr
	for my $key2 (sort keys %{ $fbcw{$key1} } ) { # family
	for my $key3 (sort keys %{ $fbcw{$key1}{$key2} } ) { # population
	my @temp1 = @{$fbcw{$key1}{$key2}{$key3}};
	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$fbcw_out{$key1}{$key2}{$key3} = $result;
	}
	}
	}

	for my $key1 (sort keys %fbcw_out) { # chr
	for my $key2 (sort keys %{ $fbcw_out{$key1} } ) { # family
	for my $key3 (sort keys %{ $fbcw_out{$key1}{$key2} } ) { # population
	push @{$fbgw_n{$key2}{$key3}}, $fbcw_out{$key1}{$key2}{$key3};
	}
	}
	}


	for my $key1 (sort keys %fbgw_n) { # family
	for my $key2 (sort keys %{ $fbgw_n{$key1} } ) { # population
	my @temp1 = @{$fbgw_n{$key1}{$key2}};
	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$fbgw_n_out{$key1}{$key2} = $result;
	}
	}


	if ($global_adm_comp == 1 ) {
	for my $key1 (sort keys %fbgw_n_out) { # family
	for my $key2 (sort keys %{ $fbgw_n_out{$key1} }) { # population
	push @{$global_n{$key2}}, $fbgw_n_out{$key1}{$key2};
	}
	}


	for my $key1 (sort keys %global_n) { # population
	my @temp1 = @{$global_n{$key1}};

	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$global_n_out{$key1} = $result;
	}
	}

	}


	if ($average_type == 1) {
	for my $key1 (sort keys %fbgw) { # family
	for my $key2 (sort keys %{ $fbgw{$key1} }) { # population
	my @temp1 = @{$fbgw{$key1}{$key2}};

	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$fbgw_out{$key1}{$key2} = $result;
	}
	}


	for my $key1 (sort keys %fbgw) { # family
	for my $key2 (sort keys %{ $fbgw{$key1} }) { # population
	push @{$global{$key2}}, @{$fbgw{$key1}{$key2}};
	}
	}


	if ($global_adm_comp == 1) {
	for my $key1 (sort keys %global) { # population
	my @temp1 = @{$global{$key1}};

	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$global_out{$key1} = $result;
	}
	}

	}

	return \(%fbgw_out, %global_out) if ($average_type == 1 && $global_adm_comp == 1);
	return \(%fbgw_n_out, %global_n_out) if ($average_type == 2 && $global_adm_comp == 1);

	return \%fbgw_out if ($average_type == 1 && $global_adm_comp == 0);
	return \%fbgw_n_out if ($average_type == 2 && $global_adm_comp == 0);

}


sub adms_fbcw_b {

	my ($adm_ref) = @_;
	my %adm = %$adm_ref;

	my %fbcw = ();
	my %fbcw_out = ();

	for my $key1 (sort keys %adm) { # chr
	for my $key2 (sort keys %{ $adm{$key1} }) { # bp position
	for my $key3 (sort keys %{ $adm{$key1}{$key2} }) { # population
	for my $key4 (sort keys %{ $adm{$key1}{$key2}{$key3} }) { # family
	push @{$fbcw{$key1}{$key4}{$key3}}, $adm{$key1}{$key2}{$key3}{$key4};
	}
	}
	}
	}


	for my $key1 (sort keys %fbcw) { # chr
	for my $key2 (sort keys %{ $fbcw{$key1} } ) { # family
	for my $key3 (sort keys %{ $fbcw{$key1}{$key2} } ) { # population
	my @temp1 = @{$fbcw{$key1}{$key2}{$key3}};
	my $sum;
	foreach (@temp1) {
	$sum += $_;
	}
	my $result = $sum/ scalar @temp1;
	$fbcw_out{$key1}{$key2}{$key3} = $result;
	}
	}
	}
	return \%fbcw_out;
}


sub adms_local {

	my ($parameters_ref, $adm_ref, $bp_pos_ref, $chr_ref, $family_ref, $extrapolation_ref) = @_; # $extrapolation = 0 (yes, default value) or 1 (no extrapolation) 
	my %parameters = %$parameters_ref;
	my %adm = %$adm_ref;
	my $chr = $$chr_ref;
	my $family = $$family_ref;
	my @bp_pos = @$bp_pos_ref;
	@bp_pos =  sort @bp_pos;
	my $extrapolation = $$extrapolation_ref;

	my @aoa =();
	my %local_adms = ();

	for my $key2 (sort {$a <=> $b} keys %{ $adm{$chr} }) { # bp_position
	my @row_aoa =();
	push @row_aoa, $key2;
	for my $key3 (sort {$a cmp $b} keys %{ $adm{$chr}{$key2} }) { # population
	push @row_aoa, $adm{$chr}{$key2}{$key3}{$family};

	}
	push @aoa, [@row_aoa];

	}


	my @subpops;
 	my %subpop = ();
	my $count1 = 0;
	for my $key2 (sort {$a <=> $b} keys %{ $adm{$chr} }) { # bp_position
	for my $key3 (sort {$a cmp $b} keys %{ $adm{$chr}{$key2} }) { # population
	push @subpops, $key3;
	$count1++;
	$subpop{$key3} = $count1;
	}
	last;
	}
	my $count2 = scalar @subpops;


	my %adm_blocks = ();
	for my $i (0..($#aoa-1)) {
	for my $j (1..$#{$aoa[$i]}) {
	if ($aoa[$i+1][$j] != $aoa[$i][$j]) {
	push @{$adm_blocks{$aoa[$i][0]}}, (@{ $aoa[$i] } [ 1..($count2)]) if (!defined $adm_blocks{$aoa[$i][0]});
	push @{$adm_blocks{$aoa[$i+1][0]}}, (@{ $aoa[$i+1] } [ 1..($count2)]) if (!defined $adm_blocks{$aoa[$i+1][0]});
	goto LABEL1;
	}
	}
	LABEL1:
	}


	push @{$adm_blocks{$aoa[0][0]}}, (@{ $aoa[0] } [ 1..($count2)]) if (!defined $adm_blocks{$aoa[0][0]});
	push @{$adm_blocks{$aoa[$#aoa][0]}}, (@{ $aoa[$#aoa] } [ 1..($count2)]) if (!defined $adm_blocks{$aoa[$#aoa][0]});


	my @bp_adm_blocks;
	for my $key ( sort {$a <=> $b} keys %adm_blocks) {
	push @bp_adm_blocks, $key;
	}


	for my $bp (@bp_pos) {
	if (!defined $local_adms{$bp}) {

	if (defined $adm_blocks{$bp}) {
	$local_adms{$bp} = [ @{$adm_blocks{$bp}} ];
	}

	elsif ($bp < $bp_adm_blocks[0] && $extrapolation == 0) { #change if no extrapolation
	$local_adms{$bp} = [ @{$adm_blocks{$bp_adm_blocks[0]}} ];
	}

	elsif ($bp > $bp_adm_blocks[$#bp_adm_blocks] && $extrapolation == 0) { #change if no extrapolation
	$local_adms{$bp} = [ @{$adm_blocks{$bp_adm_blocks[$#bp_adm_blocks]}} ];
	}

	}
	}


	for my $bp (@bp_pos) {
	if (!defined $local_adms{$bp}) {
	for (my $i=0; $i<$#bp_adm_blocks; $i++) {
	if (($bp < $bp_adm_blocks[$i+1]) && ($bp > $bp_adm_blocks[$i])) {
	my @array1 = @{$adm_blocks{$bp_adm_blocks[$i]}};
	my @array2 = @{$adm_blocks{$bp_adm_blocks[$i+1]}};

	my $array1_ref = \@array1;
	my $array2_ref = \@array2;
	my $result1_ref = array_comp ($array1_ref, $array2_ref);
	my $result1 = $$result1_ref;

	if ($result1 == 0) {
	$local_adms{$bp} = [ @{$adm_blocks{$bp_adm_blocks[$i]}} ];
	}

	elsif ($result1 == 1) {
	if ($bp - $bp_adm_blocks[$i] < $bp_adm_blocks[$i+1] - $bp) {
	$local_adms{$bp} = [ @{$adm_blocks{$bp_adm_blocks[$i]}} ];
	}

	elsif ($bp - $bp_adm_blocks[$i] >= $bp_adm_blocks[$i+1] - $bp) { 
	$local_adms{$bp} = [ @{$adm_blocks{$bp_adm_blocks[$i+1]}} ];
	}

	}


	}
	}
	}
	}

	return \(%local_adms, %subpop);
}


sub array_comp { # element-wise comparison of two arrays; input: two arrays to compare, output: 0 = arrays are equal, 1 = array elements are differentat at least at one position

	my ($array_1_ref, $array_2_ref)=@_;
	my @array_1 = @$array_1_ref;
	my @array_2 = @$array_2_ref;
	my $result = 0;

	for (my $i=0; $i<=$#array_1; $i++) {
	if (defined $array_1[$i] && defined $array_2[$i]) {
	if ($array_1[$i] == $array_2[$i]) {
	#do nothing
	} else {
	$result = 1;
	}
	}
	}

	return \$result;
}
### end of section for computation of admixture proportions



### section for computation of admixture model specific allele frequencies
sub frq_template_geno {

	my ($parameters_ref, $pop_frqs_ref, $rs_to_bp_pop_ref, $bp_rs_pos_for_frq_ref, $chr_num_ref) = @_;
	my %parameters = %$parameters_ref;
	my %pop_frqs = %$pop_frqs_ref;
	my %rs_to_bp_pop = %$rs_to_bp_pop_ref;
	my %bp_rs_pos_for_frq = %$bp_rs_pos_for_frq_ref;
	my $chr_num = $$chr_num_ref;

 	my %subpop = ();
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'REFERENCEPOPULATIONNAME'}) {
	$subpop{$parameters{$key1}{'REFERENCEPOPULATIONNAME'}{1}} = 1;
	}
	}

	my @subpop_frq_template; # just to export the order of populations
	for my $sbp (sort {$a cmp $b} keys %subpop) {
	push @subpop_frq_template, $sbp;
	}

	my $pop_num;
	$pop_num = scalar @subpop_frq_template;

	my %frq_template = ();

	for my $key1 (sort {$a <=> $b} keys %bp_rs_pos_for_frq) { # bp position in %bp_rs_pos_for_frq hash
	if (defined $rs_to_bp_pop{$chr_num}{$bp_rs_pos_for_frq{$key1}}) {
	for my $sbp (sort {$a cmp $b} keys %subpop) { # population
	push @{$frq_template{$chr_num}{$key1}{1}}, $pop_frqs{$chr_num}{$rs_to_bp_pop{$chr_num}{$bp_rs_pos_for_frq{$key1}}}{$sbp}{'REF_ALLELE'}, $pop_frqs{$chr_num}{$rs_to_bp_pop{$chr_num}{$bp_rs_pos_for_frq{$key1}}}{$sbp}{'REF_ALLELE_FRQ'};
	$frq_template{$chr_num}{$key1}{2} = $pop_frqs{$chr_num}{$rs_to_bp_pop{$chr_num}{$bp_rs_pos_for_frq{$key1}}}{$sbp}{'ALT_ALLELE'};
	}
	} else {
	print LOGFILE2 "GENO file: rs number $bp_rs_pos_for_frq{$key1} for variant at base pair position $key1 on chromosome $chr_num from the geno file was not matched to the corresponding rs number in reference population dataset.\n";
	}
	}


	for my $key1 (sort {$a <=> $b} keys %bp_rs_pos_for_frq) { # bp position in %bp_rs_pos_for_frq hash
	if (!defined $frq_template{$chr_num}{$key1} && defined $pop_frqs{$chr_num}{$key1}) {
	for my $sbp (sort {$a cmp $b} keys %subpop) { # population
	push @{$frq_template{$chr_num}{$key1}{1}}, $pop_frqs{$chr_num}{$key1}{$sbp}{'REF_ALLELE'}, $pop_frqs{$chr_num}{$key1}{$sbp}{'REF_ALLELE_FRQ'};
	$frq_template{$chr_num}{$key1}{2} = $pop_frqs{$chr_num}{$key1}{$sbp}{'ALT_ALLELE'};
	}
	}

	elsif (!defined $frq_template{$chr_num}{$key1} && !defined $pop_frqs{$chr_num}{$key1}) {
	print LOGFILE2 "GENO file: a geno file variant at base pair position $key1 on chromosome $chr_num cannot be matched to reference population dataset neither by its rs number nor base pair position.\n";
	}

	}


	for my $key1 (sort {$a <=> $b} keys %bp_rs_pos_for_frq) { # bp position in %bp_rs_pos_for_frq hash
	if (defined $frq_template{$chr_num}{$key1}) {

	my @alleles;
	my @frqs;
	for (my $el1=1; $el1<=$pop_num; $el1++) {
	my $temp_al = ($el1 * 2) - 2;
	my $temp_frq = ($el1 * 2) - 1;
	push @alleles, $frq_template{$chr_num}{$key1}{1}[$temp_al] if (exists $frq_template{$chr_num}{$key1}{1}[$temp_al]);
	push @frqs, $frq_template{$chr_num}{$key1}{1}[$temp_frq] if (exists $frq_template{$chr_num}{$key1}{1}[$temp_frq]);
	}

	my %check_al = map { $_, 1 } @alleles; # checks that all reference alleles are the same
	if (keys %check_al != 1 || scalar @alleles != 3) {
	error ("Error: not all reference alleles across specified reference populations for the variant on chromosome $chr_num at $key1 base pair position are the same and/or frequency data for the above mentioned variant are not available across all reference populations! A panel of markers in the 'geno' file must contain variants with MAF>0.25 and frequency data available across all reference populations!");
	}

	if (scalar @frqs != 3) {
	error ("Error: allele frequency data for the variant on $chr_num at $key1 base pair position are not available across specified reference populations!");
	}


	} else {
error ("Error: there are no reference population allele frequency data for a marker at $key1 base pair position!");
	}
	}

	return \%frq_template, \@subpop_frq_template;
}


sub frq_template {

	my ($parameters_ref, $sample_frqs_ref, $pop_frqs_ref, $bp_pos_for_frq_ref, $chr_num_ref, $sample_bp_rs_ref) = @_;
	my %parameters = %$parameters_ref;
	my %sample_frqs = %$sample_frqs_ref;
	my %pop_frqs = %$pop_frqs_ref;
	my @bp_pos_for_frq = @$bp_pos_for_frq_ref;
	my $chr_num = $$chr_num_ref;
	my %sample_bp_rs = %$sample_bp_rs_ref if ($sample_bp_rs_ref != 7);

	my $rs_check = 0;
	$rs_check = 1 if ($sample_bp_rs_ref != 7);
 	my %subpop = ();
	for my $key1 (sort keys %parameters) {
	if (defined $parameters{$key1}{'REFERENCEPOPULATIONNAME'}) {
	$subpop{$parameters{$key1}{'REFERENCEPOPULATIONNAME'}{1}} = 1;
	}
	}

	my @subpop_frq_template; # just to export the order of populations
	for my $sbp (sort {$a cmp $b} keys %subpop) {
	push @subpop_frq_template, $sbp;
	}

	my %frq_template = ();

	for my $el1 (sort @bp_pos_for_frq) {
	if (defined $sample_frqs{$chr_num}{$el1}) {
	for my $sbp (sort {$a cmp $b} keys %subpop) {

	if (defined $pop_frqs{$chr_num}{$el1}{$sbp} && defined $frq_template{$chr_num}{$el1}{1}) {

	if ($rs_check == 1 && defined $sample_bp_rs{$chr_num}{$el1} && defined $pop_frqs{$chr_num}{$el1}{$sbp}{'RS_NUMBER'} && (lc($sample_bp_rs{$chr_num}{$el1}) ne 'na')) {
	if ($sample_bp_rs{$chr_num}{$el1} eq $pop_frqs{$chr_num}{$el1}{$sbp}{'RS_NUMBER'}) {
# do nothing
	} else {
	print LOGFILE2 "rs number $sample_bp_rs{$chr_num}{$el1} for variant at $el1 base pair position on chromosome $chr_num was not matched to the rs number at the same base pair position in $sbp reference population dataset because either your sample and reference population datasets are based on different builds or different patches of the same build.\n";
	}
	}

	my $alt_frq = 1 - $pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE_FRQ'};
	if ($pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE'} eq @{$frq_template{$chr_num}{$el1}{1}}[0]) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE'};
	push @{$frq_template{$chr_num}{$el1}{1}}, $pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE_FRQ'};
	}
	elsif ($pop_frqs{$chr_num}{$el1}{$sbp}{'ALT_ALLELE'} eq @{$frq_template{$chr_num}{$el1}{1}}[0]) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $pop_frqs{$chr_num}{$el1}{$sbp}{'ALT_ALLELE'};
	push @{$frq_template{$chr_num}{$el1}{1}}, $alt_frq;
	} else {
	error ("Error: reference and alternative alleles for the variant at $el1 base pair position on chromosome $chr_num in $sbp reference population are different from those observed in your sample!");
	}
	}


	elsif (defined $frq_template{$chr_num}{$el1}{1} && !defined $pop_frqs{$chr_num}{$el1}{$sbp}) {
	if ($sample_frqs{$chr_num}{$el1}{'al1'} eq @{$frq_template{$chr_num}{$el1}{1}}[0] && $sample_frqs{$chr_num}{$el1}{'al1_frq'} >= 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al1'}, 0.99; 
	}
	elsif ($sample_frqs{$chr_num}{$el1}{'al1'} eq @{$frq_template{$chr_num}{$el1}{1}}[0] && $sample_frqs{$chr_num}{$el1}{'al1_frq'} < 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al1'}, 0.01; 
	}
	elsif ($sample_frqs{$chr_num}{$el1}{'al2'} eq @{$frq_template{$chr_num}{$el1}{1}}[0] && $sample_frqs{$chr_num}{$el1}{'al2_frq'} >= 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al2'}, 0.99; 
	}
	elsif ($sample_frqs{$chr_num}{$el1}{'al2'} eq @{$frq_template{$chr_num}{$el1}{1}}[0] && $sample_frqs{$chr_num}{$el1}{'al2_frq'} < 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al2'}, 0.01; 
	} else {
	error ("Error: reference and alternative alleles for the variant at $el1 base pair position on chromosome $chr_num in $sbp reference population are different from those observed in your sample!");
	}
	}


	elsif (defined $pop_frqs{$chr_num}{$el1}{$sbp} && !defined $frq_template{$chr_num}{$el1}{1}) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE'};
	push @{$frq_template{$chr_num}{$el1}{1}}, $pop_frqs{$chr_num}{$el1}{$sbp}{'REF_ALLELE_FRQ'};
	$frq_template{$chr_num}{$el1}{2} = $pop_frqs{$chr_num}{$el1}{$sbp}{'ALT_ALLELE'};
	}


	elsif (!defined $pop_frqs{$chr_num}{$el1}{$sbp} && !defined $frq_template{$chr_num}{$el1}{1}) {
	if ($sample_frqs{$chr_num}{$el1}{'al1_frq'} >= 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al1'}, 0.99;
	$frq_template{$chr_num}{$el1}{2} = $sample_frqs{$chr_num}{$el1}{'al2'};
	}
	elsif ($sample_frqs{$chr_num}{$el1}{'al2_frq'} > 0.5) {
	push @{$frq_template{$chr_num}{$el1}{1}}, $sample_frqs{$chr_num}{$el1}{'al2'}, 0.99;
	$frq_template{$chr_num}{$el1}{2} = $sample_frqs{$chr_num}{$el1}{'al1'};
	}
	}


	}


	}
	}


return \%frq_template, \@subpop_frq_template;
}


sub global_adm_frqs {

	my ($frq_template_ref, $chr_num_ref, $subpop_frq_template_ref, $global_adm_ref) = @_;
	my @subpop_frq_template = @$subpop_frq_template_ref;
	my $chr_num = $$chr_num_ref;
	my %frq_template = %$frq_template_ref;
	my %global_adm = %$global_adm_ref;

	my %global_frqs = ();
	my $contr_flow = 0;

	for my $key1 (sort keys %{ $frq_template{$chr_num} }) { # key1 = bp position
	push @{$global_frqs{$chr_num}{$key1}{1}}, $frq_template{$chr_num}{$key1}{1}[0];
	push @{$global_frqs{$chr_num}{$key1}{2}}, $frq_template{$chr_num}{$key1}{2};

	my @temp1;
	for (my $i=1; $i<=$#subpop_frq_template+1; $i++) {
	my $temp2 = ($i-1) + $i;
	if (defined $frq_template{$chr_num}{$key1}{1}[$temp2] && defined $global_adm{$subpop_frq_template[$i-1]}) {
	my $temp3 = $frq_template{$chr_num}{$key1}{1}[$temp2] * $global_adm{$subpop_frq_template[$i-1]};
	push @temp1, $temp3;
	} else {
	warning ("Frequency at $key1 base pair position on chromosome $chr_num and for $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $frq_template{$chr_num}{$key1}{1}[$temp2];
	warning ("Global admixture proportion for $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $global_adm{$subpop_frq_template[$i-1]};
	$contr_flow = 1;
	goto LABEL1;
	}
	}
	my $temp4;
	for (@temp1) {
	$temp4 += $_;
	}
	$temp4 = sprintf ("%.6f", $temp4);
	push @{$global_frqs{$chr_num}{$key1}{1}}, $temp4 if ($temp4 < 1);
	push @{$global_frqs{$chr_num}{$key1}{1}}, 0.999000 if ($temp4 >= 1 && $temp4 <= 1.02);
	warning ("MAF for variant at base pair position $key1 on chromosome $chr_num is > 1.02. Something might be wrong with either admixture proportions or reference allele frequencies for this variant") if ($temp4 > 1.02);

	if ($global_frqs{$chr_num}{$key1}{1}[1] == 0.999000) {
	push @{$global_frqs{$chr_num}{$key1}{2}}, 0.001000;
	} else {
	my $temp5 = 1 - $temp4;
	$temp5 = sprintf ("%.6f", $temp5) ;
	push @{$global_frqs{$chr_num}{$key1}{2}}, $temp5;
	}

	LABEL1:
	if ($contr_flow == 1) {
	push @{$global_frqs{$chr_num}{$key1}{1}}, 'NA';
	push @{$global_frqs{$chr_num}{$key1}{2}}, 'NA';
	$contr_flow = 0;
	}

	}


	return \%global_frqs;
}


sub fbgw_frqs {

	my ($frq_template_ref, $subpop_frq_template_ref, $fbgw_ref, $family_ref, $chr_num_ref) = @_;

	my %frq_template = %$frq_template_ref;
	my @subpop_frq_template = @$subpop_frq_template_ref;
	my %fbgw = %$fbgw_ref;
	my $family = $$family_ref;
	my $chr_num = $$chr_num_ref;

	my %fbgw_frqs = ();
	my $contr_flow = 0;

	for my $key1 (sort keys %{ $frq_template{$chr_num} }) { # key1 = bp position
	push @{$fbgw_frqs{$chr_num}{$key1}{1}}, $frq_template{$chr_num}{$key1}{1}[0];
	push @{$fbgw_frqs{$chr_num}{$key1}{2}}, $frq_template{$chr_num}{$key1}{2};
	my @temp1;
	for (my $i=1; $i<=$#subpop_frq_template+1; $i++) {
	my $temp2 = ($i-1) + $i;
	if (defined $frq_template{$chr_num}{$key1}{1}[$temp2] && defined $fbgw{$family}{$subpop_frq_template[$i-1]}) {
	my $temp3 = $frq_template{$chr_num}{$key1}{1}[$temp2] * $fbgw{$family}{$subpop_frq_template[$i-1]};
	push @temp1, $temp3;
	} else {
	warning ("Frequency at $key1 base pair position on chromosome $chr_num and for $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $frq_template{$chr_num}{$key1}{1}[$temp2];
	warning ("FBGW admixture proportion for family $family and $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $fbgw{$family}{$subpop_frq_template[$i-1]};
	$contr_flow = 1;
	goto LABEL1;
	}
	}
	my $temp4;
	for (@temp1) {
	$temp4 += $_;
	}
	$temp4 = sprintf ("%.6f", $temp4);
	push @{$fbgw_frqs{$chr_num}{$key1}{1}}, $temp4 if ($temp4 < 1);
	push @{$fbgw_frqs{$chr_num}{$key1}{1}}, 0.999000 if ($temp4 >= 1 && $temp4 <= 1.02);
	warning ("MAF for variant at base pair position $key1 on chromosome $chr_num is > 1.02. Something might be wrong with either admixture proportions or reference allele frequencies for this variant") if ($temp4 > 1.02);

	if ($fbgw_frqs{$chr_num}{$key1}{1}[1] == 0.999000) {
	push @{$fbgw_frqs{$chr_num}{$key1}{2}}, 0.001000;
	} else {
	my $temp5 = 1 - $temp4;
	$temp5 = sprintf ("%.6f", $temp5) ;
	push @{$fbgw_frqs{$chr_num}{$key1}{2}}, $temp5;
	}

	LABEL1:
	if ($contr_flow == 1) {
	push @{$fbgw_frqs{$chr_num}{$key1}{1}}, 'NA';
	push @{$fbgw_frqs{$chr_num}{$key1}{2}}, 'NA';
	$contr_flow = 0;
	}

	}

	return \%fbgw_frqs;
}


sub fbcw_frqs {

	my ($frq_template_ref, $subpop_frq_template_ref, $fbcw_ref,  $family_ref, $chr_num_ref) = @_;

	my %frq_template = %$frq_template_ref;
	my @subpop_frq_template = @$subpop_frq_template_ref;
	my %fbcw = %$fbcw_ref;
	my $family = $$family_ref;
	my $chr_num = $$chr_num_ref;

	my %fbcw_frqs = ();
	my $contr_flow = 0;

	for my $key1 (sort keys %{ $frq_template{$chr_num} }) { # key1 = bp position
	push @{$fbcw_frqs{$chr_num}{$key1}{1}}, $frq_template{$chr_num}{$key1}{1}[0];
	push @{$fbcw_frqs{$chr_num}{$key1}{2}}, $frq_template{$chr_num}{$key1}{2};
	my @temp1;
	for (my $i=1; $i<=$#subpop_frq_template+1; $i++) {
	my $temp2 = ($i-1) + $i;
	if (defined $frq_template{$chr_num}{$key1}{1}[$temp2] && defined $fbcw{$chr_num}{$family}{$subpop_frq_template[$i-1]}) {
	my $temp3 = $frq_template{$chr_num}{$key1}{1}[$temp2] * $fbcw{$chr_num}{$family}{$subpop_frq_template[$i-1]};
	push @temp1, $temp3;
	} else {
	warning ("Frequency at $key1 base pair position on chromosome $chr_num and for $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $frq_template{$chr_num}{$key1}{1}[$temp2] ;
	warning ("FBCW admixture proportion for family $family, chromosome $chr_num, and $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $fbcw{$chr_num}{$family}{$subpop_frq_template[$i-1]};
	$contr_flow = 1;
	goto LABEL1;
	}
	}
	my $temp4;
	for (@temp1) {
	$temp4 += $_;
	}
	$temp4 = sprintf ("%.6f", $temp4);
	push @{$fbcw_frqs{$chr_num}{$key1}{1}}, $temp4 if ($temp4 < 1);
	push @{$fbcw_frqs{$chr_num}{$key1}{1}}, 0.999000 if ($temp4 >= 1 && $temp4 <= 1.02);
	warning ("MAF for variant at base pair position $key1 on chromosome $chr_num is > 1.02. Something might be wrong with either admixture proportions or reference allele frequencies for this variant") if ($temp4 > 1.02);

	if ($fbcw_frqs{$chr_num}{$key1}{1}[1] == 0.999000) {
	push @{$fbcw_frqs{$chr_num}{$key1}{2}}, 0.001000;
	} else {
	my $temp5 = 1 - $temp4;
	$temp5 = sprintf ("%.6f", $temp5) ;
	push @{$fbcw_frqs{$chr_num}{$key1}{2}}, $temp5;
	}

	LABEL1:
	if ($contr_flow == 1) {
	push @{$fbcw_frqs{$chr_num}{$key1}{1}}, 'NA';
	push @{$fbcw_frqs{$chr_num}{$key1}{2}}, 'NA';
	$contr_flow = 0;
	}

	}

	return \%fbcw_frqs;
}


sub local_frqs {

	my ($frq_template_ref, $subpop_frq_template_ref, $local_adms_ref,  $family_ref, $chr_num_ref, $subpop_local_ref, $geno_adm_l) = @_; # $geno_adm_l = 0 (no Local adm for geno files) or 1 (Local adm for geno files)

	my %frq_template = %$frq_template_ref;
	my %subpop_local = %$subpop_local_ref;
	my @subpop_frq_template = @$subpop_frq_template_ref;
	my %local_adms = %$local_adms_ref;
	my $family = $$family_ref;
	my $chr_num = $$chr_num_ref;

	my %local_frqs = ();
	my %geno_adm = ();
	my $contr_flow = 0;

	for my $key1 (sort {$a <=> $b} keys %{ $frq_template{$chr_num} }) { # key1 = bp position
	push @{$local_frqs{$chr_num}{$key1}{1}}, $frq_template{$chr_num}{$key1}{1}[0];
	push @{$local_frqs{$chr_num}{$key1}{2}}, $frq_template{$chr_num}{$key1}{2};
	my @temp1;
	for (my $i=1; $i<=$#subpop_frq_template+1; $i++) {
	my $temp2 = ($i-1) + $i;
	if (defined $frq_template{$chr_num}{$key1}{1}[$temp2] && defined $local_adms{$key1}[$subpop_local{$subpop_frq_template[$i-1]}-1]) {
	my $temp3 = $frq_template{$chr_num}{$key1}{1}[$temp2] * $local_adms{$key1}[$subpop_local{$subpop_frq_template[$i-1]}-1];
	$geno_adm{$chr_num}{$key1}{$subpop_frq_template[$i-1]} = sprintf ("%.6f", $local_adms{$key1}[$subpop_local{$subpop_frq_template[$i-1]}-1]) if ($geno_adm_l == 1);
	push @temp1, $temp3;
	} else {
	warning ("Frequency at $key1 base pair position on chromosome $chr_num and for $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $frq_template{$chr_num}{$key1}{1}[$temp2];
	warning ("Local admixture proportion at $key1 base pair position on chromosome $chr_num for family $family and $subpop_frq_template[$i-1] reference population is not defined.\n") if !defined $local_adms{$key1}[$subpop_local{$subpop_frq_template[$i-1]}-1];
	$contr_flow = 1;
	goto LABEL1;
	}
	}
	my $temp4;
	for (@temp1) {
	$temp4 += $_;
	}
	$temp4 = sprintf ("%.6f", $temp4);
	push @{$local_frqs{$chr_num}{$key1}{1}}, $temp4 if ($temp4 < 1);
	push @{$local_frqs{$chr_num}{$key1}{1}}, 0.999000 if ($temp4 >= 1 && $temp4 <= 1.02);
	warning ("MAF for variant at base pair position $key1 on chromosome $chr_num is > 1.02. Something might be wrong with either admixture proportions or reference allele frequencies for this variant") if ($temp4 > 1.02);

	if ($local_frqs{$chr_num}{$key1}{1}[1] == 0.999000) {
	push @{$local_frqs{$chr_num}{$key1}{2}}, 0.001000;
	} else {
	my $temp5 = 1 - $temp4;
	$temp5 = sprintf ("%.6f", $temp5) ;
	push @{$local_frqs{$chr_num}{$key1}{2}}, $temp5;
	}

	LABEL1:
	if ($contr_flow == 1) {
	push @{$local_frqs{$chr_num}{$key1}{1}}, 'NA';
	push @{$local_frqs{$chr_num}{$key1}{2}}, 'NA';
	$contr_flow = 0;
	}

	}

	return \%local_frqs if ($geno_adm_l == 0);
	return \%local_frqs, \%geno_adm if ($geno_adm_l == 1);
}
### end of section for computation of admixture model specific allele frequencies



sub script_info {

	my $parameters_ref = $_[0];

	my $comp_type = 4;
	my $comp_type_ref = \$comp_type;
	my $dir_ref = make_dir ($parameters_ref, 0, 0, 0, $comp_type_ref);
	my $dir = $$dir_ref;

	my $date_time_2 = localtime;
	$date_time_2 =~ s/\s+/./g;
	$date_time_2 =~ s/\:+/./g;
	my $rand_n = sprintf ("%.0f", rand (999));
	open (LOGFILE, ">", "$dir/admixfrq.$date_time_2.r$rand_n.log");
	open (LOGFILE2, ">", "$dir/rs_num_missmatch.$date_time_2.r$rand_n.log");
	my $date_time = localtime;
	print "\n",
	"  \@----------------------------------------------------------\@\n",
	 "  |    ADMIXFRQ!  |     v1.00      |   14/May/2018           |\n",
	 "  |----------------------------------------------------------|\n",
	 "  |Copyright (C) 2018  Rafael A. Nafikov  GNU GPLv3          |\n",
	 "  |This program comes with ABSOLUTELY NO WARRANTY;           |\n",
	 "  |for details see <https://www.gnu.org/licenses/>.          |\n",
	 "  |This is free software, and you are welcome to redistribute|\n",
	 "  |it under certain conditions; for details see              |\n",
	 "  |<https://www.gnu.org/licenses/>.                          |\n",
	 "  |----------------------------------------------------------|\n",
	 "  | For documentation, citation & bug-reporting instructions:|\n",
	 "  |        https://github.com/RafPrograms/ADMIXFRQ           |\n",
	"  \@----------------------------------------------------------\@\n",
	 "\n",
	 "	Analysis started at: $date_time\n\n";

	print LOGFILE "\n",
	"  \@----------------------------------------------------------\@\n",
	 "  |    ADMIXFRQ!  |     v1.00      |   14/May/2018           |\n",
	 "  |----------------------------------------------------------|\n",
	 "  |Copyright (C) 2018  Rafael A. Nafikov  GNU GPLv3          |\n",
	 "  |This program comes with ABSOLUTELY NO WARRANTY;           |\n",
	 "  |for details see <https://www.gnu.org/licenses/>.          |\n",
	 "  |This is free software, and you are welcome to redistribute|\n",
	 "  |it under certain conditions; for details see              |\n",
	 "  |<https://www.gnu.org/licenses/>.                          |\n",
	 "  |----------------------------------------------------------|\n",
	 "  | For documentation, citation & bug-reporting instructions:|\n",
	 "  |        https://github.com/RafPrograms/ADMIXFRQ           |\n",
	"  \@----------------------------------------------------------\@\n",
	 "\n",
	 "	Analysis started at: $date_time\n\n";

	return *LOGFILE, *LOGFILE2;
}

