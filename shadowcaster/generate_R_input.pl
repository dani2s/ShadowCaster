#!/usr/bin/perl -w

############################################################################
#    This file is part of ShadowCaster.
#    Copyright (C) 2018  Daniela Sanchez and Aminael Sanchez
#
#    ShadowCaster is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ShadowCaster is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

use strict;
use Bio::SeqIO;
use Bio::SearchIO;

use File::Spec::Functions 'catfile';


my $input_file = $ARGV[0];
my $output_file = "alien_svm.csv";
my $blastFolder = $ARGV[1];
my $soft = $ARGV[2];
my $protPath = $ARGV[3];

opendir(DIR, $protPath);
my @species = grep(/\.(fasta|faa)$/,readdir(DIR));
closedir(DIR);

open (OUTPUT, ">".$output_file) || die "Can't open $output_file, \n";
my @vector = ();
my $frac_id;
my $qlength;
my $seqIOobj = Bio::SeqIO->new(-file=> $input_file);
while((my $seqobj = $seqIOobj->next_seq())) {
	my $id = $seqobj->id;
	$id =~ s/\.\d//;
	my $tmp_fasta = $blastFolder.$id.".fasta";
	open (FASTA, ">".$tmp_fasta) || die "Can't open $tmp_fasta";
	print FASTA ">".$id."\n".$seqobj->seq."\n";
	close FASTA;
	@vector = ();
	push(@vector, $id);
	foreach my $sp (@species){
		#print "Blasting ".$id." versus ".$sp."...\n";
		my $spDatabase = catfile($protPath, $sp);	
				
		system ($soft." -query ".$tmp_fasta." -db ".$spDatabase." -out ".$blastFolder.$id."_vs_".$sp.".blastout");
		
		#print "DONE! Output stored at ".$id."_vs_".$sp.".blastout\n";
		my $blastReport = Bio::SearchIO->new(-format => 'blast',
		                                     -file => $blastFolder.$id."_vs_".$sp.".blastout");
		my $found = 0;
		while (my $result = $blastReport->next_result){
			$qlength = $result->query_length;
			while (my $hit = $result->next_hit){
				$found++;
				my $hit_name = $hit->name();
				my $len = $hit->length();
				my $frac_query_ali = $hit->frac_aligned_query();
				while (my $hsp = $hit->next_hsp){
					$frac_id = $hsp->frac_identical('total');
					if ($frac_query_ali < 0.80){
						$frac_id = 0;
					}
					#print OUTPUT "Escherichia_coli\t".$id."\t".$qlength."\t".$frac_query_ali."\t".$frac_id."\t".$sp."\t".$hit_name."\n";
					last;
				}
				last;
			}
		}
		if ($found == 0){
			$frac_id = 0;
			#print OUTPUT "Escherichia_coli\t".$id."\t0\t0\t0\t".$sp."\tNONE\n";
		}
		push(@vector, $frac_id)
	}
	if ($qlength >= 100){
		print OUTPUT join("\t",@vector)."\n";
	}
}

close OUTPUT;