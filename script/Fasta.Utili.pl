#!/usr/bin/env perl 

###################################################
### This script iclude the useful command function# 
### for handle the fasta file 					  #
###################################################

use strict;
use warnings;
use utf8;
use 5.010;
use Data::Dumper;
use Bio::SeqIO;
use Getopt::Long;
use IO::File;
use Switch;
use List::Util qw/sum min max/;
use List::MoreUtils qw/uniq/;
use lib "/home/wudi/lib";
require "ABSyS.pl";

my ($acc,@file,$moulde,$cutofflenth,$outfile,$HELP);

GetOptions(
	"i:s{,}" =>\@file,
	"m:s" =>\$moulde,
	"cut:s" => \$cutofflenth,
	"o:s" =>\$outfile,
	"acc:s" => \$acc,
	"h|?"=>sub {
	say "\nUSAGE:";
	say "Fasta.utili.pl -m <method> -acc <gene acc | regex of acc | gff file > -i <input file> -cut <len cutoff> -o <output file> -h <show this USAGE>";
	say "   include those method:";
	say "     mean : FASTA file mean len";
	say "     max : FASTA file max len";
	say "     min : FASTA file min len";
	say "     N50 : FASTA file N50 value";
	say "     N90 : FASTA file N90 value";
	say "     GC : GC content percentage in the Sequences";
	say "     GC_Nfree : GC content percentage in the Sequences (without N)";
	say "     seqcount : count of seqences in FASTA file";
	say "     total_len : FASTA file total length";
	say "     filter : FASTA file length cutoff (need output file specified)";
	say "     muti2uniq : dump several fasta files into one (need output file specified)";
	say "     do_meta_Scf_report : Print Meta-Genomics Scaffold summary ";
	say "     do_meta_CDS_report : Print Meta-Genomics CDS summary";
	say "     cutfasta : cut fasta file into small ones ";
	say "     delete_gap : delete the '-' in FASTA file";
	say "     grep : grep seq_id with regular expression";
	say "     meta_scaffold_title : print meta scaffold title ";
	say "     Dump_Circos_Karyotype : dump the karyotype file for Circos ";
	say "     Fetch_Seqs : Fetch sequences with accessions split with \",\"";
	say "     Filter_Seqs : Filter sequences with accessions file";
	say "     make_chrom_file : make chromesome file for bedtools ";
	say "     make_length_distribution : give any sequence in fasta file length"; 
	say "     max_seq : print a max seq of ONE fasta file"; 
	say "     meta_CDS_title : print meta CDS title ";
	say "     cut_inlen : cut fasta into same length";
	say "         Exmaple : perl $0 -m cut_inlen -i <fa file> -cut 700 -o <out fasta file >\n";
	say "     cut_inrange : cut fasta with range ";
	say "         Exmaple : perl $0 -m cut_inlen -i <fa file> -cut < range file > -o <out fasta file >\n";
	say "     r_cut_inrange : cut fasta with rever range ";
	say "         Exmaple : perl $0 -m r_cut_inlen -i <fa file> -cut < range file > -o <out fasta file >\n";
	say "     amphi_n : head and tail N number ";
	say "     pretty_print : as name says... ";
	say "     gap2bed : as name says... ";

	say "     rename : rename the any sequence title with adding prefix \n\n";
	exit(0);
	}
);



switch($moulde){
case "mean"								 	{say mean_len(@file)}
case "max"									{say max_len(@file)}
case "min"									{say min_len(@file)}
case "seqcount" 							{print `grep -c ">" @file`}
case "N50"									{say N50(@file)}
case "N90"									{say N90(@file)}
case "total_len"							{say Total_len(@file)}
case "do_meta_Scf_report"					{say do_meta_Scf_report(@file)}
case "do_meta_CDS_report"					{say do_meta_CDS_report(@file)}
case "filter"								{filter(@file,$outfile)}
case "muti2uniq"							{mutiple2uniq($outfile,@file)}
case "GC"									{GC(@file)}
case "GC_Nfree"								{GC_Nfree(@file)}
case "cutfasta"								{cut_fasta_file(@file,$outfile,$cutofflenth)}
case "meta_scaffold_title"					{meta_scaffold_title()}
case "meta_CDS_title"       				{meta_CDS_title()}
case "delete_gap"							{delete_gap(@file,$outfile)}
case "grep"									{grep_by_name(@file,$outfile,$acc)}
case "make_chrom_file"                      {make_chrom_file($outfile ,@file)}
case "pos"									{posi($acc,	@file)}
case "Dump_Circos_Karyotype"				{Circos_dump($acc,@file)}
case "Fetch_Seqs"							{Fetch_Seq($acc,$outfile,@file)}
case "Filter_Seqs"							{Filter_Seq($acc,$outfile,@file)}
case "make_length_distribution"				{Distribution(@file)}
case "max_seq"								{say max_seq(@file)}
case "rename"								{ReName($acc,@file )}
case "cut_inlen"							{cut_inlen($cutofflenth,$outfile,@file)}
case "cut_inrange"							{cut_inrange($cutofflenth,$outfile,@file)}
case "r_cut_inrange"						{revers_cut_inrange($cutofflenth,$outfile,@file)}
case "amphi_n"								{STATIC_AMPHI_N(@file)}
case "pretty_print"							{pretty_print(@file)}
case "gap2bed"								{gapTObed(@file)}
}

sub gapTObed{
	my $file = shift;
	my ($in,$out);
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
 	 if($outfile){
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	 }else{
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDOUT);
	 }
	 while(my $seqobj = $in->next_seq() ){
		my $SEQ = $seqobj->seq();
		my ($result,$flag,$j) ;
		$flag = 0 ; $j = 0;
		for(my $i = 1 ; $i <= length($SEQ) ; $i++){
			if($flag == 0){
				if( substr($SEQ,$i-1,1) eq "N" ){
					$result->[$j]->{"start"} = $i ;
					$flag = 1; 
				}
			}else{
				 if(substr($SEQ,$i-1,1) eq "N"){
					 next;
				 }else{
					 $result->[$j]->{"end"} = $i -1 ;
					 $j ++ ;
					 $flag = 0 ;
				 }			 
			}
		}
		foreach(@{$result}){
			say join "\t",($seqobj->primary_id,$_->{"start"},$_->{"end"} )
			#say Dumper $_;
		}

	}
	exit(0);

}

sub pretty_print{
	my $file = shift;
	my $outfile = shift;

	my ($in,$out);
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
 	 if($outfile){
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	 }else{
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDOUT);
	 }
	 while(my $seqobj = $in->next_seq() ){
			$out->write_seq($seqobj)	;
	}
	exit(0);
}

sub revers_cut_inrange{
	my $cut = shift;
	my $outpath = shift;
	my $file = shift;
	open OUT,">$outpath";	
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 my $range;
	 open CUT,$cut || die "Can't open range file!\n";
   	while(<CUT>){
		chomp;
		my @a = split "\t";
		my @b = split ":",$a[1];
		$range->{$a[0]}->{start} = $b[0];
		$range->{$a[0]}->{end} = $b[1];
	}
	close CUT;

	 while(my $seqobj = $in->next_seq() ){
		if(not $range->{$seqobj->primary_id} ) {
			say OUT ">",$seqobj->primary_id;
			say OUT $seqobj->seq;
		}else{
			my $name = $range->{$seqobj->primary_id};
			my $s = $seqobj->seq;
			
			my @sub_s = Mycut($s,$name->{start} -1 ,$name->{end} -1); 
			my $i = 1;
#			substr($s , 
#			$name->{start} -1 ,
#			$name->{end} - $name->{start} +1 ) = "";
#			if(length() == 0){next}
			foreach(@sub_s){
				say OUT ">",$seqobj->primary_id," ",$name->{start},":",$name->{end},"_$i";
				say OUT $_;
				$i ++;
			}
		}	
	}
	close OUT;
	exit(0)

}

sub cut_inrange{
	my $cut = shift;
	my $outpath = shift;
	my $file = shift;
	open OUT,">$outpath";	
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 my $range;
	 open CUT,$cut || die "Can't open range file!\n";

   	while(<CUT>){
		chomp;
		my @a = split "\t";
		my @b = split ":",$a[1];
		$range->{$a[0]}->{start} = $b[0];
		$range->{$a[0]}->{end} = $b[1];
	}

	close CUT;
	my $i = 0;
	 while(my $seqobj = $in->next_seq() ){
		if(not $range->{$seqobj->primary_id} ) {next} 
		my $name = $range->{$seqobj->primary_id};
		say OUT ">",$seqobj->primary_id," ",$name->{start},":",$name->{end},"_$i";
		say OUT substr($seqobj->seq , 
			$name->{start} -1 ,
		$name->{end} - $name->{start} +1 )
	}
	close OUT;
	exit(0)

}

sub Mycut{
	my $line = shift;
	my $start = shift;
	my $end = shift;
	my $ele1 = substr( $line,0,$start ) ;
	my $ele2 = substr($line,$end +1 );
	return ($ele1,$ele2);
}


sub cut_inlen{
	my $cut = shift;
	my $outpath = shift;
	my $file = shift;
	open OUT,">$outpath";	
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 
	 while(my $seqobj = $in->next_seq() ){
		 
		my @array = unpack ("(A$cut)*",$seqobj->seq());
		foreach my $i ( 0 .. $#array){
			my $start = $i*$cut + 1;
			my $end = length($array[$i]) <= $cut ? $i * $cut +  length($array[$i]) :($i +1)*$cut;
			say OUT ">",$seqobj->primary_id."_$i","_".$start.":".$end;
			say OUT $array[$i];
		}
	}
	exit(0)
}

sub ReName{
	my $acc =shift;
	my $file = shift;
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 while(my $seqobj = $in->next_seq() ){
		say ">",$seqobj->primary_id."_".$acc,"\n",$seqobj->seq;
	}
	exit(0)

}

sub max_seq {
	my $file = shift;
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	my $max_len = 0;
	my $max_Seq = 0;
	while(my $seqobj = $in->next_seq() ){
		$seqobj->length() >= $max_len ? ( $max_len = $seqobj->length() and $max_Seq = $seqobj ) : next;
	}
	say ">",$max_Seq->primary_id,"\n",$max_Seq->seq;
	exit(0)
}





sub Distribution{
	my $fasta = shift;
	my ($in,$out);
	if( $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	while(my $seqobj = $in->next_seq() ){
		say $seqobj->primary_id,"\t",$seqobj->length;
	}
	exit(0);
}

sub Filter_Seq{
	my $name = shift;
	my $outfile = shift;
	my $fasta = shift;

	open FH,$name ; 
	my @names = grep {chomp} <FH>;
	my ($in,$out);
	if( $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }

 	 if($outfile){
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	 }else{
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDOUT);
	 }

	 while(my $seqobj = $in->next_seq() ){
		if (not($seqobj->primary_id ~~ @names)){
			$out->write_seq($seqobj)	
		}
	}
	exit(0);
}





sub Fetch_Seq{
	my $name = shift;
	my $outfile = shift;
	my $fasta = shift;
	my @names;
	if(-e $name ) { 	
		open FH,$name ; 
		@names = grep {chomp} <FH>;
	}else{
		@names = split ",",$name;
	}

	my ($in,$out);
	if( $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }

 	 if($outfile){
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	 }else{
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDOUT);
	 }

	 while(my $seqobj = $in->next_seq() ){
		if ($seqobj->primary_id ~~ @names){
			$out->write_seq($seqobj)	
		}
	}
	exit(0);
}

sub Circos_dump {
	my $name = shift;
	my $fasta = shift ;
	my $in;
	my $i = 0;
	if( $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 while(my $seqobj = $in->next_seq() ){
			 say join " ",("chr","-", $seqobj->primary_id,$seqobj->primary_id ,"0",$seqobj->length(),$name.($i+1));
			 $i ++ ;
	}
	exit(0);
}

sub posi{
	my $pos = shift;
	my $fasta = shift;
	my $in;
	if(-e $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 while(my $seqobj = $in->next_seq() ){
		say substr($seqobj->seq,$pos,1)
	}
	return 0;
	exit(0);
}

sub make_chrom_file{
	my $out = shift;
	my $fasta = shift;
	my $in;
	if(-e $fasta){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $fasta);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	 open OUT , ">$out" || die "Can't open out file for write\n";
	 while(my $seqobj = $in->next_seq() ){
		say OUT $seqobj->primary_id ,"\t",$seqobj->length()
	}
	close OUT;
	return 0;
	exit(0);
}

sub get_NCBI_seq{
	my $outfile = shift;
	my @acc = @_;
	my $fh = new IO::File "> $outfile" || die "can't open the output file\n";
	foreach(@acc){
		my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/".
		"eutils/efetch.fcgi?db=nuccore&id=$_&rettype=fasta"	;
		my $fasta  = get($url);
		print $fh $fasta;
	}
	$fh->close();
	exit(0);
}
sub GC{
	my $file = shift;
	my $in;
	my ($Gcontent,$Ccontent,$totlength) = (0,0,0);
	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }

	while(my $seqobj = $in->next_seq() ){
		$totlength += $seqobj->length;
		$Gcontent +=() = $seqobj->seq =~ m/G/g;
		$Ccontent +=() = $seqobj->seq =~ m/C/g;
	}
	say sprintf("%.2f",($Gcontent+$Ccontent)*100/$totlength),"%";

}

sub GC_Nfree{
	my $file = shift;
	my $in;
	my ($Gcontent,$Ccontent,$totlength,$totN) = (0,0,0,0);
	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }

	while(my $seqobj = $in->next_seq() ){
		$totlength += $seqobj->length;
		$totN +=() = $seqobj->seq =~ m/N/g ;
		$Gcontent +=() = $seqobj->seq =~ m/G/g;
		$Ccontent +=() = $seqobj->seq =~ m/C/g;
	}
	say sprintf("%.2f",($Gcontent+$Ccontent)*100/($totlength-$totN)),"%";
}

sub do_meta_Scf_report{

	my $file = shift;
	my ($tol_len,$N50,$max,$seq_count,$average_len,$in);
	$seq_count = 0;
	$max = 0;
	$cutofflenth = 0 unless defined $cutofflenth ;
	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }


	my @tot;

	while(my $seqobj = $in->next_seq() ){
		next unless $seqobj->length() >= $cutofflenth;
		push @tot,$seqobj->length() ;
		$seq_count++;

		$seqobj->length() >= $max ?  $max = $seqobj->length(): next ;
	}
	
	@tot = sort {$a <=> $b } map {$_ >= 0 ? $_:() } @tot;
	$tol_len = sum(@tot);
	$average_len = sprintf("%.2f", $tol_len / $seq_count);
	my $thr = 50*$tol_len/100;	
	my $pos = 0;
	for(@tot){
		$pos += $_;
		if($pos >= $thr){
			$N50 = $_;
			last;
		}
	}
	say join "\t",("Sample ID","Scaffold Number	Total Length(bp)","N50 Scaffold(bp)","Max Scaffold(bp)","Mean Length(bp)");
	return join "\t", ("",$seq_count,$tol_len,$N50,$max,$average_len);
  
}

sub do_meta_CDS_report{

	my $file = shift;
	my ($tol_len,$seq_count,$average_len,$in);
	$seq_count = 0;
	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	my @tot;
	while(my $seqobj = $in->next_seq() ){
		push @tot,$seqobj->length() ;
		$seq_count++;
	}
	$tol_len = sum(@tot);
	$average_len = sprintf("%.2f", $tol_len / $seq_count);
	return join "\t", ($seq_count,$tol_len,$average_len);
}

sub N50{
	my $file = shift;
	my $in;
	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	my @tot;
	while(my $seqobj = $in->next_seq() ){
	push @tot,$seqobj->length() ;
	}
	@tot = sort {$a <=> $b } map {$_ >= 0 ? $_:() } @tot;
	my $thr = 50*sum(@tot)/100;	
	my $pos = 0;
	for(@tot){
		$pos += $_;
		if($pos >= $thr){
			return $_;
			last;
		}
	}
}

sub N90{
	my $file = shift;
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	my @tot;
	while(my $seqobj = $in->next_seq() ){
	push @tot,$seqobj->length() ;
	}
	@tot = sort {$a <=> $b } map {$_ >= 0 ? $_:() } @tot;
	my $thr = 90*sum(@tot)/100;	
	my $pos = 0;
	for(@tot){
		$pos += $_;
		if($pos >= $thr){
			return $_;
			last;
		}
	}
}

sub Total_len {
	my $file = shift;
	my $in;
     if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }  
	my $tot;
	while(my $seqobj = $in->next_seq() ){
		$tot += $seqobj->length() ;
	}
	return $tot;
}

sub max_len {
	my $file = shift;
	my $in;
    if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
	my $max_len = 0;
	while(my $seqobj = $in->next_seq() ){
		$seqobj->length() >= $max_len ?  $max_len = $seqobj->length() : next;
	}
	return $max_len;
}



sub min_len {
	my $file = shift;
	my $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	my $max_len  = 100000;
	while(my $seqobj = $in->next_seq() ){

		$seqobj->length() <= $max_len ?  $max_len = $seqobj->length() : next;
	}
	return $max_len;
}

sub mean_len{
	my $file = shift;
	my $in; 
   	if($file){
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	 }else{
		 $in = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDIN);
	 }
    my $lab = 1;
	my $sum = 0;
	while(my $seqobj = $in->next_seq() ){
		$sum += $seqobj->length;
		$lab ++;
	}
	return sprintf ("%.2f",$sum / $lab);
}



sub filter{
	my $file = shift;
	my $outfile = shift;
	my $out;
	my $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	$cutofflenth = 0 unless defined $cutofflenth;
	
	if($outfile){
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	 }else{
		 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh => \*STDOUT);
	 }
	while(my $seqobj = $in->next_seq() ){
		next unless $seqobj->length >= $cutofflenth;
		$out->write_seq($seqobj)	
	}
	exit(0);
}


sub mutiple2uniq{

#### UDATE #########################
#
# to speed the program up, use fastx 
#
#####################################

	my $fastx = "/share/software/software/fastx_toolkit-0.0.14_install/bin/fastx_collapser";
	my $outfile = shift;
	my @files = @_;
	#my (%hash,$in,$out);
	#if (defined $outfile){
	#	 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
	#}else{
	#	 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh =>  \*STDOUT );
	#}	
	#foreach(@files){
	#	$in = Bio::SeqIO->new( -format => 'fasta' ,-file => $_);	
	#	while(my $seqobj = $in->next_seq() ){
	#		if (exists $hash{$seqobj->seq}){next}
	#		$out->write_seq($seqobj);
	#		$hash{$seqobj->seq} = 1;	
	#	}
	#}
	
	my $Files = join " ",@files;
	`cat $Files | $fastx > $outfile `;
	exit(0);
}


######## CUT Fasta file  ###########

=head1 
B<cut_fasta_file()>

usage:
	cut(<input file path>,<output file path>,<the seq # in the div file>)
	
	cut the fasta file into fitfull size 
	return 0 if work done 

=cut 

sub cut_fasta_file {
	my $file = shift;
	$file =	 ABSOLUTE_DIR($file);
	my $outpath = shift;

   if( defined $outpath && -d $outpath){

		$outpath = ABSOLUTE_DIR($outpath);
	}elsif(defined $outpath){
		MKDIR($outpath);
		$outpath = ABSOLUTE_DIR($outpath);
	}

#	die "Output path doesn't exists ,Please Check \n" unless -d $outpath ;
	my $one_file_size = shift;
	$one_file_size = 1000 unless $one_file_size;
	my $fasta = basename($file);
	my $flag = 1;
	my $in = Bio::SeqIO->new( -format => "fasta", -file => $file);
	my $init = ABSOLUTE_DIR($outpath)."/1"."$fasta";
	my $out = Bio::SeqIO->new(-format => "fasta" , -file => ">$init");
	while(my $seqobj = $in->next_seq() ){
		
		if ( $flag % ($one_file_size) != 0){
			$out->write_seq($seqobj);
			$flag ++;
			next;
		}else{
			my $num = $flag/($one_file_size);
 			$out = Bio::SeqIO->new(-format => "fasta" ,
			   		-file => ">".ABSOLUTE_DIR($outpath)."/".$num.$fasta) ;
			$out->write_seq($seqobj);
			$flag ++;
		}
	}	
	say "Cut Fasta Work Done!";
	exit 0;
}
sub meta_scaffold_title{
	say "Sample ID	Scaffold Number	Total Length(bp)	N50 Scaffold(bp)	Max Scaffold(bp)	Mean Length(bp)";	
	exit 0;
}

sub meta_CDS_title{
	say "Sample ID	CDSs Number	Total Length(bp)	Average Length(bp)";
	exit 0 ;
}

sub grep_by_name{
	my $file = shift;
	my $outfile = shift;
	my $query = shift;
	my $out ;
	my $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);

	if (defined $outfile){
			 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
		}else{
			 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh =>  \*STDOUT );
	}

	while(my $seqobj = $in->next_seq()){
			$out->write_seq($seqobj)	if	$seqobj->display_id =~ /$query/ 
	 }
	undef $in;
}





sub delete_gap{
	my $file = shift;
	my $outfile = shift;
	my $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);

	my %hash;
	my @array2;

	while(my $seqobj = $in->next_seq()){
		$hash{$seqobj->length} = 1;
		die "The Sequences is not in same length in this fasta File\n" if scalar(keys %hash ) == 2;
		my @array = split "",$seqobj->seq();
		my $lab = 0;
		foreach(@array){
			push @array2,$lab unless $_ ne "." and $_ ne "-";
			$lab++;
		}
		@array2 = uniq(@array2);
	 }
	undef $in;
	$in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	my $out;
	if (defined $outfile){
			 $out = Bio::SeqIO->new( -format => 'fasta' ,-file => ">$outfile");
		}else{
			 $out = Bio::SeqIO->new( -format => 'fasta' ,-fh =>  \*STDOUT );
	}
	while(my $seqobj = $in->next_seq()){
		my @array = split "",$seqobj->seq();
		my $name = $seqobj->display_id();
		my $seq_ungap = "";
		my $i = -1;
		foreach (@array){
			$i ++ &&  next if ($i ~~ @array2);
			$seq_ungap .= $_;
			$i++;
		}
		my $ungap =Bio::Seq->new(
			-display_id => $name,
			-seq => $seq_ungap
		);
		$out->write_seq($ungap);
	 }
	 exit 0;
}

sub STATIC_AMPHI_N{
	my $file = shift;
	my $in = Bio::SeqIO->new( -format => 'fasta' ,-file => $file);
	while(my $seqobj = $in->next_seq() ){
		my $tail = substr($seqobj->seq(),-101);
		my $head = substr($seqobj->seq(),0,100);
		my @t = $tail =~ /N/ig;
		my @h = $head =~ /N/ig;
		say join "\t",($seqobj->display_id,scalar(@t),scalar(@h));
	}
	exit(0);
}





