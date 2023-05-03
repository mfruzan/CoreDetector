#!/usr/bin/perl -w
use strict;

# open map 
my $file = $ARGV[0];
open(IN,'<',$file) || die "Could not open $file: $!\n";

my %map=();
while(my $line = <IN>){
	chomp $line;
	my @arr = split(/	/, $line);
	my $org = $arr[0];
	my $new = $arr[1];
	$map{$org} = $new;
}

close(IN);
	
my $file2 = $ARGV[1];
open(IN2,'<',$file2) || die "Could not open $file: $!\n";
foreach my $lin ( <IN2>){

	while (my( $key,$value) = each %map){
		$lin =~ s/${key}/$value/;
		#print " $key and $value\n";
	}
	print "$lin";
}
close(IN2);
