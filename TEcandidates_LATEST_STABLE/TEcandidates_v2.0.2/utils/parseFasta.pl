#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Seq;

my $fastaFile = shift or die "$!\n";

my $fastaSeqIO = new Bio::SeqIO(-file=>"$fastaFile",-format=>"fasta");

while (my $seqObject = $fastaSeqIO->next_seq){
	print ">".$seqObject->id." ".$seqObject->desc."\n";
	print $seqObject->seq."\n";
}
