#!/usr/bin/perl

#read file in from input line
$infile = $ARGV[0];
open(TXT, "<$infile");

#read in the DNA string using the fasta subfunction
$DNA = &read_fasta();

$len = length($DNA);
print "\n DNA Length is: $len \n";

$numG=0;
$numC=0;
$numT=0;
$numA=0;

@bases=split(//,$DNA);

foreach $bp(@bases)
{
    if($bp =~ m/G/i){$numG++};
    if($bp =~ m/C/i){$numC++};
    if($bp =~ m/T/i){$numT++};
    if($bp =~ m/A/i){$numA++};
}

print "\n Number of G bases: $numG";
print "\n Number of C bases: $numC";
print "\n Number of T bases: $numT";
print "\n Number of A bases: $numA";

$GC_content = (($numG+$numC)/$len)*100;

print "\n\n GC Content is: $GC_content % \n";

close(TXT);

sub read_fasta
{
    $sequence = "";
    while(<TXT>)
    {   
        $line = $_;
        print " $line \n";

       #remove newline characters 
       chomp($line);
        # discard fasta header line
        if($line =~ />/){ next }
       # append the line to the DNA sequence
        else{ $sequence .= $line }
    } 
    return($sequence);
}

