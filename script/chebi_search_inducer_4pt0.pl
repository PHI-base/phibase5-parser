#!/usr/bin/perl
use strict;
use warnings;
use LWP::UserAgent;
use SOAP::Lite;

# parse text file that contains the inducers from the spreadsheet
# saving the value in an array
open (INDUCER_FILE, "../input/inducer_chemicals_in_phibase_4pt0.tsv") or die "Error opening input file\n";

# iterators to gather statistics
my $inducer_count = 0;
my $chebi_found_count = 0;
my $possible_chebi_count = 0;
my $chebi_not_found_count = 0;

# open output file
my $all_inducers_filename = '../mapping/inducers_4pt0_with_chebi.tsv';
my $chebi_found_filename = '../mapping/inducers_4pt0_chebi_found.tsv';
my $possible_chebi_filename = '../mapping/inducers_4pt0_possible_chebi.tsv';
my $no_chebi_filename = '../error/inducers_4pt0_no_chebi.tsv';
open (INDUCER_OUTPUT_FILE, "> $all_inducers_filename") or die "Error opening output file\n";
open (CHEBI_FILE, "> $chebi_found_filename") or die "Error opening output file\n";
open (POSSIBLE_CHEBI_FILE, "> $possible_chebi_filename") or die "Error opening output file\n";
open (NO_CHEBI_FILE, "> $no_chebi_filename") or die "Error opening output file\n";

print "Printing to inducer output file $all_inducers_filename...\n";

# Setup SOAP web service, using the ChEBI WSDL XML as description of available fields
my $WSDL = 'http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl';
my $nameSpace = 'http://www.ebi.ac.uk/webservices/chebi';
my $soap = SOAP::Lite
-> uri($nameSpace)
-> proxy($WSDL);

# Setup method to get ChEBI entity data
my $method = SOAP::Data->name('getLiteEntity')
             ->attr({xmlns => $nameSpace});

# each row of the file contains an inducer from the spreadsheet
# for each, check if the chemical name can be found in ChEBI
while (<INDUCER_FILE>) {

   chomp;
   my $inducer = $_;

   # increment number of inducers
   $inducer_count++;

   # Set parameters, searching using inducer name, then call the SOAP method
   my @params = ( SOAP::Data->name(search => $inducer), 
                  SOAP::Data->name(maximumResults => 1),
                  SOAP::Data->name(searchCategory => 'ALL NAMES') 
                );
   my $som = $soap->call($method => @params);

   # Retrieve ChEBI ID and name,
   # based on metadata fields available from the ChEBI WSDL XML definitions
   my @ids = $som->valueof('//chebiId');
   my @names = $som->valueof('//chebiAsciiName');

   # check if ChEBI ID is available
   if (@ids) {

      my $chebi_id = shift @ids;
      my $chebi_name = shift @names;

      # only regarded as the same chemical if
      # one name is a substring of the other
      if ($chebi_name=~/\Q$inducer\E/i or $inducer=~/\Q$chebi_name\E/i) {

        $chebi_found_count++;

        # print the ChEBI ID, names & synonyms
        print CHEBI_FILE "$inducer\t$chebi_id\t$chebi_name\n";
        # also add these details to all inducers file
        print INDUCER_OUTPUT_FILE "$inducer\t$chebi_id\t$chebi_name\n";

      } else {

        # The chemicals here were returned as the top result by ChEBI,
        # but the ChEBI and Inducer names do not match,
        # so they need to be looked through manually

        # print the ChEBI ID, names & synonyms
        $possible_chebi_count++;
        print POSSIBLE_CHEBI_FILE "$inducer\t$chebi_id\t$chebi_name\n";

        # since ChEBI ID needs confirmation,
        # only add inducer name to all inducer file
        print INDUCER_OUTPUT_FILE "$inducer\n";

      }

   } else { # ChEBI entry not found

      $chebi_not_found_count++;
      print NO_CHEBI_FILE "$inducer\n";
      print INDUCER_OUTPUT_FILE "$inducer\n";

   } 

   # print message for every 20th inducer processed
   print "Inducers processed:$inducer_count\n" unless ($inducer_count % 20);

} # end while rows

close (INDUCER_FILE);
close (INDUCER_OUTPUT_FILE);
close (CHEBI_FILE);
close (POSSIBLE_CHEBI_FILE);
close (NO_CHEBI_FILE);

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique inducers from PHI-base 3.6 spreadsheet:$inducer_count, available with ChEBI IDs in output file $all_inducers_filename\n";
print "ChEBI entries retrieved:$chebi_found_count, available in output file $chebi_found_filename\n";
print "Possible ChEBI matches - need manual checking:$possible_chebi_count, available in output file $possible_chebi_filename\n";
print "ChEBI IDs not found:$chebi_not_found_count, availalbe in file $no_chebi_filename\n";

