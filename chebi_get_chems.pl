#!/usr/bin/perl
use strict;
use warnings;
use LWP::UserAgent;
use DBI; # load perl postgresql module
use SOAP::Lite;

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all (non-empty) chebi identifiers from PHI-base
my $sql_query = qq(SELECT chebi_id FROM chemical WHERE chebi_id != '';);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# iterators to gather statistics
my $chebi_count = 0;
my $chebi_found_count = 0;
my $chebi_not_found_count = 0;

# open output file
my $chebi_filename = './output/chebi_chem_details.txt';
my $no_chebi_filename = './output/chebi_no_chem_details.txt';
open (CHEBI_FILE, "> $chebi_filename") or die "Error opening output file\n";
open (NO_CHEBI_FILE, "> $no_chebi_filename") or die "Error opening output file\n";

print "Printing to ChEBI output file $chebi_filename...\n";

# Setup SOAP web service, using the ChEBI WSDL XML as description of available fields
my $WSDL = 'http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl';
my $nameSpace = 'http://www.ebi.ac.uk/webservices/chebi';
my $soap = SOAP::Lite
-> uri($nameSpace)
-> proxy($WSDL);

# Setup method to get ChEBI entity data
my $method = SOAP::Data->name('getCompleteEntity')
             ->attr({xmlns => $nameSpace});

# iterate through each ChEBI ID
while (my @row = $sql_stmt->fetchrow_array()) {

   my $chebi_id = shift @row;

   # increment number of ChEBI ids
   $chebi_count++;

   # Set ChEBI ID as parameter, then call the SOAP method
   my @params = ( SOAP::Data->name(chebiId => $chebi_id));
   my $som = $soap->call($method => @params);

   # Retrieve chemical data, based on metadata fields available from the ChEBI WSDL XML definitions
   my @names = $som->valueof('//chebiAsciiName');
   my @synonyms = $som->valueof('//Synonyms//data');
   my @registry = $som->valueof('//RegistryNumbers//data');
   my @registry_types = $som->valueof('//RegistryNumbers//type');

   # check if names are available
   # if they are, print the values to file, 
   # if no names are available, assume ChEBI entry was not found
   if (@names) {

      $chebi_found_count++;

      # print the ChEBI ID, names & synonyms
      print CHEBI_FILE "\nCHEBI:$chebi_id\n";
      print CHEBI_FILE "Name: $_\n" foreach @names; 
      print CHEBI_FILE "Synonym: $_\n" foreach @synonyms;

      # only interested in the CAS registry, not other registry numbers
      # so test the registry type, then print the corresponding CAS registry number
      # (which will have an identical index to the registry type)
      my $registry_index = 0;
      foreach my $registry_type (@registry_types) {
        if ($registry_type eq 'CAS Registry Number') {
          print CHEBI_FILE "CAS Registry: $registry[$registry_index]\n";
          # since CAS reg has been found, can not exit loop
          last;
        }
        $registry_index++;
      }

   } else { # ChEBI entry not found
      $chebi_not_found_count++;
      print STDERR "CHEBI:$chebi_id not found\n";
      print NO_CHEBI_FILE "CHEBI:$chebi_id\n";
   } 

   # print message for every 20th ChEBI entry processed
   print "ChEBI entries processed:$chebi_count\n" unless ($chebi_count % 20);

} # end while rows

close (CHEBI_FILE);
close (NO_CHEBI_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Total ChEBI IDs from chemical table:$chebi_count\n";
print "ChEBI entries retrieved:$chebi_found_count, available in output file $chebi_filename\n";
print "ChEBI IDs not found:$chebi_not_found_count, availalbe in file $no_chebi_filename\n";


