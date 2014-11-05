use strict;
use warnings;
use LWP::UserAgent;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all UniProt accessions from PHI-base
my $sql_query = qq(SELECT DISTINCT uniprot_accession FROM pathogen_gene_mutant;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# iterator to count rows
my $uniprot_count = 0;

# open output file
my $TEST_filename = './output/uniprot_TEST.tsv';
open (TEST_FILE, "> $TEST_filename") or die "Error opening output file\n";

print "Printing to TEST output file $TEST_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=database(UniParc)";

     # execute query and process response
     my $TEST_response = query_uniprot($query);
     my @TEST_plus_header = split ("\n",$TEST_response); # split into header & EMBL IDs
     my $TEST_string = $TEST_plus_header[1]; # the EMBL IDs string is second element, after the header
#     my @TEST = split (";",$TEST_string); # split into array of individual EMBL IDs, which are delimited by semi-colon

#     foreach my $TEST (@TEST) 
#     {
#       print TEST_FILE "$TEST_id\t";
        print TEST_FILE "$uniprot_acc\t";
        #print TEST_FILE "$TEST_string\n";
        print TEST_FILE "$TEST_response\n";
#     }
#     print TEST_FILE "\n";

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (TEST_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "Output file: $TEST_filename\n";

