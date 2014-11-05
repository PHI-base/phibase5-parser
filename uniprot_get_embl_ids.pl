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

# iterators to gather statistics
my $uniprot_count = 0;
my $embl_id_count = 0;
my $no_embl_id_count = 0;

# open output files
my $embl_filename = './output/uniprot_embl_ids.tsv';
my $no_embl_filename = './output/uniprot_no_embl_ids.tsv';
open (EMBL_FILE, "> $embl_filename") or die "Error opening output file\n";
open (NO_EMBL_FILE, "> $no_embl_filename") or die "Error opening output file\n";

print "Printing EMBL identifiers to output file $embl_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     # RESTful URL query to get EMBL IDs for the current UniProt accession
     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=database(EMBL)";

     # execute query and process response
     my $embl_ids_response = query_uniprot($query);
     my @embl_ids_plus_header = split ("\n",$embl_ids_response); # split into header & EMBL IDs
     my $embl_ids_string = $embl_ids_plus_header[1]; # the EMBL IDs string is second element, after the header
     if (defined $embl_ids_string) {
       print EMBL_FILE "$uniprot_acc\t";
       my @embl_ids = split (";",$embl_ids_string); # split into array of individual EMBL IDs, wrich are delimited by semi-colon
       # print each EMBL ID to tab-delimited file (for most UniProt entries there will only be one EMBL ID)
       foreach my $embl_id (@embl_ids) 
       {
         $embl_id_count++;
         print EMBL_FILE "$embl_id\t";
       }
       print EMBL_FILE "\n";
     } else {
       # if not defined, there were no EMBL IDs for the UniProt entry
       $no_embl_id_count++;
       print "No EMBL ID found for UniProt accession $uniprot_acc\n";
       print NO_EMBL_FILE "$uniprot_acc\n";
     }

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (EMBL_FILE);
close (NO_EMBL_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "EMBL IDs retrieved:$embl_id_count, available in output file $embl_filename\n";
print "UniProt accessions without an EMBL ID:$no_embl_id_count, availalbe in file $no_embl_filename\n";

