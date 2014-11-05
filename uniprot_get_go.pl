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
my $go_term_count = 0;
my $no_go_term_count = 0;

# open output files
my $go_filename = './output/uniprot_go_terms.tsv';
my $no_go_filename = './output/uniprot_no_go_terms.tsv';
open (GO_FILE, "> $go_filename") or die "Error opening output file\n";
open (NO_GO_FILE, "> $no_go_filename") or die "Error opening output file\n";

print "Printing Gene Ontology terms to output file $go_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs
     #print GO_FILE "$uniprot_acc\t";

     # RESTful URL query to get Gene Ontology terms for the current UniProt accession (both the GO ID and the term itself)
     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=go-id,go";

     # execute query and process response
     my $go_terms_response = query_uniprot($query);
     my @go_terms_plus_header = split ("\n",$go_terms_response); # split into header & go terms
     my $go_terms_string = $go_terms_plus_header[1]; # the go terms string is second element, after the header

     if (defined $go_terms_string) {
 
       $go_terms_string =~ s/\s+$//; # remove blank space from end of string
       $go_terms_string =~ s/^\s+//; # remove blank space from start of string

       if ($go_terms_string ne "") {
         print GO_FILE "$uniprot_acc\n";

         # split into array of go ids (first element) and go terms (second element), delimited by tab
         my @go_term_set = split ("\t",$go_terms_string);

         # first element will be all of the GO IDs, each of which are separated by semi-colon
         my $go_id_set = shift @go_term_set;
         my @go_ids = split ("; ",$go_id_set);

         # second element will be all of the GO term names, each of which are separated by semi-colon
         my $go_term_name_set = shift @go_term_set;
         my @go_terms = split ("; ",$go_term_name_set);

         # loop through the list of GO IDs, then print the GO ID & its corresponding GO term name
         while (my $go_id = shift @go_ids) {
           $go_term_count++;
           my $go_term = shift @go_terms;
           print GO_FILE "$go_id\t$go_term\n";
         }
         print GO_FILE "\n";

       } else {
         # no GO terms found for the UniProt entry
         $no_go_term_count++;
         print "No GO term found for UniProt accession $uniprot_acc\n";
         print NO_GO_FILE "$uniprot_acc\n";
       }

     } else {
       # no data found for UniProt entry
       $no_go_term_count++;
       print "No data found for UniProt accession $uniprot_acc\n";
       print NO_GO_FILE "$uniprot_acc\n";
     }

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (GO_FILE);
close (NO_GO_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "Gene Ontology terms retrieved:$go_term_count, available in output file $go_filename\n";
print "UniProt accessions without a Gene Ontology term:$no_go_term_count, availalbe in file $no_go_filename\n";

