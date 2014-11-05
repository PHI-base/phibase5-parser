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
my $protein_name_count = 0;
my $no_protein_name_count = 0;

# open output files
my $protein_filename = './output/uniprot_protein_names.tsv';
my $no_protein_filename = './output/uniprot_no_protein_names.tsv';
open (PROTEIN_NAME_FILE, "> $protein_filename") or die "Error opening output file\n";
open (NO_PROTEIN_NAME_FILE, "> $no_protein_filename") or die "Error opening output file\n";

print "Printing protein names to output file $protein_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     # RESTful URL query to get protein names for the current UniProt accession
     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=protein%20names";

     # execute query and process response
     my $protein_names_response = query_uniprot($query);
     my @protein_names_plus_header = split ("\n",$protein_names_response); # split into header & protein names
     my $protein_names_string = $protein_names_plus_header[1]; # the protein names string is second element, after the header

     if (defined $protein_names_string) {
        print PROTEIN_NAME_FILE "$uniprot_acc\t";
        my @protein_names = split ("\t",$protein_names_string); # split into array of individual protein names (not sure what the delimiter should be - so far only one name given for each UniProt entry)
        foreach my $protein_name (@protein_names) 
        {
          $protein_name_count++;
          print PROTEIN_NAME_FILE "$protein_name\t";
        }
        print PROTEIN_NAME_FILE "\n";
     } else {
       # if not defined, there were no EMBL IDs for the UniProt entry
       $no_protein_name_count++;
       print "No protein names found for UniProt accession $uniprot_acc\n";
       print NO_PROTEIN_NAME_FILE "$uniprot_acc\n";
     }

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (PROTEIN_NAME_FILE);
close (NO_PROTEIN_NAME_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "Protein names retrieved:$protein_name_count, available in output file $protein_filename\n";
print "UniProt accessions without an EMBL ID:$no_protein_name_count, availalbe in file $no_protein_filename\n";

