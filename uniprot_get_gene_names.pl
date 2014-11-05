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
my $gene_name_count = 0;
my $no_gene_name_count = 0;

# open output files
my $gene_filename = './output/uniprot_gene_names.tsv';
my $no_gene_filename = './output/uniprot_no_gene_names.tsv';
open (GENE_NAME_FILE, "> $gene_filename") or die "Error opening output file\n";
open (NO_GENE_NAME_FILE, "> $no_gene_filename") or die "Error opening output file\n";

print "Printing gene names to output file $gene_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     # RESTful URL query to get gene names for the current UniProt accession
     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=genes";

     # execute query and process response
     my $gene_names_response = query_uniprot($query);
     my @gene_names_plus_header = split ("\n",$gene_names_response); # split into header & gene names
     my $gene_names_string = $gene_names_plus_header[1]; # the gene names string is second element, after the header
     if (defined $gene_names_string) {
        print GENE_NAME_FILE "$uniprot_acc\t";
        my @gene_names = split (" ",$gene_names_string); # split into array of individual gene names
        foreach my $gene_name (@gene_names) 
        {
           $gene_name_count++;
           print GENE_NAME_FILE "$gene_name\t";
        }
        print GENE_NAME_FILE "\n";
     } else {
       # if not defined, there were no gene names for the UniProt entry
       $no_gene_name_count++;
       print "No gene names found for UniProt accession $uniprot_acc\n";
       print NO_GENE_NAME_FILE "$uniprot_acc\n";
     }
     
  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (GENE_NAME_FILE);
close (NO_GENE_NAME_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "Gene names retrieved:$gene_name_count, available in output file $gene_filename\n";
print "UniProt accessions without a gene name:$no_gene_name_count, availalbe in file $no_gene_filename\n";

