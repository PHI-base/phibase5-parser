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

# counters to gather statistics
my $uniprot_count = 0;
my $taxon_id_count = 0;
my $taxon_sci_name_count = 0;
my $taxon_strain_count = 0;
my $taxon_synonym_count = 0;
my $no_taxon_count = 0;

# open output files
my $taxon_filename = '../output/uniprot_taxon_details.tsv';
my $no_taxon_filename = '../error/uniprot_no_taxon_details.tsv';
open (TAXON_FILE, "> $taxon_filename") or die "Error opening output file\n";
open (NO_TAXON_FILE, "> $no_taxon_filename") or die "Error opening output file\n";

print "Printing taxon details to output file $taxon_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     # RESTful URL query to get taxon details for the current UniProt accession
     my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=organism-id,organism"; # can also use "taxon", in place of organism

     # execute query and process response
     my $taxon_details_response = query_uniprot($query);
     my @taxon_details_plus_header = split ("\n",$taxon_details_response); # split into header & taxon details
     my $taxon_details_string = $taxon_details_plus_header[1]; # the taxon details string is second element, after the header

     if (defined $taxon_details_string) {

        # split into array of individual taxon id & names,
        # delimited by either tab or open brace
        my @taxon_items = split (/[\t\(]/,$taxon_details_string);

        # the first element of the array will always be the taxon id
        my $taxon_id = shift @taxon_items;
        $taxon_id_count++;
        print TAXON_FILE "$uniprot_acc\nNCBI Taxon ID: $taxon_id\n";

        # the second element of the array will always be the scientific name of the species (genus followed by species name)
        my $taxon_sci_name = shift @taxon_items;
        $taxon_sci_name_count++;
        print TAXON_FILE "Scientific Name: $taxon_sci_name\n";
        
        # the remaining items will be a variety of strain and synonyms
        foreach my $taxon_name (@taxon_items) 
        {
          # remove close bracket ")" and trailing spaces from taxon synonyms
          $taxon_name = substr( $taxon_name, 0, index($taxon_name, ')') );

          # if the name begins with "strain" then it contains strain details
          # otherwise, it will be a species synonym
          if ($taxon_name =~ /^strain/) {
            $taxon_strain_count++;
            print TAXON_FILE "Strain: $taxon_name\n";
          } else {
            # assume this is a synonym of the species name
            $taxon_synonym_count++;
            print TAXON_FILE "Synonym: $taxon_name\n";
          }
        }
        print TAXON_FILE "\n";

     } else {
       # if not defined, there were no taxon details for the UniProt entry
       $no_taxon_count++;
       print STDERR "No taxon details found for UniProt accession $uniprot_acc\n";
       print NO_TAXON_FILE "$uniprot_acc\n";
     }

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (TAXON_FILE);
close (NO_TAXON_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "Taxon IDs retrieved:$taxon_id_count, available in output file $taxon_filename\n";
print "Taxon scientific names retrieved:$taxon_sci_name_count\n";
print "Taxon strain names retrieved:$taxon_strain_count\n";
print "Taxon synonym names retrieved:$taxon_synonym_count\n";
print "UniProt accessions without taxon details:$no_taxon_count, availalbe in file $no_taxon_filename\n";

