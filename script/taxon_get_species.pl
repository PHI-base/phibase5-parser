use strict;
use warnings;

use LWP::Simple;
use XML::Twig;
use DBI;

use phibase_subroutines qw(connect_to_phibase); # load PHI-base functions
my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all pathogen taxon IDs from PHI-base
my $sql_query = qq(SELECT DISTINCT ncbi_taxon_id FROM pathogen_gene);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $path_taxon_count = 0;
my $path_species_count = 0;
my $path_no_species_count = 0;

# array for list of unique species
my @unique_species;

# open output file
my $taxon_species_filename = '../output/path_taxon_species.tsv';
open (TAXON_SPECIES_FILE,"> $taxon_species_filename") or die "Error opening output file\n";
print "Printing pathogen taxon species to output file $taxon_species_filename...\n";

# iterate through each taxon ID
while (my @row = $sql_stmt->fetchrow_array()) {

  $path_taxon_count++;

  foreach my $path_taxon_id (@row) {

     $path_taxon_id =~ s/^\s+//; # remove blank space from start of taxon IDs
     $path_taxon_id =~ s/\s+$//; # remove blank space from end of taxon IDs

     my $query = "http://www.ebi.ac.uk/ena/data/view/Taxon:$path_taxon_id&display=xml";
     my $xml_response = get $query;
     die "Error getting $query" unless defined $xml_response;

     # use XML twig to parse the XML data
     my $xml_twig = XML::Twig->new();
     $xml_twig->parse($xml_response);

     # parse the XML data to get the relevant pathogen taxon info
     my $path_taxon = $xml_twig->root->first_child('taxon');
     my $path_taxon_rank = $path_taxon->{'att'}->{'rank'};
     my $path_taxon_name = $path_taxon->{'att'}->{'scientificName'};

     # if the pathogen taxon was not found,
     # display an error message and move to next pathogen 
     if (not defined $path_taxon_name) {
       print STDERR "ERROR: Pathogen taxon ID $path_taxon_id not found\n";
       $path_no_species_count++;
       next;
     }

     # print the taxon ID and name to the output file
     print TAXON_SPECIES_FILE "$path_taxon_id\t$path_taxon_name\t";

     # variables for the species taxon ID and name
     my $species_taxon_id;
     my $species_taxon_name;

     # if the given pathogen taxon happens to be a species,
     # then simply assign the taxon ID and name to the species,
     # otherwise, search through the lineage for the current pathogen taxon
     # until the species rank is found
     if (defined $path_taxon_rank and $path_taxon_rank eq "species") {

       $species_taxon_id = $path_taxon_id;
       $species_taxon_name = $path_taxon_name;

     } else {

       # get all the taxon ids for the lineage of the current taxon
       my @lineage_taxons = $path_taxon->first_child('lineage')->children('taxon');

       # check each of the lineage taxon ranks to find the species rank
       # then print the corresponding taxon ID and name
       foreach my $lineage_taxon (@lineage_taxons) {

          my $lineage_taxon_rank = $lineage_taxon->{'att'}->{'rank'};
          my $lineage_taxon_id = $lineage_taxon->{'att'}->{'taxId'};
          my $lineage_taxon_name = $lineage_taxon->{'att'}->{'scientificName'};

          if (defined $lineage_taxon_rank and $lineage_taxon_rank eq "species") {
             $species_taxon_id = $lineage_taxon_id;
             $species_taxon_name = $lineage_taxon_name;
             last;
          }

       } # end foreach lineage taxon

     } # end else not species

     # unless the pathogen species ID is already in the species list
     # add it to the list and increment the counter
     unless ($species_taxon_id ~~ @unique_species) {
        push (@unique_species,$species_taxon_id); 
        $path_species_count++;
     }

     # print the species taxon ID and ID to file 
     print TAXON_SPECIES_FILE "$species_taxon_id\t$species_taxon_name\n";

  } # end foreach pathogen taxon

} # end while rows

close (TAXON_SPECIES_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Unique pathogen taxon IDs: $path_taxon_count\n";
print "Number of pathogen species: $path_species_count\n";
print "Pathogen taxa not found: $path_no_species_count\n";

