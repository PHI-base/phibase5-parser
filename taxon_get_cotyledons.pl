use strict;
use warnings;

use LWP::Simple;
use XML::Twig;
use DBI;

use phibase_subroutines qw(connect_to_phibase); # load PHI-base functions
my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all host taxon IDs from PHI-base
my $sql_query = qq(SELECT DISTINCT ncbi_taxon_id FROM interaction_host;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $host_taxon_count = 0;
my $monocot_count = 0;
my $dicot_count = 0;
my $non_plant_count = 0;

# open output file
my $cotyledon_filename = './output/cotyledon_details.tsv';
open (COTYLEDON_FILE,"> $cotyledon_filename") or die "Error opening output file\n";
print "Printing cotyledon details to output file $cotyledon_filename...\n";

# iterate through each taxon ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $host_taxon_count++;
  foreach my $host_taxon_id (@row) {

     $host_taxon_id =~ s/^\s+//; # remove blank space from start of taxon IDs
     $host_taxon_id =~ s/\s+$//; # remove blank space from end of taxon IDs

     my $query = "http://www.ebi.ac.uk/ena/data/view/Taxon:$host_taxon_id&display=xml";
     my $xml_response = get $query;
     die "Error getting $query" unless defined $xml_response;

     # use XML twig to parse the XML data
     my $xml_twig = XML::Twig->new();
     $xml_twig->parse($xml_response);

     # parse the XML data to get the relevant host taxon info
     my $host_taxon = $xml_twig->root->first_child('taxon');
     my $host_taxon_id = $host_taxon->{'att'}->{'taxId'};
     my $host_taxon_name = $host_taxon->{'att'}->{'scientificName'};

     # need to check if common name exists for this taxon
     my $host_taxon_common_name;
     if ($host_taxon->{'att'}->{'commonName'}) {
       $host_taxon_common_name = "$host_taxon->{'att'}->{'commonName'}"
     } else {
       $host_taxon_common_name = "none";
     }

     # print the taxon info to file
     print COTYLEDON_FILE "$host_taxon_id\t$host_taxon_name\t$host_taxon_common_name\t";

     # get all the taxon ids for the lineage of the current taxon
     my @lineage_taxons = $host_taxon->first_child('lineage')->children('taxon');

     # flag indicating if mono or dicot found
     my $mono_dicot_found = 0;

     # check each of the lineage taxon IDs against the taxon IDs for monocot (4447) or dicot (71240)
     # then print the appropriate label for the host taxon
     foreach my $lineage_taxon (@lineage_taxons) {

        my $lineage_taxon_id = $lineage_taxon->{'att'}->{'taxId'};

        if ($lineage_taxon_id eq "4447") {
           print COTYLEDON_FILE "Monocot\n";
           $mono_dicot_found = 1;
           $monocot_count++;
           last;
        }
        elsif ($lineage_taxon_id eq "71240") {
           print COTYLEDON_FILE "Dicot\n";
           $mono_dicot_found = 1;
           $dicot_count++;
           last;
        }

     } # end foreach lineage taxon

     # if neither monocot or dicot found in lineage,
     # then assume the host is not a plant species
     unless ($mono_dicot_found) {
        print COTYLEDON_FILE "n\\a\n";
        $non_plant_count++;
     }

  } # end foreach host taxon

} # end while rows

close (COTYLEDON_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Unique host taxon IDs: $host_taxon_count\n";
print "Monocot hosts: $monocot_count\n";
print "Dicot hosts: $dicot_count\n";
print "Non-plant hosts: $non_plant_count\n";

