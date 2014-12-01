#!/usr/bin/perl
use strict;
use warnings;

#-----------------------------------------------#
#  Retrieve data from the Plant Ontology      #
#  using their Browser Rest API                 #
#                                               #
#  Results are returned in JSON format          #
#-----------------------------------------------#

use LWP::Simple;
use JSON;
use Data::Dumper;

# counters to gather statistics
my $po_term_count = 0;
my $po_term_found_count = 0;
my $po_term_not_found_count = 0;

# open output files
my $po_term_filename = './output/plant_ontology_terms.tsv';
my $no_po_term_filename = './output/plant_ontology_terms_not_found.tsv';
open (PO_TERM_FILE,"> $po_term_filename") or die "Error opening output file\n";
open (NO_PO_TERM_FILE,"> $no_po_term_filename") or die "Error opening output file\n";
print "Printing Plant Ontology terms to output file $po_term_filename...\n";

# plant ontology API url
my $url = "http://palea.cgrb.oregonstate.edu/services/PO_web_service.php?request_type=term_detail&accession_id=";

# array of identifiers of plant ontology terms (POIDs)
my @list_of_poids =
  ("PO:0000252"
  );

# iterate through each of the IDs in the list
foreach my $poid (@list_of_poids) {

  $po_term_count++;

  # append the PO ID onto the URL to make the REST query
  my $rest_query = $url.$poid;

  # use LWP simple to retrieve JSON output from REST query
  my $json_result = get($rest_query);
  #my $query_result;

  # check if the plant ontology term has been found
  if (defined $json_result) {

    $po_term_found_count++;

    # decode the JSON format to text
    my $query_result = decode_json $json_result;

    # get the name and definition of the plant ontology term
    my $name = $query_result->{'PO_term_detail_response'}[0]{'name'};
    my $definition = $query_result->{'PO_term_detail_response'}[0]{'definition'};

    # output the plant ontology details to file
    print PO_TERM_FILE "$poid\nPlant Ontology term:$name\n";
    print PO_TERM_FILE "Definition:$definition\n\n";

  } else { # Plant Ontology term not found

    $po_term_not_found_count++;
    print STDERR "ERROR: Plant Ontology ID $poid not found\n";
    print NO_PO_TERM_FILE "$poid\n";

  } # end if plant ontology ID found

  # print message for every 20th plant ontology term processed
  print "Plant Ontology terms processed:$po_term_count\n" unless ($po_term_count % 20);

} # end foreach Plant Ontology identifier

close (PO_TERM_FILE);
close (NO_PO_TERM_FILE);

print "Total plant ontology IDs: $po_term_count\n";
print "Plant Ontology IDs found: $po_term_found_count, output file $po_term_filename\n";
print "Plant Ontology IDs not found: $po_term_not_found_count, output file $no_po_term_filename\n";

