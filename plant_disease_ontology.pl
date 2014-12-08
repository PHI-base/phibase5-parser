#!/usr/bin/perl
use strict;
use warnings;

#-----------------------------------------------#
#  Retrieve data from the Disease Ontology      #
#  using their Browser Rest API                 #
#                                               #
#  Results are returned in JSON format          #
#-----------------------------------------------#

use LWP::Simple;
use JSON;
use Data::Dumper;


# parse tab-separated file that maps phenotype outcome values of the spreadsheet
# to the identifier in the phenotype outcome ontology
# saving the value and identifier as key/value pairs in a hash
open (PLANT_DISEASE_FILE, "ontologies/Disease/PlantDisease/plant_disease_ontology.obo") or die "Error opening input file\n";

# hash to map phenotype outcome text to ontology identifier
my %phenotype_outcome_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "phenotype outcome ontology identifier", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes the ontology identifier
while (<PLANT_DISEASE_FILE>) {
  chomp;
  my ($phen_outcome_value,$phen_outcome_ontology_id_list) = split(/\t/,$_);
  $phenotype_outcome_mapping{$phen_outcome_value} = $phen_outcome_ontology_id_list;
}
close (PLANT_DISEASE_FILE);



# counters to gather statistics
my $disease_count = 0;
my $disease_found_count = 0;
my $disease_not_found_count = 0;

# open output files
my $disease_filename = './output/disease_details.tsv';
my $no_disease_filename = './output/disease_not_found.tsv';
open (DISEASE_FILE,"> $disease_filename") or die "Error opening output file\n";
open (NO_DISEASE_FILE,"> $no_disease_filename") or die "Error opening output file\n";
print "Printing Disease Ontology terms to output file $disease_filename...\n";

# disease ontology API url
my $url = "http://www.disease-ontology.org/api/metadata/";

# array of identifiers of disease ontology terms (DOIDs)
# the examples here are human and animal diseases found in PHI-base
my @list_of_doids =
  ("DOID:1498","DOID:8913","DOID:9588","DOID:1508","DOID:9965",
   "DOID:13564","DOID:0050153","DOID:10554","DOID:9065",
   "DOID:12140","DOID:1731","DOID:2326","DOID:13250","DOID:0060185",
   "DOID:1770","DOID:399","DOID:11729","DOID:9065","DOID:6132",
   "DOID:3488","DOID:14262","DOID:552","DOID:11573","DOID:8778",
   "DOID:10754","DOID:9471","DOID:1579","DOID:37","DOID:2275",
   "DOID:8712","DOID:14115","DOID:10314","DOID:848","DOID:12053",
   "DOID:10398","DOID:10554","DOID:2123","DOID:12385"
  );

# iterate through each of the IDs in the list
foreach my $doid (@list_of_doids) {

  $disease_count++;

  # append the DO ID onto the URL to make the REST query
  my $rest_query = $url.$doid;

  # use LWP simple to retrieve JSON output from REST query
  my $json_result = get($rest_query);
  #my $query_result;

  # check if the disease has been found
  if (defined $json_result) {

    $disease_found_count++;

    # decode the JSON format to text
    my $query_result = decode_json $json_result;

    # get the name of the disease
    my $name = $query_result->{'name'};

    # output the disease details to file
    print DISEASE_FILE "$doid\nDisease:$name\n";

    # get the definition string
    my $def_string = $query_result->{'definition'};

    # check if a definition is supplied
    if (defined $def_string) {

       # separate the string into defintion part & url part
       # based on double quote delimiter
       my @def_parts = split ("\"",$def_string);

       # actual definition should be second element in array
       my $definition = $def_parts[1];

       # add the definition to file
       print DISEASE_FILE "Definition:$definition\n\n";
    } else {
       # indicate no definition is given
       print DISEASE_FILE "Definition: None given\n\n";
    }

  } else { # Disease not found

    $disease_not_found_count++;
    print STDERR "ERROR: Disease Ontology ID $doid not found\n";
    print NO_DISEASE_FILE "$doid\n";

  } # end if disease ID found

  # print message for every 20th disease processed
  print "Disease entries processed:$disease_count\n" unless ($disease_count % 20);

} # end foreach Disease Ontology identifier

close (DISEASE_FILE);
close (NO_DISEASE_FILE);

print "Total disease IDs: $disease_count\n";
print "Disease IDs found: $disease_found_count, output file $disease_filename\n";
print "Disease IDs not found: $disease_not_found_count, output file $no_disease_filename\n";

