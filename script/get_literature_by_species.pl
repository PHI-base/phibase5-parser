#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module
use LWP::Simple;
use JSON;
use XML::Twig;

# load PHI-base functions
use phibase_subroutines 
  qw(connect_to_phibase 
     query_uniprot 
     ontology_mapping
    );

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open output file
my $db_data_filename = '../output/literature_by_species.tsv';  
open (DATABASE_DATA_FILE, "> $db_data_filename") or die "Error opening output file\n";

# counters to gather statistics
my $path_taxon_count = 0;
my $path_species_count = 0;
my $path_no_species_count = 0;

# array for list of unique species
my @unique_species;

# print the headers for the output file
print DATABASE_DATA_FILE 
"Pathogen Taxon ID\tPathogen Taxon Name\tPathogen Species Taxon ID\tPathogen Species Name\tPubMed IDs\n";

# first, get NCBI taxon ID for all pathogen taxon
my $sql_stmt = qq(SELECT DISTINCT ncbi_taxon_id FROM pathogen_gene);

my $sql_result = $db_conn->prepare($sql_stmt);
$sql_result->execute() or die $DBI::errstr;

my $pathogen_taxa_count = 0;

print "Parsing pathogen taxa...\n";

while (my @row = $sql_result->fetchrow_array()) {

  $pathogen_taxa_count++;

  my $path_taxon_id = shift @row;

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
  print DATABASE_DATA_FILE "$path_taxon_id\t$path_taxon_name\t";

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
  print DATABASE_DATA_FILE "$species_taxon_id\t$species_taxon_name\t";

  # get the PubMed IDs for the current pathogen taxon 
  my $sql_stmt2 = qq(SELECT DISTINCT interaction_literature.pubmed_id
                    FROM interaction,
                         interaction_literature,
                         interaction_pathogen_gene_mutant,
                         pathogen_gene_mutant
                   WHERE pathogen_gene_mutant.ncbi_taxon_id =  $path_taxon_id
                     AND interaction_pathogen_gene_mutant.pathogen_gene_mutant_id = pathogen_gene_mutant.id
                     AND interaction_pathogen_gene_mutant.interaction_id = interaction.id
                     AND interaction.id = interaction_literature.interaction_id
                ;);

  my $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for literature
  my $pubmed_output_string = "";

  # since there may be multiple PubMed articles,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (my @row2 = $sql_result2->fetchrow_array()) {

    my $pubmed_id = shift @row2;
    $pubmed_output_string .= "$pubmed_id;";

  }

  # remove the final semi-colon from end of the string
  $pubmed_output_string =~ s/;$//;
  # print the list of pubmed ids to file
  print DATABASE_DATA_FILE "$pubmed_output_string\n";

}

close (DATABASE_DATA_FILE);

$sql_result->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "\nProcess completed successfully.\n";
print "Total pathogen taxa:$pathogen_taxa_count\n";
print "Tab-separated file of pathogen species and PubMed IDs: $db_data_filename\n\n";

