#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use JSON;
use DBI;
use OBO::Parser::OBOParser;

# open output file
my $phenotype_filename = './output/phenotype_term_details.tsv';
open (PHENOTYPE_FILE,"> $phenotype_filename") or die "Error opening output file\n";
print "Printing phenotype term details to output file $phenotype_filename...\n";

# counter for statistics
my $phenotype_term_count = 0;

# load and parse the ontology file
my $obo_parser = OBO::Parser::OBOParser->new;
my $ontology = $obo_parser->work("phenotype_outcome.obo");

# get the ontology term, based on the identifier
#my $ontology_term = $ontology->get_term_by_id("ID:0000001");

my @ontology_terms = @{$ontology->get_terms()}; # get all the terms in the ontology

foreach my $term (@ontology_terms) {

  # get the details of the ontology term
  my $term_id = $term->id;
  my $term_name = $term->name;
  my $term_definition = $term->def->text();

  # increment counter
  $phenotype_term_count++;

  # print the details of the term
  print PHENOTYPE_FILE "$term_id\n$term_name\n$term_definition\n\n";

}

close (PHENOTYPE_FILE);

print "Number of phenotype terms:$phenotype_term_count\n";
print "Phenotype term details output file $phenotype_filename\n";

