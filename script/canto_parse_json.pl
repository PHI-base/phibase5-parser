#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use JSON;
use DBI;
use JSON::Parse 'json_file_to_perl';

#use phibase_subroutines qw(connect_to_phibase); 

#my $db_conn = connect_to_phibase();

# run query to get all PubMed IDs from PHI-base
#my $sql_query = qq(SELECT DISTINCT pubmed_id FROM interaction_literature;);
#my $sql_stmt = $db_conn->prepare($sql_query);
#my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
#my $article_count = 0;
#my $article_found_count = 0;
#my $article_not_found_count = 0;

# open output files
my $json_filename = '../input/canto/approved_annotation_2015-03-06.json';
my $output_filename = '../output/json_output.tsv';
open (JSON_OUTPUT_FILE,"> $output_filename") or die "Error opening output file\n";
print "Printing JSON details to output file $json_filename...\n";

my $text_response = json_file_to_perl($json_filename);

my $creation_date = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'creation_date'};
my $curator_name = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'curator'}{'name'};
my $gene_id = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'genes'}{'Fusarium oxysporum FOXG_00076'}{'uniquename'};
my $pathogen_species = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'genes'}{'Fusarium oxysporum FOXG_00076'}{'organism'};
my $pubmed_id = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'publication'};
my $ontology = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'type'};
my $ontology_term = $text_response->{'curation_sessions'}{'36e8ad26889b5b01'}{'annotations'}[0]{'term'};
print "PubMed ID: $pubmed_id\n";
print "Pathogen Species: $pathogen_species\n";
print "Pathogen Gene ID: $gene_id\n";
print "Ontology: $ontology\n";
print "Ontology Term: $ontology_term\n";
print "Curator: $curator_name\n";
print "Creation Date: $creation_date\n";
close (JSON_OUTPUT_FILE);

#$sql_stmt->finish() or die "Failed to finish SQL statement\n";
#$db_conn->disconnect() or die "Failed to disconnect database\n";

#print "Unique PubMed IDs:$article_count\n";
#print "PubMed articles found: $article_found_count, output file $pubmed_filename\n";
#print "PubMed articles not found: $article_not_found_count, output file $no_pubmed_filename\n";

