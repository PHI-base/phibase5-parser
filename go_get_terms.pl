#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use XML::Twig;
use DBI;

use phibase_subroutines qw(connect_to_phibase); # load PHI-base functions
my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all GO IDs from PHI-base
my $sql_query = qq(SELECT DISTINCT go_id FROM interaction_go_term;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $go_term_count = 0;
my $go_found_count = 0;
my $go_not_found_count = 0;

# open output files
my $go_filename = './output/go_terms.tsv';
my $no_go_filename = './output/go_terms_not_found.tsv';
open (GO_FILE,"> $go_filename") or die "Error opening output file\n";
open (NO_GO_FILE,"> $no_go_filename") or die "Error opening output file\n";
print "Printing Gene Ontology terms to output file $go_filename...\n";

# iterate through each GO ID
while (my @row = $sql_stmt->fetchrow_array()) {

  $go_term_count++;

  my $go_id = shift @row;

  $go_id =~ s/^\s+//; # remove blank space from start of taxon IDs
  $go_id =~ s/\s+$//; # remove blank space from end of taxon IDs

  my $query = "http://www.ebi.ac.uk/QuickGO/GTerm?id=$go_id&format=oboxml";
  my $xml_response = get $query;
  #die "Error getting $query" unless defined $xml_response;

  # use XML twig to parse the XML data
  my $xml_twig = XML::Twig->new();

  if (defined $xml_response) {

     # parse the XML data to get the GO term name
     $xml_twig->parse($xml_response);
     my $go = $xml_twig->root->first_child('term')->field('name');

     # print to output file
     $go_found_count++;
     print GO_FILE "$go_id\t$go\n";

  } else {
  
     $go_not_found_count++;
     print STDERR "ERROR: Gene Ontology term not found for $go_id\n";
     print NO_GO_FILE "$go_id\n";

  }

} # end while GO terms

close (GO_FILE);
close (NO_GO_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Unique GO IDs: $go_term_count\n";
print "GO IDs found: $go_found_count, output file $go_filename\n";
print "GO IDs not found: $go_not_found_count, output file $no_go_filename\n";

