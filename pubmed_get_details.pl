#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use JSON;
use DBI;

use phibase_subroutines qw(connect_to_phibase); 

my $db_conn = connect_to_phibase();

# run query to get all PubMed IDs from PHI-base
my $sql_query = qq(SELECT DISTINCT pubmed_id FROM interaction_literature;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $article_count = 0;
my $article_found_count = 0;
my $article_not_found_count = 0;

# open output files
my $pubmed_filename = './output/pubmed_details.tsv';
my $no_pubmed_filename = './output/pubmed_not_found.tsv';
open (PUBMED_FILE,"> $pubmed_filename") or die "Error opening output file\n";
open (NO_PUBMED_FILE,"> $no_pubmed_filename") or die "Error opening output file\n";
print "Printing PubMed details to output file $pubmed_filename...\n";

# iterate through each PubMed ID
while (my @row = $sql_stmt->fetchrow_array()) {

  $article_count++;

  my $pubmed_id = shift @row; # get ID from row
  $pubmed_id =~ s/^\s+//; # remove blank space from start of PubMed IDs
  $pubmed_id =~ s/\s+$//; # remove blank space from end of PubMed IDs

  # run REST query and get JSON response
  my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=EXT_ID:$pubmed_id&format=json";
  my $json_response = get $url;
  my $text_response = decode_json($json_response);

  # parse each of the relevant parameters from the JSON text
  my $authors = $text_response->{'resultList'}{'result'}[0]{'authorString'};
  my $year    = $text_response->{'resultList'}{'result'}[0]{'pubYear'};
  my $title   = $text_response->{'resultList'}{'result'}[0]{'title'};
  my $journal = $text_response->{'resultList'}{'result'}[0]{'journalTitle'};
  my $volume  = $text_response->{'resultList'}{'result'}[0]{'journalVolume'};
  my $issue   = $text_response->{'resultList'}{'result'}[0]{'issue'};
  my $pages   = $text_response->{'resultList'}{'result'}[0]{'pageInfo'};
  my $doi     = $text_response->{'resultList'}{'result'}[0]{'doi'};

  # if title is empty or undefined, then assume the article has not been found
  if (defined $title and $title ne "") {
     $article_found_count++;
     # print article details in citation format
     # note that warnings about non-ascii characters is suppressed,
     # but these characters may not display as desired.
     { no warnings; print PUBMED_FILE "$pubmed_id\t$authors ($year). \"$title\" $journal $volume($issue): $pages. $doi.\n" };
  } else { # article not found
     $article_not_found_count++;
     print STDERR "ERROR:PubMed ID $pubmed_id not found\n";
     print NO_PUBMED_FILE "$pubmed_id\n";
  }

} # end while PubMed IDs

close (PUBMED_FILE);
close (NO_PUBMED_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Unique PubMed IDs:$article_count\n";
print "PubMed articles found: $article_found_count, output file $pubmed_filename\n";
print "PubMed articles not found: $article_not_found_count, output file $no_pubmed_filename\n";

