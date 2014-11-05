use strict;
use warnings;

use LWP::Simple;
use JSON;
use DBI;

use phibase_subroutines qw(connect_to_phibase); 

my $db_conn = connect_to_phibase();

print "Connection established.\n";

# run query to get all PubMed IDs from PHI-base
my $sql_query = qq(SELECT DISTINCT pubmed_id FROM interaction_literature;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# iterator to count rows
my $row_count = 0;

# open output file
open (PUBMED_FILE,"> ./output/pubmed_details.tsv") or die "Error opening output file\n";
print "Printing PubMed details to output file /output/pubmed_details.tsv...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $row_count++;
  foreach my $pubmed_id (@row) {

     $pubmed_id =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $pubmed_id =~ s/\s+$//; # remove blank space from end of UniProt IDs
     print PUBMED_FILE "$pubmed_id\t";

     my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=EXT_ID:$pubmed_id&format=json";
     my $json_response = get $url;
     die "Error getting $url" unless defined $json_response;

     my $text_response = decode_json($json_response);
     my $authors = $text_response->{'resultList'}{'result'}[0]{'authorString'};
     my $year    = $text_response->{'resultList'}{'result'}[0]{'pubYear'};
     my $title   = $text_response->{'resultList'}{'result'}[0]{'title'};
     my $journal = $text_response->{'resultList'}{'result'}[0]{'journalTitle'};
     my $volume  = $text_response->{'resultList'}{'result'}[0]{'journalVolume'};
     my $issue   = $text_response->{'resultList'}{'result'}[0]{'issue'};
     my $pages   = $text_response->{'resultList'}{'result'}[0]{'pageInfo'};
     my $doi     = $text_response->{'resultList'}{'result'}[0]{'doi'};

     print PUBMED_FILE "$authors ($year). \"$title\" $journal $volume($issue): $pages. $doi.\n\n";
  }
}
print "Unique PubMed IDs:$row_count\n";
close (PUBMED_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

