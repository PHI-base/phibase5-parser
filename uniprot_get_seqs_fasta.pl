use strict;
use warnings;
use LWP::UserAgent;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all UniProt accessions from PHI-base
my $sql_query = qq(SELECT DISTINCT uniprot_accession FROM pathogen_gene_mutant;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $uniprot_count = 0;
my $fasta_seq_count = 0;
my $no_fasta_seq_count = 0;

# open output file
my $fasta_filename = './output/uniprot_seq.fasta';
my $no_fasta_filename = './output/uniprot_no_seq.tsv';
open (FASTA_FILE, "> $fasta_filename") or die "Error opening output file\n";
open (NO_FASTA_FILE, "> $no_fasta_filename") or die "Error opening output file\n";

print "Printing to FASTA output file $fasta_filename...\n";

# iterate through each UniProt ID
while (my @row = $sql_stmt->fetchrow_array()) {
  $uniprot_count++;
  foreach my $uniprot_acc (@row) {

     $uniprot_acc =~ s/^\s+//; # remove blank space from start of UniProt IDs
     $uniprot_acc =~ s/\s+$//; # remove blank space from end of UniProt IDs

     my $query = "http://www.uniprot.org/uniprot/$uniprot_acc.fasta";

     my $contact = 'alistair.irvine@rothamsted.ac.uk'; # email address requested by UniProt in case of problems.
     my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
     push @{$agent->requests_redirectable}, 'POST';
     my $response = $agent->post("$query");

     while (my $wait = $response->header('Retry-After')) {
       print STDERR "Waiting ($wait)...\n";
       sleep $wait;
       $response = $agent->get($response->base);
     }

     if ($response->is_success) {
        $fasta_seq_count++;
        print FASTA_FILE $response->content."\n";
     }  else {
        $no_fasta_seq_count++;
        print STDERR "Unable to get FASTA sequence for UniProt ID $uniprot_acc, tried URL ".$response->request->uri."\n";
        print NO_FASTA_FILE "$uniprot_acc\n";
        next;
     }

  } # end foreach UniProt entry

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_count\n" unless ($uniprot_count % 20);

} # end while rows

close (FASTA_FILE);
close (NO_FASTA_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

# print statistics and output filenames
print "Process completed successfully.\n";
print "Unique UniProt accessions:$uniprot_count\n";
print "FASTA sequences retrieved:$fasta_seq_count, available in output file $fasta_filename\n";
print "UniProt accessions without a FASTA sequence:$no_fasta_seq_count, availalbe in file $no_fasta_filename\n";

