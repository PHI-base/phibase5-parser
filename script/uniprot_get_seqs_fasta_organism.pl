use strict;
use warnings;
use LWP::UserAgent;

my $taxon = $ARGV[0]; # Taxonomy identifier of organism.

my $query = "http://www.uniprot.org/uniprot/?query=organism:$taxon&format=fasta&include=yes";
my $file = '../output/taxon/aa_sequences_taxon_'.$taxon.'.fasta';

my $contact = 'alistair.irvine@rothamsted.ac.uk'; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
my $response = $agent->mirror($query, $file);

if ($response->is_success) {
  my $results = $response->header('X-Total-Results');
  my $release = $response->header('X-UniProt-Release');
  my $date = sprintf("%4d-%02d-%02d", HTTP::Date::parse_date($response->header('Last-Modified')));
  print "Downloaded $results entries of UniProt release $release ($date) to file $file\n";
}
elsif ($response->code == HTTP::Status::RC_NOT_MODIFIED) {
  print "Data for taxon $taxon is up-to-date.\n";
}
else {
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
}

