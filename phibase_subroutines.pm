package phibase_subroutines;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(connect_to_phibase query_uniprot);

sub connect_to_phibase
{
  print "Connecting to PHI-base...\n";

  # define credentials for phibase database
  my $db_name = "phibase";
  my $db_host = "localhost";
  my $db_user = "postgres";
  my $db_pw = "";

  # connect to database
  my $db_conn = DBI->connect("DBI:Pg:dbname=$db_name;host=$db_host","$db_user","$db_pw");
}

sub query_uniprot
{
  my $query = shift;  # uniprot RESTful query URL should be first argument
  
  my $contact = 'alistair.irvine@rothamsted.ac.uk'; # email address requested by UniProt in case of problems.
  my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
  push @{$agent->requests_redirectable}, 'POST';
  my $response = $agent->post("$query");

  while (my $wait = $response->header('Retry-After')) {
    print STDERR "Waiting ($wait)...\n";
    sleep $wait;
    $response = $agent->get($response->base);
  }
  $response->is_success ? $response->content : die 'Failed, got '.$response->status_line.' for '.$response->request->uri."\n";
}

return 1;  # return true to calling function

