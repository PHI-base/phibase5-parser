package phibase_subroutines;

use strict;
use warnings;
use Exporter qw(import);
use OBO::Parser::OBOParser;

our @EXPORT_OK = qw(connect_to_phibase query_uniprot print_ontology_terms ontology_mapping);

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

sub print_ontology_terms
{
  my $ontology_name = shift;
  my $obo_filename = shift;

  # open output file
  my $term_filename = '../output/'.$ontology_name.'_term_details.tsv';
  open (TERM_FILE,"> $term_filename") or die "Error opening output file\n";

  # counter for statistics
  my $ontology_term_count = 0;
  my $term_with_def_count = 0;
  my $term_without_def_count = 0;

  # load and parse the ontology file
  my $obo_parser = OBO::Parser::OBOParser->new;
  my $ontology = $obo_parser->work($obo_filename);

  # get the ontology term, based on the identifier
  #my $ontology_term = $ontology->get_term_by_id("ID:0000001");

  my @ontology_terms = @{$ontology->get_terms()}; # get all the terms in the ontology

  foreach my $term (@ontology_terms) {

     # increment counter
     $ontology_term_count++;

     # get the ID and name of the ontology term
     my $term_id = $term->id;
     my $term_name = $term->name;

     # print to output file
     print TERM_FILE "$term_id\nTerm:$term_name\n";

     # get definition of the term
     my $term_definition = $term->def->text();

     # check if a definition was given
     if ($term_definition) {
        # print the definition
        $term_with_def_count++;
        print TERM_FILE "Definition:$term_definition\n\n";
     } else {
        $term_without_def_count++;
        print TERM_FILE "Definition:None given\n\n";
     }
  }

  close (TERM_FILE);

  print "Number of $ontology_name terms:$ontology_term_count\n";
  print "Number of terms defined:$term_with_def_count\n";
  print "Number of terms not defined:$term_without_def_count\n";
  print "$ontology_name term details output file $term_filename\n";

} # end print_ontology_terms sub


sub ontology_mapping
{
  my $obo_filename = shift;

  # open output file
  my %ontology_hash;

  # load and parse the ontology file
  my $obo_parser = OBO::Parser::OBOParser->new;
  my $ontology = $obo_parser->work($obo_filename);

  my @ontology_terms = @{$ontology->get_terms()}; # get all the terms in the ontology

  foreach my $term (@ontology_terms) {

     # get the ID and name of the ontology term
     my $term_id = $term->id;
     my $term_name = $term->name;

     # create hash to map term name to the
     # corresponding term ID of the ontology
     $ontology_hash{$term_name} = $term_id;
     
  }

  return %ontology_hash;

}


return 1;  # return true to calling function

