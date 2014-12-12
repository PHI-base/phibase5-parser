#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use List::MoreUtils qw(uniq);
use DBI;

use phibase_subroutines qw(connect_to_phibase); # load PHI-base functions
my $db_conn = connect_to_phibase(); # connect to PHI-base database

# run query to get all UniProt accessions from PHI-base
my $sql_query = qq(SELECT DISTINCT uniprot_accession FROM pathogen_gene_mutant;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# counters to gather statistics
my $uniprot_acc_count = 0;
my $go_function_count = 0;
my $uniprot_no_go_function_count = 0;

# open output files
my $go_filename = '../output/go_function_annots.tsv';
my $no_go_filename = '../error/go_function_annots_not_found.tsv';
open (GO_FILE,"> $go_filename") or die "Error opening output file\n";
open (NO_GO_FILE,"> $no_go_filename") or die "Error opening output file\n";
print "Printing Gene Ontology function annotations to output file $go_filename...\n";

# iterate through each UniProt entry
while (my @row = $sql_stmt->fetchrow_array()) {

  $uniprot_acc_count++;

  my $uniprot_acc = shift @row;

  $uniprot_acc =~ s/^\s+//; # remove blank space from start of taxon IDs
  $uniprot_acc =~ s/\s+$//; # remove blank space from end of taxon IDs

  # run QuickGO REST query, using UniProt accession as parameter,
  # returning only the Molecular Function annotations
  my $query = "http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&protein=$uniprot_acc&aspect=F&col=goID,goName,evidence";
  my $response = get $query;
    die "Error getting $query" unless defined $response;

  # each go annotation will be a separate line, starting with a header
  my @go_annotations = split("\n",$response);
  # remove the header from the annotations
  shift @go_annotations;
  
  # since some GO annotations are repeated (from different sources),
  # we want to get only the unique values
  @go_annotations = uniq(@go_annotations);

  # if any GO terms exist, print them to output file
  if (@go_annotations) {

     print GO_FILE "\nUniProt accession:$uniprot_acc\n";
     print GO_FILE "GO ID\tGo Term\tGO Evidence Code\n";
     
     foreach my $go_annotation (@go_annotations) {

        $go_function_count++;

        # separate the go_annotation into individual fields,
        # which are tab-separated
        my @go_fields = split("\t",$go_annotation);

	# print each field of the GO annotation
        foreach my $go_field (@go_fields) {
          print GO_FILE "$go_field\t";
        }
        print GO_FILE "\n";

     } # foreach go annotation

  } else { # no go annotations for UniProt entry 
     $uniprot_no_go_function_count++;
     print NO_GO_FILE "$uniprot_acc\n";
  }

  # print message for every 20th UniProt entry processed
  print "UniProt entries processed:$uniprot_acc_count\n" unless ($uniprot_acc_count % 20);

} # end while UniProt entries

close (GO_FILE);
close (NO_GO_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Unique UniProt accessions: $uniprot_acc_count\n";
print "Total GO Molecular Function annotations found: $go_function_count, output file $go_filename\n";
print "UniProt entries without GO Molecular Function annotation: $uniprot_no_go_function_count, output file $no_go_filename\n";

