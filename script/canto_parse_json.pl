#!/usr/bin/perl
use strict;
use warnings;

#use LWP::Simple;
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
my $json_filename = '../input/canto/canto_phibase_extensions.json';
#my $json_filename = '../input/canto/approved_annotation_2015-03-06.json';
my $output_filename = '../output/json_output.tsv';
open (JSON_OUTPUT_FILE,"> $output_filename") or die "Error opening output file\n";
print "Printing JSON details to output file $output_filename...\n";

my $text_response = json_file_to_perl($json_filename);

my @session_ids = keys $text_response->{'curation_sessions'};
my $session_id = $session_ids[0];
print "Session ID: $session_id\n";

# annotations are an array of hashes
my @annotations = @{ $text_response->{'curation_sessions'}{$session_id}{'annotations'} };

# iterate through the array of annotations
foreach my $annot_index (0 .. $#annotations) {

   my %annotation = %{ $text_response->{'curation_sessions'}{$session_id}{'annotations'}[$annot_index] };

   my $creation_date = $annotation{'creation_date'};
   my $curator_name = $annotation{'curator'}{'name'};

   my @gene_organism_list = keys $annotation{'genes'};
   my $gene_organism_name = $gene_organism_list[0];

   my $gene_id = $annotation{'genes'}{$gene_organism_name}{'uniquename'};
   my $pathogen_species = $annotation{'genes'}{$gene_organism_name}{'organism'};

   my $pubmed_id = $annotation{'publication'};

   my $annot_extension_list = $annotation{'annotation_extension'};
   my @annot_extensions = split(/,/,$annot_extension_list);

   my $host;
   my $tissue;
   my $interaction_partner;

   # identify each annotation extension in the list
   foreach my $annot_ext (@annot_extensions) {

     # if the annotation extension begins with 'pathogen_of', then assign value between brackets to host
     if ($annot_ext =~ /^pathogen_of/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $host = $annot_ext[1];
     }

     # if the annotation extension begins with 'occurs_in', then assign value between brackets to tissue
     if ($annot_ext =~ /^occurs_in/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $tissue = $annot_ext[1];
     }

     # if the annotation extension begins with 'interaction_partner',
     # then assign value between brackets to interaction partner
     # (TODO: THERE CAN BE MULTIPLE INTERACTION PARTNERS)
     if ($annot_ext =~ /^interaction_partner/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $interaction_partner = $annot_ext[1];
     }

   }

   my $ontology = $annotation{'type'};
   my $ontology_term = $annotation{'term'};

   print "\nANNOTATION ".($annot_index+1)."\n";
   print "PubMed ID: $pubmed_id\n";
   print "Pathogen Species: $pathogen_species\n";
   print "Pathogen Gene ID: $gene_id\n";
   print "Ontology: $ontology\n";
   print "Ontology Term: $ontology_term\n";
   print "Curator: $curator_name\n";
   print "Creation Date: $creation_date\n";
   print "Annotation Extensions: $annot_extension_list\n";
   print "Host ID: $host\n";
   print "Tissue ID: $tissue\n";

   if (defined $interaction_partner) {
     print "Interaction Partner: $interaction_partner\n";
   }

}

close (JSON_OUTPUT_FILE);

#$sql_stmt->finish() or die "Failed to finish SQL statement\n";
#$db_conn->disconnect() or die "Failed to disconnect database\n";

#print "Unique PubMed IDs:$article_count\n";
#print "PubMed articles found: $article_found_count, output file $pubmed_filename\n";
#print "PubMed articles not found: $article_not_found_count, output file $no_pubmed_filename\n";

