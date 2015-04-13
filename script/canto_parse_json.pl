#!/usr/bin/perl
use strict;
use warnings;

#use LWP::Simple;
use JSON;
use DBI;
use JSON::Parse 'json_file_to_perl';
use Data::Compare;

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

# create a hash of hashes to store the interactions
# where the key will be an identifier for each interaction
# and the value will be a hash of the values that identify each unique interaction
my %interaction_profiles;
my $interaction_count = 0;

# iterate through the array of annotations
foreach my $annot_index (0 .. $#annotations) {

   # create temporary hash to identify the interaction to which the annotation should belong
   # this should include all the info that specifies the unqiue interaction,
   # such as UniProt accessions, interaction partners, and annotation extensions
   my %annotation_interaction_hash = ();

   my %annotation = %{ $text_response->{'curation_sessions'}{$session_id}{'annotations'}[$annot_index] };

   my @gene_organism_list = keys $annotation{'genes'};
   my $gene_organism_name = $gene_organism_list[0];

   my $gene_id = $annotation{'genes'}{$gene_organism_name}{'uniquename'};

   # add uniprot accession to help identify the correct interaction
   $annotation_interaction_hash{'pathogen_gene_id'} = $gene_id;

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
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_taxon_id'} = $host;
     }

     # if the annotation extension begins with 'occurs_in', then assign value between brackets to tissue
     if ($annot_ext =~ /^occurs_in/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $tissue = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_tissue'} = $tissue;
     }

     # if the annotation extension begins with 'interaction_partner',
     # then assign value between brackets to interaction partner
     # (TODO: THERE CAN BE MULTIPLE INTERACTION PARTNERS)
     if ($annot_ext =~ /^interaction_partner/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $interaction_partner = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'interaction_partner'} = $interaction_partner;
     }

   }


   # declare flag to indicate if the interaction associated with
   # the current annotation has been found
   my $interaction_found = 0;

   # iterate through each of the existing interactions to find out
   # if the current annotation profile matches an existing interaction
   # if so, then add the corresponding interaction id to the annotation
   # otherwise, create a new interaction profile for the annotation
   foreach my $interaction_id (keys %interaction_profiles) {

     # get the interaction details
     my %interaction_hash = %{ $interaction_profiles{$interaction_id} };

     # compare the current annotation profile to the interaction profile
     # NOTE USE OF REFERENCE OPERATOR (BACKSLASH) TO PASS THE REFERENCES OF HASHES TO THE COMPARE OPERATOR
     if ( Compare(\%annotation_interaction_hash,\%interaction_hash) ) {
        # if the annotation interaction profile matches an existing interaction,
        # then that annotation should be assigned the associated iteraction ID
        $annotations[$annot_index]{"interaction_id"} = $interaction_id;
        $interaction_found = 1;
        last;
     }

   }

   # if a matching interaction profile was not found,
   # then a new interaction profile should be created based on the annotation profile,
   # with the new interaction ID assigned to the annotation
   if (not $interaction_found) {
      $interaction_count++;
      $interaction_profiles{$interaction_count} =  { %annotation_interaction_hash } ;
      $annotations[$annot_index]{"interaction_id"} = $interaction_count;
   }

}


# iterate through each interaction and retrieve all of the annotations
# associated with the interaction
foreach my $int_id (1 .. $interaction_count) {

  print "\nINTERACTION:$int_id\n";

  my $first_annot_of_int_flag = 1;

  # iterate through all annotations, to find out which ones
  # are associated with the current interaction
  foreach my $annot_index (0 .. $#annotations) {

    # get the interaction ID assocaiated with the annotation
    my $annot_int_id = $annotations[$annot_index]{"interaction_id"};

    # check if this interaction ID matches the current interaction
    if ($annot_int_id == $int_id) {

      # get the hash of annotation details
      my %annotation = %{ $annotations[$annot_index] };

      # since most of the information about an interaction is repeated
      # in all of the annotations, we can just retrieve the relevant data
      # from the first annotation
      if ($first_annot_of_int_flag) {

         $first_annot_of_int_flag = 0;

         my $pubmed_id = $annotation{'publication'};
         print "PubMed ID: $pubmed_id\n";

         my @gene_organism_list = keys $annotation{'genes'};
         my $gene_organism_name = $gene_organism_list[0];
         my $gene_id = $annotation{'genes'}{$gene_organism_name}{'uniquename'};
         print "Annot Attr Path Gene:$gene_id\n";

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

         } # end foreach annotation extension

         print "Host ID: $host\n";
         print "Tissue ID: $tissue\n";

         if (defined $interaction_partner) {
           print "Interaction Partner: $interaction_partner\n";
         }

         my $curator_name = $annotation{'curator'}{'name'};
         print "Curator: $curator_name\n";

         my $pathogen_species = $annotation{'genes'}{$gene_organism_name}{'organism'};
         print "Pathogen Species: $pathogen_species\n";

         my $creation_date = $annotation{'creation_date'};
         print "Creation Date:$creation_date\n";

      } # end if first annotation

      print "Annot Attr Term:$annotations[$annot_index]{'type'}\n";
      print "Annot Attr Value:$annotations[$annot_index]{'term'}\n";

    } # end if annotation belongs to interaction

  } # end foreach annotation

} # end foreach interaction


print "\nTotal number of interactions:$interaction_count\n";
print "Total number of annotations :".($#annotations+1)."\n";

close (JSON_OUTPUT_FILE);

#$sql_stmt->finish() or die "Failed to finish SQL statement\n";
#$db_conn->disconnect() or die "Failed to disconnect database\n";

#print "Unique PubMed IDs:$article_count\n";
#print "PubMed articles found: $article_found_count, output file $pubmed_filename\n";
#print "PubMed articles not found: $article_not_found_count, output file $no_pubmed_filename\n";

