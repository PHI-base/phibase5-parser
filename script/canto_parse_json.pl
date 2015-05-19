#!/usr/bin/perl
use strict;
use warnings;

#use LWP::Simple;
use JSON;
use DBI;
use JSON::Parse 'json_file_to_perl';
use Data::Compare;

use phibase_subroutines qw(connect_to_phibase); 

my $db_conn = connect_to_phibase();

# counters to gather statistics
my $annotation_count = 0;
my $interaction_count = 0;
my $go_annotation_count = 0;
my $host_response_count = 0;
my $phenotype_outcome_count = 0;
my $disease_count = 0;
# CURRENTLY UNUSED
#my $defect_count = 0;
#my $anti_infective_count = 0;
#my $inducer_count = 0;
#my $exp_spec_count = 0;

my $json_filename = '../input/canto/paper_with_double_triple_mutants.json';
#my $json_filename = '../input/canto/canto_with_go.json';
#my $json_filename = '../input/canto/canto_triple_mutant.json';
#my $json_filename = '../input/canto/canto_phibase_extensions.json';
#my $json_filename = '../input/canto/approved_annotation_2015-03-06.json';

#my $output_filename = '../output/json_output.tsv';
#open (JSON_OUTPUT_FILE,"> $output_filename") or die "Error opening output file\n";
#print "Printing JSON details to output file $output_filename...\n";

my $text_response = json_file_to_perl($json_filename);

my @session_ids = keys $text_response->{'curation_sessions'};
my $session_id = $session_ids[0];
print "Session ID: $session_id\n";

# get data about the approver of the session
my $approver_name = $text_response->{'curation_sessions'}{$session_id}{'metadata'}{'approver_name'};
my $approver_email = $text_response->{'curation_sessions'}{$session_id}{'metadata'}{'approver_email'};
print "Session Approver: $approver_name, $approver_email\n";

# get the pathogen taxon for the session (assumes only one pathogen for the entire session - use the first one available)
my @pathogen_taxon_ids = keys $text_response->{'curation_sessions'}{$session_id}{'organisms'};
my $pathogen_taxon_id = $pathogen_taxon_ids[0];
print "Pathogen Taxon ID: $pathogen_taxon_id\n";

# annotations are an array of hashes
my @annotations = @{ $text_response->{'curation_sessions'}{$session_id}{'annotations'} };

# create a hash of hashes to store the interactions
# where the key will be an identifier for each interaction
# and the value will be a hash of the values that identify each unique interaction
my %interaction_profiles;

# Read in the experiment specification ontology
# so that the identifier can be matched from the
# available exp specification strring (from the evidence code)
print "Reading ontology files...\n";
my $obo_parser = OBO::Parser::OBOParser->new;
my $exp_spec_ontology = $obo_parser->work("../ontology/phibase/experiment_specification.obo");

# iterate through the array of annotations
foreach my $annot_index (0 .. $#annotations) {

   $annotation_count++;

   # create temporary hash to identify the interaction to which the annotation should belong
   # this should include all the info that specifies the unqiue interaction,
   # such as UniProt accessions, interaction partners, and annotation extensions
   my %annotation_interaction_hash = ();

   # array to hold all partner gene uniprot ids
   # (for multiple gene interactions)
   my @interaction_partners;

   my %annotation = %{ $text_response->{'curation_sessions'}{$session_id}{'annotations'}[$annot_index] };

   my @gene_organism_list = keys $annotation{'genes'};
   my $gene_organism_name = $gene_organism_list[0];

   my $uniprot_acc = $annotation{'genes'}{$gene_organism_name}{'uniquename'};

   # add uniprot accession to help identify the correct interaction
   $annotation_interaction_hash{'pathogen_uniprot_acc'} = $uniprot_acc;

   my $annot_extension_list = $annotation{'annotation_extension'};
   my @annot_extensions = split(/,/,$annot_extension_list);

   my $host_taxon;
   my $tissue;
   my $interaction_partner;

   # identify each annotation extension in the list
   foreach my $annot_ext (@annot_extensions) {

     # if the annotation extension begins with 'pathogen_of', then assign value between brackets to host
     if ($annot_ext =~ /^pathogen_of/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $host_taxon = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_taxon_id'} = $host_taxon;
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
     # (note that there can be multiple interaction partners)
     if ($annot_ext =~ /^interaction_partner/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $interaction_partner = $annot_ext[1];
        # add current interaction partner to the array of partners
        # (for cases where there are three are more genes involved in interaction)
        push(@interaction_partners,$interaction_partner);
        # add the interaction partners array to help identify the correct interaction
        $annotation_interaction_hash{'interaction_partner'} = @interaction_partners;
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

  # flag indicating first annotation in the interaction
  my $first_annot_of_int_flag = 1;

  # declare variables
  my $uniprot_acc;
  my $pubmed_id;
  my $host_taxon;
  my $host_taxon_id;
  my $tissue;
  my $interaction_partner;
  my $curator_name;
  my $curator_email;
  my $pathogen_species;
  my $creation_date;
  my $pathogen_taxon;
#  my $pathogen_taxon_id;
  my $phenotype_outcome;
  my $disease;
  my $go_id;
  my $go_evid_code;
  my $host_response_id;
  my $interaction_id;
  my $experiment_spec;

  # declare array to hold all pathogen gene
  # involved in the interaction
  # (useful for multiple gene interactions)
  my @pathogen_genes;

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

         $pubmed_id = $annotation{'publication'};
         print "PubMed Article: $pubmed_id\n";

=pod
         # the first annotation type must be the pathogen taxon ID
         my $annot_type = $annotations[$annot_index]{'type'};
         my $annot_value = $annotations[$annot_index]{'term'};

         # extract the pathogen taxon ID from the annotation
         if ($annot_type eq 'pathogen_taxon') {
            $pathogen_taxon = $annot_value;
            print "Pathogen Taxon: $pathogen_taxon\n";
            # the pathogen taxon is split between the name and value, separated by colon
            # need to extract the taxon value
            my @pathogen_taxon_parts = split(/[:]/,$pathogen_taxon);
            $pathogen_taxon_id = $pathogen_taxon_parts[1];
         }
=cut

         my @gene_organism_list = keys $annotation{'genes'};
         my $gene_organism_name = $gene_organism_list[0];
         $uniprot_acc = $annotation{'genes'}{$gene_organism_name}{'uniquename'};
         # add the uniprot accession for this gene to the array
         # of genes, which includes all interaction partner genes
         # (for multiple gene interaction)
         push(@pathogen_genes,$uniprot_acc);

         my $annot_extension_list = $annotation{'annotation_extension'};
         my @annot_extensions = split(/,/,$annot_extension_list);

         # identify each annotation extension in the list
         foreach my $annot_ext (@annot_extensions) {

            # if the annotation extension begins with 'pathogen_of', then assign value between brackets to host
            if ($annot_ext =~ /^pathogen_of/) {
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $host_taxon = $annot_ext[1];
            }

            # if the annotation extension begins with 'occurs_in', then assign value between brackets to tissue
            if ($annot_ext =~ /^occurs_in/) {
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $tissue = $annot_ext[1];
            }

            # if the annotation extension begins with 'causes_disease', then assign value between brackets to disease
            if ($annot_ext =~ /^causes_disease/) {
print "Disease Annot Ext:$annot_ext\n";
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $disease = $annot_ext[1];
print "Disease ID:$disease\n";
            }

            # if the annotation extension begins with 'host_responsee', then assign value between brackets to host_response
            if ($annot_ext =~ /^host_response/) {
print "Host Response Annot Ext:$annot_ext\n";
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $host_response_id = $annot_ext[1];
print "Host Response ID:$host_response_id\n";
            }

            # if the annotation extension begins with 'interaction_partner',
            # then assign value between brackets to interaction partner
            # (note that there can multiple interaction partners)
            if ($annot_ext =~ /^interaction_partner/) {

              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $interaction_partner = $annot_ext[1];
	      # UniProt accession is the interaction partner with 'UniProt:' prefix removed
	      my $interaction_partner_id = substr($interaction_partner,8);

	      # add the uniprot accession for this gene to the array
	      # of genes, which includes all interaction partner genes
	      # (for multiple gene interaction) and the original gene
	      push(@pathogen_genes,$interaction_partner_id);

            }


         } # end foreach annotation extension

         foreach my $pathogen_gene (@pathogen_genes) {
           print "Pathogen Gene: $pathogen_gene\n";
         }

         print "Host Taxon: $host_taxon\n";
         print "Tissue: $tissue\n";
         print "Disease: $disease\n";
         print "Host Response: $host_response_id\n";

         $curator_name = $annotation{'curator'}{'name'};
         $curator_email = $annotation{'curator'}{'email'};
         print "Curator: $curator_name, $curator_email\n";
         print "Approver: $approver_name, $approver_email\n";

         $pathogen_species = $annotation{'genes'}{$gene_organism_name}{'organism'};
         print "Pathogen Species: $pathogen_species\n";

         $creation_date = $annotation{'creation_date'};
         print "Creation Date: $creation_date\n";


         # declare array to hold the list of IDs from the pathogen_gene_mutant table 
         my @pathogen_gene_mutant_ids;

	 # insert data into the pathogen_gene table, for each pathogen gene
	 # if it does not exist already (based on combination of taxon id and gene name)
	 # FOR CANTO, USE UNIPROT ID, INSTEAD OF GENE NAME
	 # (NOTE THAT AT THE MOMENT UNIPROT ID IS USED FOR BOTH pathogen_gene AND pathogen_gene_mutant
         # MORE CLARITY IS REQUIRED TO DEFINE THE RELATIONSHIP BETWEEN pathogen_gene
         # (the reference gene from which the amino acid sequence if derived? the wild type gene?)
         # AND pathogen_gene_mutant (the gene mutant on which experiments are performed? the mutated gene?) )
	 foreach my $pathogen_gene_uniprot_acc (@pathogen_genes) {

	    my $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, gene_name) 
				    SELECT $pathogen_taxon_id,'$uniprot_acc'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM pathogen_gene
				      WHERE ncbi_taxon_id = $pathogen_taxon_id
				      AND gene_name = '$pathogen_gene_uniprot_acc'
				   ));

	    my $sql_result2 = $db_conn->prepare($sql_statement2);
	    $sql_result2->execute() or die $DBI::errstr;

	    # get the unique identifier for the inserted pathogen_gene record
	    my $sql_statement4 = qq(SELECT id FROM pathogen_gene
				    WHERE ncbi_taxon_id = $pathogen_taxon_id
				    AND gene_name = '$pathogen_gene_uniprot_acc');

	    my $sql_result4 = $db_conn->prepare($sql_statement4);
	    $sql_result4->execute() or die $DBI::errstr;
	    my @row4 = $sql_result4->fetchrow_array();
	    my $pathogen_gene_id = shift @row4;

	    # insert data into pathogen_gene_mutant table, including foreign key to pathogen_gene table 
	    my $sql_statement3 = qq(INSERT INTO pathogen_gene_mutant (pathogen_gene_id,ncbi_taxon_id,uniprot_accession) 
				     SELECT $pathogen_gene_id,$pathogen_taxon_id,
					    '$pathogen_gene_uniprot_acc'
				     WHERE NOT EXISTS (
				       SELECT 1 FROM pathogen_gene_mutant
				       WHERE pathogen_gene_id = $pathogen_gene_id
				       AND ncbi_taxon_id = $pathogen_taxon_id
				       AND uniprot_accession = '$pathogen_gene_uniprot_acc'
				     )
				   );

	    my $sql_result3 = $db_conn->prepare($sql_statement3);
	    $sql_result3->execute() or die $DBI::errstr;

	    # get the unique identifier for the inserted pathogen_gene_mutant record
	    my $sql_statement5 = qq(SELECT id FROM pathogen_gene_mutant
				    WHERE pathogen_gene_id = $pathogen_gene_id
				    AND ncbi_taxon_id = $pathogen_taxon_id
				    AND uniprot_accession = '$pathogen_gene_uniprot_acc');

	    my $sql_result5 = $db_conn->prepare($sql_statement5);
	    $sql_result5->execute() or die $DBI::errstr;
	    my @row5 = $sql_result5->fetchrow_array();
	    my $pathogen_gene_mutant_id = shift @row5;

	    # add the current pathogen gene mutant ID to the list for all genes
	    # (in case of multiple gene interaction)
	    push(@pathogen_gene_mutant_ids,$pathogen_gene_mutant_id);

         } # end foreach pathogen gene UniProt identifier


	 # get the largest available value for phi_base_accession,
	 # so that we know to increment this number
	 # (ACTUALLY, WE ONLY NEED TO DO THIS ONCE - AT VERY BEGINNING - NOT WITHIN THIS LOOP)
	 my $sql_statement1 = qq(SELECT phi_base_accession FROM interaction
				 WHERE id = 
				   (SELECT MAX(id) FROM interaction)
				);

	 my $sql_result1 = $db_conn->prepare($sql_statement1);
	 $sql_result1->execute() or die $DBI::errstr;
	 my @row1 = $sql_result1->fetchrow_array();
	 my $largest_phibase_accession = shift @row1;

	 # PHI-base accession number is the accession with 'PHI:I' prefix removed
	 my $interaction_num = substr($largest_phibase_accession,5);

	 # increment interaction counter, to become new PHI-base accession
	 $interaction_num++;
	 my $phi_base_accession = "PHI:I".$interaction_num;

	 # insert a new record into the interaction table, returning the interaction identifier
	 $sql_statement1 = qq(INSERT INTO interaction (phi_base_accession,curation_date) 
				   VALUES ('$phi_base_accession','$creation_date') RETURNING id;);
	 $sql_result1 = $db_conn->prepare($sql_statement1);
	 $sql_result1->execute() or die $DBI::errstr;
	 @row1 = $sql_result1->fetchrow_array();
	 $interaction_id = shift @row1;

	 # host taxon id is the taxon with 'NCBItaxon:' prefix removed
	 $host_taxon_id = substr($host_taxon,10);
	 # PubMed id is the pubmed string with 'PMID:' prefix removed
	 my $pubmed_id_num = substr($pubmed_id,5);

	 # add records for the literature, host, and tissue tables associated with the interaction,
	 # using the interaction id as a foreign key to the interaction table (where values exist)
         if (defined $host_taxon_id) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
					   VALUES ($interaction_id,$host_taxon_id);
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }
         if (defined $tissue) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_tissue (interaction_id,brenda_tissue_id) 
					   VALUES ($interaction_id,'$tissue');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }
         if (defined $disease) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_disease (interaction_id,disease_id) 
					   VALUES ($interaction_id,'$disease');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }
         if (defined $pubmed_id_num) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_literature (interaction_id, pubmed_id)
					   VALUES ($interaction_id, '$pubmed_id_num');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }


         if (defined $host_response_id) {

	    print "Host Response: $host_response_id\n";

	    # before we can insert the host response,
	    # we need to get the identifier for the appropriate host
	    my $sql_statement = qq(SELECT id FROM interaction_host
				     WHERE interaction_id = $interaction_id
				     AND ncbi_taxon_id = $host_taxon_id
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $interaction_host_id = shift @row;

	    # insert data into interaction_host_response table,
	    # with foreign keys to the interaction table and the host_response ontology
	    $host_response_count++;
	    $sql_statement = qq(INSERT INTO interaction_host_response (interaction_host_id, host_response_id)
				  VALUES ($interaction_id, '$host_response_id');
			       );
	    $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

         } # end if host response


	 # add records for each of the pathogen gene mutants associated with the interaction,
	 # using the interaction id as a foreign key to the interaction table
         # (in the case of a multiple gene interaction there will be multiple entries)
         foreach my $pathogen_gene_mutant_id (@pathogen_gene_mutant_ids) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_pathogen_gene_mutant (interaction_id,pathogen_gene_mutant_id) 
					   VALUES ($interaction_id,$pathogen_gene_mutant_id);
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }


         # Enter curator data into database
         # First, find out if the curator already exists, based on the email address
         # (the email will be a better identifier than the curator name)
         # NOTE: curation organisation is not currently captured by Canto
	 my $sql_statement = qq(SELECT id FROM curator
				  WHERE email = '$curator_email';
			       );

	 my $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $DBI::errstr;
	 my @row = $sql_result->fetchrow_array();
	 my $curator_id = shift @row;


         # if an existing curator is not found a new curator record needs to be added
         if (not defined $curator_id) {

	    # insert a new curator into the curator table, returning the curator identifier
	    $sql_statement = qq(INSERT INTO curator (name,email) 
				      VALUES ('$curator_name','$curator_email') RETURNING id;
				);
	    $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    $curator_id = shift @row;

         }

	 # insert data into interaction_curator table,
	 # with foreign keys to the interaction table and the curator table 
	 $sql_statement = qq(INSERT INTO interaction_curator (interaction_id, curator_id)
			       VALUES ($interaction_id, $curator_id);
			    );
	 $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $DBI::errstr;


         # Enter approver data into database
         # First, find out if the approver already exists as a curator
         # (an approver is just a special type of curator),
         # based on the email address 
         # (the email will be a better identifier than the approver name)
	 $sql_statement = qq(SELECT id FROM curator
				  WHERE email = '$approver_email';
			       );

	 $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $DBI::errstr;
	 @row = $sql_result->fetchrow_array();
	 $curator_id = shift @row;


         # if an existing curator is not found,
         # then a new curator record needs to be added for the approver
         if (not defined $curator_id) {

	    # insert a new curator into the curator table, returning the curator identifier
	    $sql_statement = qq(INSERT INTO curator (name,email) 
				      VALUES ('$approver_name','$approver_email') RETURNING id;
				);
	    $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    $curator_id = shift @row;

         }

	 # insert data into interaction_approver table,
	 # with foreign keys to the interaction table and the curator table 
	 $sql_statement = qq(INSERT INTO interaction_approver (interaction_id, curator_id)
			       VALUES ($interaction_id, $curator_id);
			    );
	 $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $DBI::errstr;

      } # end if first annotation

      # for each annotation, find the type of data being annotated
      # and the corresponding value and evidence code
      my $annot_type = $annotations[$annot_index]{'type'};
      my $annot_value = $annotations[$annot_index]{'term'};
      my $annot_evid_code = $annotations[$annot_index]{'evidence_code'};

      # make appropriate assignment depending on the type of data
#      if ($annot_type eq 'pathogen_taxon') {
#         $pathogen_taxon = $annot_value;
#         print "Pathogen Taxon: $pathogen_taxon\n";
#         # the pathogen taxon is split between the name and value, separated by colon
#         # need to extract the taxon value
#         my @pathogen_taxon_parts = split(/[:]/,$pathogen_taxon);
#         $pathogen_taxon_id = $pathogen_taxon_parts[1];
#         print "Pathogen Taxon ID: $pathogen_taxon_id\n";
#      } elsif ($annot_type eq 'phenotype_outcome') {
      if ($annot_type eq 'phenotype_outcome') {

         $phenotype_outcome = $annot_value;
         print "Phenotype Outcome: $phenotype_outcome\n";

	 # insert data into the appropriate interaction_phenotype_outcome table,
	 # with for a foreign key to the phenotype_outcome ontology
	 $phenotype_outcome_count++;
	 my $sql_statement = qq(INSERT INTO interaction_phenotype_outcome
			       (interaction_id, phenotype_outcome_id)
			       VALUES ($interaction_id, '$phenotype_outcome');
			    );
	 my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;


         # get the experiment specification term associated with the phenotype outcome,
         # which is given by the evidence code
         # get the term ID from the ontology, based on the name
         print "Evidence Code:$annot_evid_code\n";
         my $exp_spec_id = $exp_spec_ontology->get_term_by_name($annot_evid_code)->id;

         if (defined $exp_spec_id) {
            print "Experiment Specification ID:$exp_spec_id\n";
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_experiment_spec (interaction_id, experiment_spec_id) 
					   VALUES ($interaction_id,'$exp_spec_id');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }

=pod
      } elsif ($annot_type eq 'disease') {

         $disease = $annot_value;
         print "Disease: $disease\n";

	 # insert data into interaction_disease table,
	 # with foreign keys to the interaction table and the disease ontology
	 $disease_count++;
	 #$sql_statement = qq(INSERT INTO interaction_disease (interaction_id, disease_id, disease_severity_id)
	 my $sql_statement = qq(INSERT INTO interaction_disease (interaction_id, disease_id)
			       VALUES ($interaction_id, '$disease');
			    );
	 my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
=cut

      } elsif ($annot_type eq 'molecular_function' or $annot_type eq 'biological_process' or $annot_type eq 'cellular_component' ) {

         $go_id = $annot_value;
         $go_evid_code = $annot_evid_code;
         print "GO Annotation: $go_id ($go_evid_code)\n";

	 # insert data into interaction_go_annotation table,
	 # with foreign keys to the interaction table and the go_evidence_code_table 
	 my $sql_statement = qq(INSERT INTO interaction_go_annotation (interaction_id, go_id, go_evidence_code)
			       VALUES ($interaction_id, '$go_id', '$go_evid_code');
			    );
	 my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
	 $go_annotation_count++;

=pod
      } elsif ($annot_type eq 'host_response') {

         $host_response_id = $annot_value;
         print "Host Response: $host_response_id\n";

         # before we can insert the host response,
         # we need to get the identifier for the appropriate host
	 my $sql_statement = qq(SELECT id FROM interaction_host
				  WHERE interaction_id = $interaction_id
				  AND ncbi_taxon_id = $host_taxon_id
			       );

	 my $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $DBI::errstr;
	 my @row = $sql_result->fetchrow_array();
	 my $interaction_host_id = shift @row;

	 # insert data into interaction_host_response table,
	 # with foreign keys to the interaction table and the host_response ontology
	 $host_response_count++;
	 $sql_statement = qq(INSERT INTO interaction_host_response (interaction_host_id, host_response_id)
			       VALUES ($interaction_id, '$host_response_id');
			    );
	 $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
=cut

      } # end elsif annotation type

    } # end if annotation belongs to interaction

  } # end foreach annotation

} # end foreach interaction

print "\nTotal interactions:$interaction_count\n";
print "Total annotations :$annotation_count\n";
print "Total phenotype outcomes: $phenotype_outcome_count\n";
#print "Total diseases: $disease_count\n";
print "Total GO annotations: $go_annotation_count\n";
#print "Total host responses: $host_response_count\n";

#close (JSON_OUTPUT_FILE);

$db_conn->disconnect() or die "Failed to disconnect database\n";

