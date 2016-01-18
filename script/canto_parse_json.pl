#!/usr/bin/perl
use strict;
use warnings;

use JSON;
use DBI;
use JSON::Parse 'json_file_to_perl';
use Data::Compare;

use phibase_subroutines qw(connect_to_phibase); 

my $db_conn = connect_to_phibase();

# initialise counters to gather statistics
my $annotation_count = 0;
my $interaction_count = 0;
my $go_annotation_count = 0;
my $effector_annotation_count = 0;
my $host_response_count = 0;
my $phi_interaction_phenotype_count = 0;
my $phi_pathogen_phenotype_count = 0;
my $phi_host_phenotype_count = 0;
my $disease_count = 0;
my $go_mol_function_count = 0;
my $go_biol_process_count = 0;
my $go_cell_location_count = 0;
my $allele_count = 0;
# CURRENTLY UNUSED
#my $anti_infective_count = 0;
#my $inducer_count = 0;
#my $exp_spec_count = 0;

my $json_filename = '../input/canto/canto_2016-01-16.json';
#my $json_filename = '../input/canto/canto_2016-01-14.json';
#my $json_filename = '../input/canto/canto_2016-01-13.json';
#my $json_filename = '../input/canto/canto_with_annot_ext.json';
#my $json_filename = '../input/canto/canto_using_alleles.json';
#my $json_filename = '../input/canto/paper_with_complex_annot_extension.json';
#my $json_filename = '../input/canto/paper_with_double_triple_mutants.json';
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

# get the pathogen genes for the session
my %pathogen_genes = %{ $text_response->{'curation_sessions'}{$session_id}{'genes'} };

# get the pathogen alleles for the session
my %pathogen_alleles = %{ $text_response->{'curation_sessions'}{$session_id}{'alleles'} };

# get the pathogen genotypes for the session
my %pathogen_genotypes = %{ $text_response->{'curation_sessions'}{$session_id}{'genotypes'} };

# annotations are an array of hashes
my @annotations = @{ $text_response->{'curation_sessions'}{$session_id}{'annotations'} };

# create a hash of hashes to store the pathogen-host interactions
# where the key will be an identifier for each interaction
# and the value will be a hash of the values that identify each unique interaction
my %interaction_profiles;

# create a separate hash of hashes to store the GO annotations
# (which are not part of the pathogen-host interaction, but a direct property of the pathogen gene)
# where the key will be an identifier for each GO annotation
# and the value will be a hash of the values for each annotation
my @go_annotations;

# create a separate hash of hashes to store the effector gene annotations
# (which are not part of the pathogen-host interaction, but a direct property of the pathogen gene)
# where the key will be an identifier for each effector annotation
# and the value will be a hash of the values for each annotation
my @effector_annotations;

# create a separate hash of hashes to store the GO annotations
# (which are not part of the pathogen-host interaction, but a direct property of the pathogen gene)
# where the key will be an identifier for each GO annotation
# and the value will be a hash of the values for each annotation
my @post_trans_mod_annotations;

# Read in the experiment specification ontology
# so that the identifier can be matched from the
# available exp specification strring (from the evidence code)
print "Reading ontology files...\n";
my $obo_parser = OBO::Parser::OBOParser->new;
my $exp_spec_ontology = $obo_parser->work("../ontology/phibase/experiment_specification.obo");

# iterate through the array of annotations
foreach my $annot_index (0 .. $#annotations) {

   $annotation_count++;

   # create temporary hash to identify the interaction to which the phenotype annotation should belong
   # (GO annotations are associated directly with a pathogen gene, not an interaction)
   # this should include all the info that specifies the unqiue interaction,
   # such as the genotype identifier and annotation extensions
   my %annotation_interaction_hash = ();

   # get the current annotation hash
   my %annotation = %{ $text_response->{'curation_sessions'}{$session_id}{'annotations'}[$annot_index] };

   # GO annotation is annotated against a gene, whereas phenotypes are annotated against a genotype
   # so need to get correct annotation type for the identifier
   my $gene_or_genotype_id;
   my $annot_type = $annotation{'type'};
   if ($annot_type eq 'phi_phenotype') {
#   if ($annot_type eq 'phi_phenotype' or $annot_type eq 'phenotype_outcome') {
#   if ($annot_type eq 'phi_interaction_phenotype' or $annot_type eq 'phi_pathogen_phenotype' or $annot_type eq 'phi_host_phenotype') {
      $gene_or_genotype_id = $annotation{'genotype'};
   } elsif ($annot_type eq 'molecular_function' or $annot_type eq 'biological_process' or $annot_type eq 'cellular_component' or $annot_type eq 'effector' or $annot_type eq 'post_translational_modification') {
      $gene_or_genotype_id = $annotation{'gene'};
   } else {
      print STDERR "Error: Invalid annotation type: $annot_type\n"; 
   }

   # add gene or genotype to help identify the correct interaction
   $annotation_interaction_hash{'gene_or_genotype_id'} = $gene_or_genotype_id;

   # get the annotation extensions for the current annotation,
   # which are separated by commas
#   my $annot_extension_list = $annotation{'annotation_extension'};
#   my @annot_extensions = split(/,/,$annot_extension_list);

   # get the annotation extensions for the current annotation,
   # which is an array of hashes (one hash for each extension)
   my @annot_extensions = @{ $annotation{'extension'} };


   my $host_taxon_id;
   my $host_strain_name;
   my $tissue;
   my $disease;
   my $inducer_chebi_id;
   my $inducer_cas_num;
   my $inducer_gene;
   my $anti_infective;
   my $host_interacting_protein;
   my $pathogen_interacting_protein;
   my $pathogen_strain_id;
   my $pathogen_strain_name;

   # identify each annotation extension in the list
#   foreach my $annot_ext (@annot_extensions) {
   foreach my $annot_ext_index (0 .. $#annot_extensions) {

     # get the hash for the annotation extension
     my %annot_ext_hash = %{ $annot_extensions[$annot_ext_index] };
     my $annot_ext_relation = $annot_ext_hash{'relation'}; 
     my $annot_ext_value = $annot_ext_hash{'rangeValue'}; 

     # if extension relation is 'experimental_host_id', then assign range value to host taxon id
     if ($annot_ext_relation =~ /^experimental_host_id/) {
        $host_taxon_id = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_taxon_id'} = $host_taxon_id;
     }

     # if extension relation is 'experimental_host_strain', then assign range value to host strain name
     if ($annot_ext_relation =~ /^experimental_host_strain/) {
        $host_strain_name = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_strain_name'} = $host_strain_name;
     }

     # if extension relation is 'infects_tissue', then assign range value to tissue
     if ($annot_ext_relation =~ /^infects_tissue/) {
        $tissue = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_tissue'} = $tissue;
     }

     # if extension relation is 'causes_disease', then assign range value to disease
     if ($annot_ext_relation =~ /^causes_disease/) {
        $disease = $annot_ext_value;
	# add annotation extension to help identify the correct interaction
	$annotation_interaction_hash{'disease'} = $disease;
     }

     # if extension relation is 'induced_by_chemical', then assign range value to inducer_chemical
     if ($annot_ext_relation =~ /^induced_by_chemical/) {

        # Need to find out if the chemical input is a ChEBI ID or a CAS Registry Number.
        # If the value begins with 'CHEBI', then it is a ChEBI ID,
        # otherwise, assume that it is a CAS Registry Number.
	# Add an appropriate annotation extension, depending on the identifier
        if ($annot_ext_value =~ /^CHEBI/) {
	   $inducer_chebi_id = $annot_ext_value;
	   $annotation_interaction_hash{'inducer_chemical_chebi_id'} = $inducer_chebi_id;
        } else {
	   $inducer_cas_num = $annot_ext_value;
	   $annotation_interaction_hash{'inducer_chemical_cas_num'} = $inducer_cas_num;
        }

     }

     # if extension relation is 'induced_by_gene', then assign range value to inducer_gene
     if ($annot_ext_relation =~ /^induced_by_gene/) {
        $inducer_gene = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'inducer_gene'} = $inducer_gene;
     }

     # if extension relation is 'anti_infective', then assign range value to anti_infective
     if ($annot_ext_relation =~ /^anti_infective/) {
        $anti_infective = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'anti_infective'} = $anti_infective;
     }

     # if extension relation is 'host_interacting_protein', then assign range value to host_interacting_protein
     if ($annot_ext_relation =~ /^host_interacting_protein/) {
        $host_interacting_protein = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_interacting_protein'} = $host_interacting_protein;
     }

     # if extension relation is 'pathogen_interacting_protein', then assign range value to pathogen_interacting_protein
     if ($annot_ext_relation =~ /^pathogen_interacting_protein/) {
        $pathogen_interacting_protein = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'pathogen_interacting_protein'} = $pathogen_interacting_protein;
     }

     # if extension relation is 'pathogen_strain_id', then assign range value to pathogen_strain_id
     if ($annot_ext_relation =~ /^pathogen_strain_id/) {
        $pathogen_strain_id = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'pathogen_strain_id'} = $pathogen_strain_id;
     }

     # if extension relation is 'pathogen_strain_name', then assign range value to pathogen_strain_name
     if ($annot_ext_relation =~ /^pathogen_strain_name/) {
        $pathogen_strain_name = $annot_ext_value;
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'pathogen_strain_name'} = $pathogen_strain_name;
     }

   } # end foreach annotation extension


   # declare flag to indicate if the interaction associated with
   # the current annotation has been found
   my $interaction_found = 0;

   # for a phenotype annotation,
   # iterate through each of the existing interactions to find out
   # if the current annotation profile matches an existing interaction
   # if so, then add the corresponding interaction id to the annotation
   # otherwise, create a new interaction profile for the annotation
   #
   # GO annotations are not linked to specific pathogen-host interactions,
   # so for a GO annotation, simply add the annotation to the array of GO annotations
   #
   if ($annot_type eq 'phi_phenotype') {
#if ($annot_type eq 'phi_phenotype' or $annot_type eq 'phenotype_outcome') {
#if ($annot_type eq 'phi_interaction_phenotype' or $annot_type eq 'phi_pathogen_phenotype' or $annot_type eq 'phi_host_phenotype') {

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

      } # end foreach interaction id

      # if a matching interaction profile was not found,
      # then a new interaction profile should be created based on the annotation profile,
      # with the new interaction ID assigned to the annotation
      if (not $interaction_found) {
	 $interaction_count++;
	 $interaction_profiles{$interaction_count} =  { %annotation_interaction_hash } ;
	 $annotations[$annot_index]{"interaction_id"} = $interaction_count;
      }

   } elsif ($annot_type eq 'molecular_function' or $annot_type eq 'biological_process' or $annot_type eq 'cellular_component') {

      # add the current annotation to the array of GO annotations;
      push (@go_annotations, $annotations[$annot_index]);

   } elsif ($annot_type eq 'effector') {

      # add the current annotation to the array of effector annotations;
      push (@effector_annotations, $annotations[$annot_index]);

   } elsif ($annot_type eq 'post_translational_modification') {

      # add the current annotation to the array of post translational modification annotations;
      push (@post_trans_mod_annotations, $annotations[$annot_index]);

   } else {
      print STDERR "Error: Invalid annotation type: $annot_type\n"; 
   }


} # end foreach annotation


# iterate through each interaction and retrieve all of the annotations
# associated with the interaction
foreach my $int_id (1 .. $interaction_count) {

  print "\nINTERACTION:$int_id\n";

  # flag indicating first annotation in the interaction
  my $first_annot_of_int_flag = 1;

  # declare variables
  my $pubmed_id;
  #my $host_taxon;
  my $host_taxon_id;
  my $host_strain_name;
  my $tissue;
  my $curator_name;
  my $curator_email;
  my $pathogen_species;
  my $pathogen_taxon;
  my $phi_interaction_phenotype;
  my $phi_pathogen_phenotype;
  my $phi_host_phenotype;
  my $disease;
  my $go_id;
  my $go_evid_code;
  my $interaction_id;
  my $inducer_chebi_id;
  my $inducer_cas_num;
  my $inducer_gene;
  my $anti_infective;
  my $host_interacting_protein;
  my $pathogen_interacting_protein;
  my $pathogen_strain_id;
  my $pathogen_strain_name;
  my $creation_date;

  my $gene_or_genotype_id;

  # declare array to hold all pathogen genes
  # involved in the interaction
  my @pathogen_genes;

  # iterate through all annotations, to find out which ones
  # are associated with the current interaction
  foreach my $annot_index (0 .. $#annotations) {

    # get the interaction ID assocaiated with the annotation
    my $annot_int_id = $annotations[$annot_index]{"interaction_id"};

    # check if this interaction ID matches the current interaction
    # note that the interaction id will only be defined if the current interaction
    # is a phenotype annotation, but not if it is a GO annotation
    if (defined $annot_int_id and $annot_int_id == $int_id) {

      # get the hash of annotation details
      my %annotation = %{ $annotations[$annot_index] };

      # since most of the information about an interaction is repeated
      # in all of the annotations, we can just retrieve the relevant data
      # from the first annotation
      if ($first_annot_of_int_flag) {

         $first_annot_of_int_flag = 0;

         $pubmed_id = $annotation{'publication'};
         print "PubMed Article: $pubmed_id\n";

         my $annot_type = $annotation{'type'};
         $gene_or_genotype_id = $annotation{'genotype'};

         # get the alleles for this interaction, as an array in the specified genotype
         my @genotype_alleles = @{ $pathogen_genotypes{$gene_or_genotype_id}{"alleles"} };

         # iterate through the array of alleles in the genotype
         # to identify the list of pathogen genes involved in the
         # current interaction
         foreach my $genotype_allele_index (0 .. $#genotype_alleles) {

           # get the allele hash from the genotype      
	   my %allele = %{ $genotype_alleles[$genotype_allele_index] };

           # retrieve the identifier for the allele
	   my $allele_id = $allele{"id"};

           # get the identifier for the gene to which the allele belongs
	   my $allele_gene = $pathogen_alleles{$allele_id}{"gene"};

           # use the gene identifer to get the UniProt accession of the gene
	   my $path_uniprot_acc = $pathogen_genes{$allele_gene}{"uniquename"};

	   # add the uniprot accession for this gene to the array
	   # of genes, which includes all pathogen genes involved 
	   # in the current interaction
	   push(@pathogen_genes,$path_uniprot_acc);

	 } # end foreach genotype allele

         # get the annotation extensions for the current annotation,
         # which are separated by commas
#         my $annot_extension_list = $annotation{'annotation_extension'};
#         my @annot_extensions = split(/,/,$annot_extension_list);

   # get the annotation extensions for the current annotation,
   # which is an array of hashes (one hash for each extension)
   my @annot_extensions = @{ $annotation{'extension'} };

         # identify each annotation extension in the list
#         foreach my $annot_ext (@annot_extensions) {
          foreach my $annot_ext_index (0 .. $#annot_extensions) {

	    # get the hash for the annotation extension
	    my %annot_ext_hash = %{ $annot_extensions[$annot_ext_index] };
	    my $annot_ext_relation = $annot_ext_hash{'relation'}; 
	    my $annot_ext_value = $annot_ext_hash{'rangeValue'}; 

	    # if extension relation is 'experimental_host_id', then assign range value to host taxon id
	    if ($annot_ext_relation =~ /^experimental_host_id/) {
	      $host_taxon_id = $annot_ext_value;
	    }

	    # if extension relation is 'experimental_host_strain', then assign range value to host strain name
	    if ($annot_ext_relation =~ /^experimental_host_strain/) {
	      $host_strain_name = $annot_ext_value;
	    }

	    # if extension relation is 'infects_tissue', then assign range value to tissue
            if ($annot_ext_relation =~ /^infects_tissue/) {
	      $tissue = $annot_ext_value;
            }

	    # if extension relation is 'causes_disease', then assign range value to disease
            if ($annot_ext_relation =~ /^causes_disease/) {
	      $disease = $annot_ext_value;
            }

#	    # if extension relation is 'induced_by_chemical', then assign range value to inducer
#	    if ($annot_ext_relation =~ /^induced_by_chemical/) {
#	      $inducer_chebi_id = $annot_ext_value;
#	    }

	    # if extension relation is 'induced_by_chemical', then assign range value to inducer_chemical
	    if ($annot_ext_relation =~ /^induced_by_chemical/) {

	       # Need to find out if the chemical input is a ChEBI ID or a CAS Registry Number.
	       # If the value begins with 'CHEBI', then it is a ChEBI ID,
	       # otherwise, assume that it is a CAS Registry Number.
	       if ($annot_ext_value =~ /^CHEBI/) {
		  $inducer_chebi_id = $annot_ext_value;
	       } else {
		  $inducer_cas_num = $annot_ext_value;
	       }

	    }

	    # if extension relation is 'induced_by_gene', then assign range value to inducer
	    if ($annot_ext_relation =~ /^induced_by_gene/) {
	      $inducer_gene = $annot_ext_value;
	    }

	    # if extension relation is 'anti_infective', then assign range value to anti_infective
	    if ($annot_ext_relation =~ /^anti_infective/) {
	      $anti_infective = $annot_ext_value;
	    }

	    # if extension relation is 'host_interacting_protein', then assign range value to host_interacting_protein
	    if ($annot_ext_relation =~ /^host_interacting_protein/) {
	       $host_interacting_protein = $annot_ext_value;
	    }

	    # if extension relation is 'pathogen_interacting_protein', then assign range value to pathogen_interacting_protein
	    if ($annot_ext_relation =~ /^pathogen_interacting_protein/) {
	       $pathogen_interacting_protein = $annot_ext_value;
	    }

	    # if extension relation is 'pathogen_strain_id', then assign range value to pathogen_strain_id
	    if ($annot_ext_relation =~ /^pathogen_strain_id/) {
	       $pathogen_strain_id = $annot_ext_value;
	    }

	    # if extension relation is 'pathogen_strain_name', then assign range value to pathogen_strain_name
	    if ($annot_ext_relation =~ /^pathogen_strain_name/) {
	       $pathogen_strain_name = $annot_ext_value;
	    }

         } # end foreach annotation extension


         print "Host Taxon: $host_taxon_id\n";
         print "Host Strain Name: $host_strain_name\n" if defined $host_strain_name;
         print "Tissue: $tissue\n" if defined $tissue;
         print "Disease: $disease\n" if defined $disease;
         print "Inducer Chemical: $inducer_chebi_id\n" if defined $inducer_chebi_id;
         print "Inducer Chemical: $inducer_cas_num\n" if defined $inducer_cas_num;
         print "Inducer Gene: $inducer_gene\n" if defined $inducer_gene;
         print "Anti-infective: $anti_infective\n" if defined $anti_infective;
         print "Host Interaction protein: $host_interacting_protein\n" if defined $host_interacting_protein;
         print "Pathogen Interaction protein: $pathogen_interacting_protein\n" if defined $pathogen_interacting_protein;
         print "Pathogen Strain ID: $pathogen_strain_id\n" if defined $pathogen_strain_id;
         print "Pathogen Strain Name: $pathogen_strain_name\n" if defined $pathogen_strain_name;

         $curator_name = $annotation{'curator'}{'name'};
         $curator_email = $annotation{'curator'}{'email'};
         print "Curator: $curator_name, $curator_email\n";
         print "Approver: $approver_name, $approver_email\n";

         $creation_date = $annotation{'creation_date'};
         print "Creation Date: $creation_date\n";


         # if a pathogen strain taxon ID has been explicitly given in the annotation extension
         # then this taxon ID should be used for all interactions, otherwise use the taxon ID
         # retrieved from the UniProt accessions supplied
         if (defined $pathogen_strain_id) {
#	    # pathogen taxon id becomes the strain taxon with 'NCBItaxon:' prefix removed
#            $pathogen_taxon_id = substr($pathogen_strain_id,10);
            $pathogen_taxon_id = $pathogen_strain_id;
         }


	 # insert data into the pathogen_gene table, for each pathogen gene
	 # if it does not exist already
         # (based on combination of taxon id, UniProt accession, and where available the pathogen strain name)
	 foreach my $pathogen_gene_uniprot_acc (@pathogen_genes) {

            print "Pathogen Gene: $pathogen_gene_uniprot_acc\n";

            my $sql_statement2;
            if (defined $pathogen_strain_name) {
	       $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, pathogen_strain_name, uniprot_accession) 
				    SELECT $pathogen_taxon_id,'$pathogen_strain_name','$pathogen_gene_uniprot_acc'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM pathogen_gene
				      WHERE ncbi_taxon_id = $pathogen_taxon_id
				      AND pathogen_strain_name = '$pathogen_strain_name'
				      AND uniprot_accession = '$pathogen_gene_uniprot_acc'
				   ));
            } else {
	       $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
				    SELECT $pathogen_taxon_id,'$pathogen_gene_uniprot_acc'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM pathogen_gene
				      WHERE ncbi_taxon_id = $pathogen_taxon_id
				      AND uniprot_accession = '$pathogen_gene_uniprot_acc'
				   ));
            }

	    my $sql_result2 = $db_conn->prepare($sql_statement2);
	    $sql_result2->execute() or die $DBI::errstr;

	    # get the unique identifier for the inserted pathogen_gene record
	    my $sql_statement4 = qq(SELECT id FROM pathogen_gene
				    WHERE ncbi_taxon_id = $pathogen_taxon_id
				    AND uniprot_accession = '$pathogen_gene_uniprot_acc');

	    my $sql_result4 = $db_conn->prepare($sql_statement4);
	    $sql_result4->execute() or die $DBI::errstr;
	    my @row4 = $sql_result4->fetchrow_array();
	    my $pathogen_gene_id = shift @row4;

            # iterate through the full list of alleles and identify those
            # that belong to the current pathogen gene, based on comparison
            # of the UniProt accessions
            # when an allele belonging to the current pathogen gene is found,
            # retrieve the relevant allele details, then insert a record into
            # the pathogen_gene_allele table
            foreach my $allele_id (keys %pathogen_alleles) {

               # get the hash of allele details
               my %pathogen_allele = %{ $pathogen_alleles{$allele_id} };

               # get the gene identifier that the allele belongs to
               my $pathogen_allele_gene_id = $pathogen_allele{"gene"};

               # get the hash of gene details from the pathogen_genes hash, using the identifier
               my %pathogen_allele_gene = %{ $pathogen_genes{$pathogen_allele_gene_id} };

               # get the UniProt accession of the gene from the gene details hash
               my $pathogen_allele_gene_uniprot_acc = $pathogen_allele_gene{"uniquename"};

               # compare the UniProt accessions of the allele's gene and the current pathogen gene
               if ($pathogen_allele_gene_uniprot_acc eq $pathogen_gene_uniprot_acc) {

                   # now we know that this allele needs to be added as
                   # a pathogen_gene_allele record for the current pathogen gene
		   $allele_count++;

                   # retreive other allele details from the allele hash
                   my $gene_allele_name = $pathogen_allele{"name"};
                   my $gene_allele_type = $pathogen_allele{"allele_type"};
                   my $gene_allele_desc = $pathogen_allele{"description"};

		   # insert data into pathogen_gene_allele table, including foreign key to pathogen_gene table
                   # where it does not already exist
		   my $sql_statement3 = qq(INSERT INTO pathogen_gene_allele (pathogen_gene_id, allele_name, allele_type, allele_description) 
					    SELECT $pathogen_gene_id, '$gene_allele_name', '$gene_allele_type', '$gene_allele_desc'
					    WHERE NOT EXISTS (
					      SELECT 1 FROM pathogen_gene_allele
					      WHERE pathogen_gene_id = $pathogen_gene_id
					      AND allele_name = '$gene_allele_name'
					      AND allele_type = '$gene_allele_type'
					    )
					  );

		   my $sql_result3 = $db_conn->prepare($sql_statement3);
		   $sql_result3->execute() or die $DBI::errstr;

               } # end if UniProt accessions match

            } # end foreach pathogen allele

         } # end foreach pathogen gene UniProt identifier


	 # get the largest available value for phi_base_accession,
	 # so that we know to increment this number
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

	 # PubMed id is the pubmed string with 'PMID:' prefix removed
	 my $pubmed_id_num = substr($pubmed_id,5);

	 # add records for the all of the tables associated with the interaction,
	 # using the interaction id as a foreign key to the interaction table

         if (defined $host_taxon_id) {

            # if the host strain and/or host interacting protein are defined, then these are also added to the interaction_host table,
            # otherwise, they are left out of the SQL insert statement
            my $inner_sql_statement;
            if (defined $host_strain_name and defined $host_interacting_protein) { # both are defined
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id,host_strain_name,first_target_uniprot_accession) 
					   VALUES ($interaction_id,$host_taxon_id,'$host_strain_name','$host_interacting_protein');
					);
            }
            elsif (defined $host_strain_name) { # only host strain defined
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id,host_strain_name) 
					   VALUES ($interaction_id,$host_taxon_id,'$host_strain_name');
					);
            }
            elsif (defined $host_interacting_protein) { # only interacting protein defined
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id,first_target_uniprot_accession) 
					   VALUES ($interaction_id,$host_taxon_id,'$host_interacting_protein');
					);
            } else { # neither defined
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
					   VALUES ($interaction_id,$host_taxon_id);
					);
            }
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

         }

         if (defined $disease) {
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_disease (interaction_id,disease_id) 
					   VALUES ($interaction_id,'$disease');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }

	 if (defined $pathogen_interacting_protein) {
	    my $inner_sql_statement = qq(
	                                 INSERT INTO pathogen_interacting_protein (interaction_id, uniprot_accession)
				           VALUES ($interaction_id, '$pathogen_interacting_protein');
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


         if (defined $tissue) {

	    # before we can insert the host tissue table,
	    # we need to get the identifier for the appropriate host
	    my $sql_statement = qq(SELECT id FROM interaction_host
				   WHERE interaction_id = $interaction_id
				   AND ncbi_taxon_id = $host_taxon_id
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $dbi::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $interaction_host_id = shift @row;

	    # insert data into interaction_host_tissue table,
	    # with foreign keys to the interaction_host table and the BRENDA host tissue ontology
	    my $inner_sql_statement = qq(
					 INSERT INTO interaction_host_tissue (interaction_host_id,brenda_tissue_id) 
					   VALUES ($interaction_host_id,'$tissue');
					);
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
         }


         if (defined $inducer_chebi_id) {

	    # insert data into the chemical_table,
	    # if it does not exist already (based on the ChEBI ID)
	    my $sql_statement = qq(INSERT INTO chemical (chebi_id) 
				      SELECT '$inducer_chebi_id'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM chemical
				      WHERE chebi_id = '$inducer_chebi_id'
				    )
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;

	    # before we can insert the inducer chemical,
	    # we need to get the identifier for the chemical
	    $sql_statement = qq(SELECT id FROM chemical
				     WHERE chebi_id = '$inducer_chebi_id'
				  );
	    $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $chemical_id = shift @row;

            # insert data into interaction_inducer_chemical table,
            # with foreign keys to the interaction table and the chemical table 
	    $sql_statement = qq(INSERT INTO interaction_inducer_chemical (interaction_id, chemical_id)
                                  VALUES ($interaction_id, $chemical_id);
                               );
	    $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

         } # end if inducer chemical ChEBI ID


         if (defined $inducer_cas_num) {

	    # insert data into the chemical_table,
	    # if it does not exist already (based on the CAS number)
	    my $sql_statement = qq(INSERT INTO chemical (cas_registry) 
				      SELECT '$inducer_cas_num'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM chemical
				      WHERE cas_registry = '$inducer_cas_num'
				    )
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;

	    # before we can insert the inducer chemical,
	    # we need to get the identifier for the chemical
	    $sql_statement = qq(SELECT id FROM chemical
				     WHERE cas_registry = '$inducer_cas_num'
				  );
	    $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $chemical_id = shift @row;

            # insert data into interaction_inducer_chemical table,
            # with foreign keys to the interaction table and the chemical table 
	    $sql_statement = qq(INSERT INTO interaction_inducer_chemical (interaction_id, chemical_id)
                                  VALUES ($interaction_id, $chemical_id);
                               );
	    $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

         } # end if inducer chemical CAS number


         if (defined $inducer_gene) {

            # MAY WANT TO CHECK IF THE UNIPROT ACCESSION IS VALID UNIPROT FORMAT

            # insert data into interaction_inducer_gene table,
            # with foreign keys to the interaction table and UniProt 
	    my $sql_statement = qq(INSERT INTO interaction_inducer_gene (interaction_id, uniprot_accession)
                                     VALUES ($interaction_id, '$inducer_gene');
                                  );
	    my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

         } # end if inducer gene


         if (defined $anti_infective) {

	    # insert data into the chemical_table,
	    # if it does not exist already (based on the ChEBI ID)
	    my $sql_statement = qq(INSERT INTO chemical (chebi_id) 
				      SELECT '$anti_infective'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM chemical
				      WHERE chebi_id = '$anti_infective'
				    )
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;

	    # before we can insert the anti-infective chemical,
	    # we need to get the identifier for the chemical
	    $sql_statement = qq(SELECT id FROM chemical
				     WHERE chebi_id = '$anti_infective'
				  );
	    $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $chemical_id = shift @row;

            # insert data into interaction_anti_infective_chemical table,
            # with foreign keys to the interaction table and the chemical table 
	    $sql_statement = qq(INSERT INTO interaction_anti_infective_chemical (interaction_id, chemical_id)
                                  VALUES ($interaction_id, $chemical_id);
                               );
	    $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

         } # end if anti-infective


         # iterate through the array of alleles in the current genotype
         # inserting a record in the interaction_pathogen_gene_allele table
         foreach my $genotype_allele_index (0 .. $#genotype_alleles) {

            # get the allele hash available from the genotype 
	    my %allele = %{ $genotype_alleles[$genotype_allele_index] };

            # get the identifier for the allele 
	    my $allele_id = $allele{"id"};

            # get the allele expression used for the current interaction
	    my $allele_expression = $allele{"expression"};

            # get the allele name and type
	    my $allele_name = $pathogen_alleles{$allele_id}{"name"};
	    my $allele_type = $pathogen_alleles{$allele_id}{"allele_type"};

	    # get the unique identifier for the relevant pathogen_gene_allele record
            # based on allele name and allele type
	    my $sql_statement5 = qq(SELECT id FROM pathogen_gene_allele
				    WHERE allele_name = '$allele_name'
				    AND allele_type = '$allele_type');

	    my $sql_result5 = $db_conn->prepare($sql_statement5);
	    $sql_result5->execute() or die $DBI::errstr;
	    my @row5 = $sql_result5->fetchrow_array();
	    my $pathogen_gene_allele_id = shift @row5;

	    # add a record for each allele associated with the interaction,
	    # using the interaction_id and pathogen_gene_allele_id as a foreign keys,
	    # (in the case of a multiple allele interaction there will be multiple entries)

            # USE DIFFERENT INSERT STATEMENTS, DEPENDING IF THE EXPRESSION WAS GIVEN
	    my $inner_sql_statement;
            if (defined $allele_expression) {
	       $inner_sql_statement = qq(
		  	                 INSERT INTO interaction_pathogen_gene_allele 
                                           (interaction_id, pathogen_gene_allele_id, allele_expression) 
					   VALUES ($interaction_id, $pathogen_gene_allele_id, '$allele_expression');
					);
            } else {
	       $inner_sql_statement = qq(
		  	                 INSERT INTO interaction_pathogen_gene_allele 
                                           (interaction_id, pathogen_gene_allele_id) 
					   VALUES ($interaction_id, $pathogen_gene_allele_id);
					);
            }
	    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

	 } # end foreach genotype allele


         # Enter curator data into database
         # First, find out if the curator already exists, based on the email address
         # (the email will be a better identifier than the curator name)
         # note: curation organisation is not currently captured by Canto
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

      # insert the appropriate record for the annotation, depending on the type of data
      if ($annot_type eq 'phi_phenotype') {
#if ($annot_type eq 'phi_interaction_phenotype') {

         $phi_interaction_phenotype = $annot_value;
         #print "PHI Interaction Phenotype: $phi_interaction_phenotype\n";

	 # insert data into the appropriate interaction_phi_interaction_phenotype table,
	 # with for a foreign key to the phi_phenotype ontology
	 $phi_interaction_phenotype_count++;
	 my $sql_statement = qq(INSERT INTO interaction_phi_interaction_phenotype
			       (interaction_id, phi_phenotype_id, phi_evidence)
			       VALUES ($interaction_id, '$phi_interaction_phenotype', '$annot_evid_code');
			    );
	 my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

      } elsif ($annot_type eq 'phi_pathogen_phenotype') {

         $phi_pathogen_phenotype = $annot_value;
         #print "PHI Pathogen Phenotype: $phi_pathogen_phenotype\n";

	 # insert data into the appropriate interaction_phi_interaction_pathogen table,
	 # with for a foreign key to the phi_phenotype ontology
	 $phi_pathogen_phenotype_count++;
	 my $sql_statement = qq(INSERT INTO interaction_phi_pathogen_phenotype
			       (interaction_id, phi_phenotype_id, phi_evidence)
			       VALUES ($interaction_id, '$phi_pathogen_phenotype', '$annot_evid_code');
			    );
	 my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

      } elsif ($annot_type eq 'phi_host_phenotype') {

         $phi_host_phenotype = $annot_value;
         #print "PHI Host Phenotype: $phi_host_phenotype\n";

	 # before we can insert the host phenotype,
	 # we need to get the identifier for the appropriate host
	 my $sql_statement = qq(SELECT id FROM interaction_host
				WHERE interaction_id = $interaction_id
				AND ncbi_taxon_id = $host_taxon_id
			       );

	 my $sql_result = $db_conn->prepare($sql_statement);
	 $sql_result->execute() or die $dbi::errstr;
	 my @row = $sql_result->fetchrow_array();
	 my $interaction_host_id = shift @row;

	 # insert data into the appropriate interaction_phi_host_phenotype table,
	 # with for a foreign key to the phi_phenotype ontology
	 $phi_host_phenotype_count++;
	 my $sql_statement2 = qq(INSERT INTO interaction_phi_host_phenotype
			           (interaction_host_id, phi_phenotype_id, phi_evidence)
			           VALUES ($interaction_host_id, '$phi_host_phenotype', '$annot_evid_code');
			        );
	 my $sql_result2 = $db_conn->do($sql_statement2) or die $DBI::errstr;

      } # end elsif annotation type

    } # end if annotation belongs to interaction

  } # end foreach annotation

} # end foreach interaction


# iterate through all GO annotations
# inserting each into the database, associated with the
# relevant pathogen gene and PubMed ID
foreach my $go_annot_ref (@go_annotations) {

  # declare variables
  my $uniprot_acc;
  my $pubmed_id;
  my $pathogen_species;
  my $pathogen_taxon;
  my $go_id;
  my $go_evid_code;
  my $pathogen_strain_id;
  my $pathogen_strain_name;

  # increment overall GO annotation count
  $go_annotation_count++;
 
  # get the hash of GO annotation details
  my %go_annotation = %{ $go_annot_ref };

  # get the type of GO annotation
  my $annot_type = $go_annotation{'type'};

  # increment the relevant counter,
  # depending on the type of GO annotation
  if ($annot_type eq 'molecular_function') {
    $go_mol_function_count++;
  } elsif ($annot_type eq 'biological_process') {
    $go_biol_process_count++;
  } elsif ($annot_type eq 'cellular_component' ) {
    $go_cell_location_count++;
  }

  # get the gene identifier
  my $gene_id = $go_annotation{'gene'};

  # get the gene from the hash of pathogen genes, using the ID
  my %gene = %{ $pathogen_genes{$gene_id} };
           
  # get the UniProt accession for the gene
  my $pathogen_gene_uniprot_acc = $gene{"uniquename"};

  # get the PubMed ID associated with the GO annotation
  $pubmed_id = $go_annotation{'publication'};

  # get the annotation extensions for the GO annotation
  # which is an array of hashes (one hash for each extension)
  my @annot_extensions = @{ $go_annotation{'extension'} };

  # from the annotation extensions in the list
  # find out if a specific pathogen strain taxon ID has been given
  # and/or a specific pathogen strain name
  foreach my $annot_ext_index (0 .. $#annot_extensions) {

    # get the hash for the annotation extension
    my %annot_ext_hash = %{ $annot_extensions[$annot_ext_index] };
    my $annot_ext_relation = $annot_ext_hash{'relation'}; 
    my $annot_ext_value = $annot_ext_hash{'rangeValue'}; 

    # if extension relation is 'pathogen_strain_id', then assign range value to pathogen_strain_id
    if ($annot_ext_relation =~ /^pathogen_strain_id/) {
       $pathogen_strain_id = $annot_ext_value;
    }

    # if extension relation is 'pathogen_strain_name', then assign range value to pathogen_strain_name
    if ($annot_ext_relation =~ /^pathogen_strain_name/) {
       $pathogen_strain_name = $annot_ext_value;
    }

  } # end foreach annotation extension

  # if a pathogen strain taxon ID has been explicitly given in the annotation extension
  # then this taxon ID should be used for all interactions, otherwise use the taxon ID
  # retrieved from the UniProt accessions supplied
  if (defined $pathogen_strain_id) {
    $pathogen_taxon_id = $pathogen_strain_id;
  }

  # insert data into the pathogen_gene table,
  # if it does not exist already
  # (based on combination of taxon id, UniProt accession, and where available the pathogen strain name)
  my $sql_statement2;
  if (defined $pathogen_strain_name) {
     $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
			  SELECT $pathogen_taxon_id,'$pathogen_strain_name','$pathogen_gene_uniprot_acc'
			  WHERE NOT EXISTS (
			    SELECT 1 FROM pathogen_gene
			    WHERE ncbi_taxon_id = $pathogen_taxon_id
			    AND pathogen_strain_name = '$pathogen_strain_name'
			    AND uniprot_accession = '$pathogen_gene_uniprot_acc'
			 ));
  } else {
     $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
			  SELECT $pathogen_taxon_id,'$pathogen_gene_uniprot_acc'
			  WHERE NOT EXISTS (
			    SELECT 1 FROM pathogen_gene
			    WHERE ncbi_taxon_id = $pathogen_taxon_id
			    AND uniprot_accession = '$pathogen_gene_uniprot_acc'
			 ));
  }

  my $sql_result2 = $db_conn->prepare($sql_statement2);
  $sql_result2->execute() or die $DBI::errstr;

  # get the unique identifier for the inserted pathogen_gene record
  my $sql_statement4 = qq(SELECT id FROM pathogen_gene
		          WHERE ncbi_taxon_id = $pathogen_taxon_id
		          AND uniprot_accession = '$pathogen_gene_uniprot_acc');

  my $sql_result4 = $db_conn->prepare($sql_statement4);
  $sql_result4->execute() or die $DBI::errstr;
  my @row4 = $sql_result4->fetchrow_array();
  my $pathogen_gene_id = shift @row4;

  $go_id = $go_annotation{'term'};
  $go_evid_code = $go_annotation{'evidence_code'};

  # insert data into pathogen_gene_go_annotation table,
  # with foreign keys to the pathogen_gene table, the go_evidence_code_table, the GO ontology, and PubMed,
  # unless it already exists 
  my $sql_statement = qq(INSERT INTO pathogen_gene_go_annotation (pathogen_gene_id, pubmed_id, go_id, go_evidence_code)
			 SELECT $pathogen_gene_id, '$pubmed_id', '$go_id', '$go_evid_code'
			 WHERE NOT EXISTS (
			   SELECT 1 FROM pathogen_gene_go_annotation
			   WHERE pathogen_gene_id = $pathogen_gene_id
                           AND pubmed_id = '$pubmed_id'
                           AND go_id = '$go_id'
			));
  my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

  # get the unique identifier for the inserted pathogen_gene_go_annotation record
  $sql_statement4 = qq(SELECT id FROM pathogen_gene_go_annotation
			  WHERE pathogen_gene_id = $pathogen_gene_id
                          AND pubmed_id = '$pubmed_id'
                          AND go_id = '$go_id'
		         );

  $sql_result4 = $db_conn->prepare($sql_statement4);
  $sql_result4->execute() or die $DBI::errstr;
  @row4 = $sql_result4->fetchrow_array();
  my $pathogen_gene_go_annotation_id = shift @row4;

  # for each annotation extension in the list
  # add a pathogen_gene_go_annot_ext record,
  # with the relevant foreign key to the pathogen_gene_go_annotation table
  foreach my $annot_ext_index (0 .. $#annot_extensions) {

    # get the hash for the annotation extension
    my %annot_ext_hash = %{ $annot_extensions[$annot_ext_index] };
    my $annot_ext_relation = $annot_ext_hash{'relation'}; 
    my $annot_ext_value = $annot_ext_hash{'rangeValue'}; 

    my $sql_statement = qq(INSERT INTO pathogen_gene_go_annot_ext (pathogen_gene_go_annotation_id, go_annot_ext_relation, go_annot_ext_value)
			 VALUES ($pathogen_gene_go_annotation_id, '$annot_ext_relation', '$annot_ext_value');
			);
    my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

  } # end foreach annotation extension

} # end foreach GO annotation




# iterate through all effector annotations
# inserting each into the database, associated with the
# relevant pathogen gene and PubMed ID
foreach my $effector_annot_ref (@effector_annotations) {

  # declare variables
  my $pubmed_id;
  my $pathogen_species;
  my $pathogen_taxon;
  my $pathogen_strain_id;
  my $pathogen_strain_name;
  my $effector_id;
  my $effector_evid_code;
  my $location_in_host_go_id;
  my $host_target_uniprot_acc;

  # increment effector annotation count
  $effector_annotation_count++;
 
  # get the hash of effector annotation details
  my %effector_annotation = %{ $effector_annot_ref };

  # get the type of effector annotation
  #my $annot_type = $effector_annotation{'type'};

  # get the pathogen gene identifier
  my $gene_id = $effector_annotation{'gene'};

  # get the gene from the hash of pathogen genes, using the ID
  my %gene = %{ $pathogen_genes{$gene_id} };
           
  # get the UniProt accession for the gene
  my $pathogen_gene_uniprot_acc = $gene{"uniquename"};

  # get the PubMed ID associated with the effector annotation
  $pubmed_id = $effector_annotation{'publication'};

  # get the annotation extensions for the effector annotation
  # which is an array of hashes (one hash for each extension)
  my @annot_extensions = @{ $effector_annotation{'extension'} };

  # from the annotation extensions in the list
  # find out if a specific pathogen strain taxon ID has been given
  # and/or a specific pathogen strain name
  # as well as identifying the other effector gene extensions
  foreach my $annot_ext_index (0 .. $#annot_extensions) {

    # get the hash for the annotation extension
    my %annot_ext_hash = %{ $annot_extensions[$annot_ext_index] };
    my $annot_ext_relation = $annot_ext_hash{'relation'}; 
    my $annot_ext_value = $annot_ext_hash{'rangeValue'}; 

    # if extension relation is 'pathogen_strain_id', then assign range value to pathogen_strain_id
    if ($annot_ext_relation =~ /^pathogen_strain_id/) {
       $pathogen_strain_id = $annot_ext_value;
    }

    # if extension relation is 'pathogen_strain_name', then assign range value to pathogen_strain_name
    if ($annot_ext_relation =~ /^pathogen_strain_name/) {
       $pathogen_strain_name = $annot_ext_value;
    }

    # if extension relation is 'location_in_host', then assign range value to location_in_host_go_id
    if ($annot_ext_relation =~ /^location_in_host/) {
       $location_in_host_go_id = $annot_ext_value;
    }

    # if extension relation is 'target_in_host', then assign range value to host_target_uniprot_acc
    if ($annot_ext_relation =~ /^target_in_host/) {
       $host_target_uniprot_acc = $annot_ext_value;
    }

  } # end foreach annotation extension

  # if a pathogen strain taxon ID has been explicitly given in the annotation extension
  # then this taxon ID should be used for all interactions, otherwise use the taxon ID
  # retrieved from the UniProt accessions supplied
  if (defined $pathogen_strain_id) {
    $pathogen_taxon_id = $pathogen_strain_id;
  }

  # insert data into the pathogen_gene table,
  # if it does not exist already
  # (based on combination of taxon id, UniProt accession, and where available the pathogen strain name)
  my $sql_statement2;
  if (defined $pathogen_strain_name) {
     $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
			  SELECT $pathogen_taxon_id,'$pathogen_strain_name','$pathogen_gene_uniprot_acc'
			  WHERE NOT EXISTS (
			    SELECT 1 FROM pathogen_gene
			    WHERE ncbi_taxon_id = $pathogen_taxon_id
			    AND pathogen_strain_name = '$pathogen_strain_name'
			    AND uniprot_accession = '$pathogen_gene_uniprot_acc'
			 ));
  } else {
     $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
			  SELECT $pathogen_taxon_id,'$pathogen_gene_uniprot_acc'
			  WHERE NOT EXISTS (
			    SELECT 1 FROM pathogen_gene
			    WHERE ncbi_taxon_id = $pathogen_taxon_id
			    AND uniprot_accession = '$pathogen_gene_uniprot_acc'
			 ));
  }

  my $sql_result2 = $db_conn->prepare($sql_statement2);
  $sql_result2->execute() or die $DBI::errstr;

  # get the unique identifier for the inserted pathogen_gene record
  my $sql_statement4 = qq(SELECT id FROM pathogen_gene
		          WHERE ncbi_taxon_id = $pathogen_taxon_id
		          AND uniprot_accession = '$pathogen_gene_uniprot_acc');

  my $sql_result4 = $db_conn->prepare($sql_statement4);
  $sql_result4->execute() or die $DBI::errstr;
  my @row4 = $sql_result4->fetchrow_array();
  my $pathogen_gene_id = shift @row4;

  $effector_id = $effector_annotation{'term'};
  $effector_evid_code = $effector_annotation{'evidence_code'};

  # insert data into effector_gene table,
  # with foreign keys to the pathogen_gene table, PubMed, the GO ontology, and UniProt,
  # depending on which values have been defined
  my $sql_statement; 
  if (defined $location_in_host_go_id and defined $host_target_uniprot_acc) {
     $sql_statement = qq(INSERT INTO effector_gene (pathogen_gene_id, pubmed_id, phi_effector_id, phi_effector_evidence_code, location_in_host_go_id, host_target_uniprot_acc)
			 VALUES ($pathogen_gene_id, '$pubmed_id', '$effector_id', '$effector_evid_code', '$location_in_host_go_id', '$host_target_uniprot_acc')
                        );
  } elsif (defined $location_in_host_go_id) {
     $sql_statement = qq(INSERT INTO effector_gene (pathogen_gene_id, pubmed_id, phi_effector_id, phi_effector_evidence_code, location_in_host_go_id)
			 VALUES ($pathogen_gene_id, '$pubmed_id', '$effector_id', '$effector_evid_code', '$location_in_host_go_id')
                        );
  } elsif (defined $host_target_uniprot_acc) {
     $sql_statement = qq(INSERT INTO effector_gene (pathogen_gene_id, pubmed_id, phi_effector_id, phi_effector_evidence_code, host_target_uniprot_acc)
			 VALUES ($pathogen_gene_id, '$pubmed_id', '$effector_id', '$effector_evid_code', '$host_target_uniprot_acc')
                        );
  } else { # neither location_in_host nor host_target defined
     $sql_statement = qq(INSERT INTO effector_gene (pathogen_gene_id, pubmed_id, effector_id, effector_evidence_code)
			 VALUES ($pathogen_gene_id, '$pubmed_id', '$effector_id', '$effector_evid_code')
			);
  }

  my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

} # end foreach effector annotation




print "\nTotal interactions:$interaction_count\n";
print "Total alleles:$allele_count\n";
print "Total annotations :$annotation_count\n";
print "Total PHI interaction phenotypes: $phi_interaction_phenotype_count\n";
print "Total PHI pathogen phenotypes: $phi_pathogen_phenotype_count\n";
print "Total PHI host phenotypes: $phi_host_phenotype_count\n";
print "Total GO annotations: $go_annotation_count\n";
print "Total effector annotations: $effector_annotation_count\n";
print "GO molecular function annotations: $go_mol_function_count\n";
print "GO biological process annotations: $go_biol_process_count\n";
print "GO cellular location annotations: $go_cell_location_count\n";

#close (JSON_OUTPUT_FILE);

$db_conn->disconnect() or die "Failed to disconnect database\n";

