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

my $json_filename = '../input/canto/canto_using_alleles.json';
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

=pod
  # iterate through all annotations, to find out which ones
  # are associated with the current interaction
  foreach my $pathogen_gene_key (keys %pathogen_genes) {

      print "Pathogen Gene Key:$pathogen_gene_key\n";

      my $path_organism = $pathogen_genes{$pathogen_gene_key}{"organism"};
      print "Pathogen Gene organism:$path_organism\n";

      my $path_unique_name = $pathogen_genes{$pathogen_gene_key}{"uniquename"};
      print "Pathogen Gene uniquename:$path_unique_name\n";

  }
=cut

# get the pathogen alleles for the session
my %pathogen_alleles = %{ $text_response->{'curation_sessions'}{$session_id}{'alleles'} };

=pod
  # iterate through all annotations, to find out which ones
  # are associated with the current interaction
  foreach my $pathogen_allele_key (keys %pathogen_alleles) {

      print "Pathogen Allele Key:$pathogen_allele_key\n";

      my $allele_identifier = $pathogen_alleles{$pathogen_allele_key}{"primary_identifier"};
      print "Pathogen Allele Identifier:$allele_identifier\n";

      my $allele_name = $pathogen_alleles{$pathogen_allele_key}{"name"};
      print "Pathogen Allele Name:$allele_name\n";

      my $allele_gene = $pathogen_alleles{$pathogen_allele_key}{"gene"};
      print "Pathogen Allele Gene:$allele_gene\n";

      my $allele_type = $pathogen_alleles{$pathogen_allele_key}{"allele_type"};
      print "Pathogen Allele Type:$allele_type\n";

      my $allele_desc = $pathogen_alleles{$pathogen_allele_key}{"description"};
      print "Pathogen Allele Description:$allele_desc\n";

  }
=cut

# get the pathogen genotypes for the session
my %pathogen_genotypes = %{ $text_response->{'curation_sessions'}{$session_id}{'genotypes'} };

=pod
  # iterate through all annotations, to find out which ones
  # are associated with the current interaction
  foreach my $pathogen_genotype_key (keys %pathogen_genotypes) {

      print "Pathogen Genotype Key:$pathogen_genotype_key\n";

      my @genotype_alleles = @{ $pathogen_genotypes{$pathogen_genotype_key}{"alleles"} };

      # iterate through the array of alleles in the genotype
      foreach my $genotype_allele_index (0 .. $#genotype_alleles) {
      
        my %allele = %{ $genotype_alleles[$genotype_allele_index] };

        my $allele_id = $allele{"id"};
        print "Genotype allele id:$allele_id\n";

        my $allele_expression = $allele{"expression"};
        print "Genotype allele expression:$allele_expression\n";

      }

  }
=cut

# annotations are an array of hashes
my @annotations = @{ $text_response->{'curation_sessions'}{$session_id}{'annotations'} };

# create a hash of hashes to store the pathogen-host interactions
# where the key will be an identifier for each interaction
# and the value will be a hash of the values that identify each unique interaction
my %interaction_profiles;

# create a separate hash of hashes to store the GO annotations
# (which are not part of the pathogen-host interaction, but a direct property of the pathogen gene)
# where the key will be an identifier for each GO annotation
# and the value will be a hash of the values that identify each unique interaction
my @go_annotations;

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
# such as the genotype identifier and annotation extensions
   my %annotation_interaction_hash = ();

   # array to hold all partner gene uniprot ids
   # (for multiple gene interactions)
#   my @interaction_partners;

   my %annotation = %{ $text_response->{'curation_sessions'}{$session_id}{'annotations'}[$annot_index] };

#   my @gene_organism_list = keys $annotation{'genes'};
#   my $gene_organism_name = $gene_organism_list[0];

#   my $uniprot_acc = $annotation{'genes'}{$gene_organism_name}{'uniquename'};
   # GO annotation is annotated against a gene, whereas phenotypes are annotated against a genotype
   # so need to get correct annotation type for the identifier
   my $gene_or_genotype_id;
   my $annot_type = $annotation{'type'};
   if ($annot_type eq 'phi_interaction_phenotype' or $annot_type eq 'phi_pathogen_phenotype' or $annot_type eq 'phi_host_phenotype') {
      $gene_or_genotype_id = $annotation{'genotype'};
   } elsif ($annot_type eq 'molecular_function' or $annot_type eq 'biological_process' or $annot_type eq 'cellular_component' ) {
      $gene_or_genotype_id = $annotation{'gene'};
   } else {
      print STDERR "Error: Invalid annotation type: $annot_type\n"; 
   }

   # add uniprot accession to help identify the correct interaction
#   $annotation_interaction_hash{'pathogen_uniprot_acc'} = $uniprot_acc;
   $annotation_interaction_hash{'gene_or_genotype_id'} = $gene_or_genotype_id;

   my $annot_extension_list = $annotation{'annotation_extension'};
   my @annot_extensions = split(/,/,$annot_extension_list);

   my $host_taxon;
   my $tissue;
   my $interaction_partner;
   my $disease;
   my $inducer;
   my $anti_infective;
   my $host_interacting_protein;
   my $pathogen_interacting_protein;
   my $experiment_spec;
   my $host_response_id;
   my $pathogen_strain_id;

   # identify each annotation extension in the list
   foreach my $annot_ext (@annot_extensions) {

     # if the annotation extension begins with 'experimental_host', then assign value between brackets to host
     if ($annot_ext =~ /^experimental_host/) {
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

     # if the annotation extension begins with 'causes_disease', then assign value between brackets to disease
     if ($annot_ext =~ /^causes_disease/) {
	my @annot_ext = split(/[\(\)]/,$annot_ext);
	$disease = $annot_ext[1];
	# add annotation extension to help identify the correct interaction
	$annotation_interaction_hash{'disease'} = $disease;
     }

=pod
# REPLACED WITH PHI PHENOTYPE ONTOLOGY FOR HOST PHENOTYPES
     # if the annotation extension begins with 'host_responsee', then assign value between brackets to host_response
     if ($annot_ext =~ /^host_response/) {
	my @annot_ext = split(/[\(\)]/,$annot_ext);
	$host_response_id = $annot_ext[1];
	# add annotation extension to help identify the correct interaction
	$annotation_interaction_hash{'host_response_id'} = $host_response_id;
     }
=cut

     # if the annotation extension begins with 'induced_by', then assign value between brackets to inducer
     if ($annot_ext =~ /^induced_by/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $inducer = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'inducer'} = $inducer;
     }

     # if the annotation extension begins with 'anti_infective', then assign value between brackets to anti_infective
     if ($annot_ext =~ /^anti_infective/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $anti_infective = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'anti_infective'} = $anti_infective;
     }

     # if the annotation extension begins with 'host_interacting_protein',
     # then assign value between brackets to host_interacting_protein
     if ($annot_ext =~ /^host_interacting_protein/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $host_interacting_protein = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'host_interacting_protein'} = $host_interacting_protein;
     }

     # if the annotation extension begins with 'pathogen_interacting_protein',
     # then assign value between brackets to pathogen_interacting_protein
     if ($annot_ext =~ /^pathogen_interacting_protein/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $pathogen_interacting_protein = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'pathogen_interacting_protein'} = $pathogen_interacting_protein;
     }

     # if the annotation extension begins with 'pathogen_strain', then assign value between brackets to pathogen_strain_id
     if ($annot_ext =~ /^pathogen_strain/) {
        my @annot_ext = split(/[\(\)]/,$annot_ext);
        $pathogen_strain_id = $annot_ext[1];
        # add annotation extension to help identify the correct interaction
        $annotation_interaction_hash{'pathogen_strain_id'} = $pathogen_strain_id;
     }


=pod
# REPLACED WITH GENOTYPE DESCRIPTION BASED ON COMBINATION OF ALLELES
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
=cut


   } # end foreach annotation extension


   # declare flag to indicate if the interaction associated with
   # the current annotation has been found
   my $interaction_found = 0;

   # iterate through each of the existing interactions to find out
   # if the current annotation profile matches an existing interaction
   # if so, then add the corresponding interaction id to the annotation
   # otherwise, create a new interaction profile for the annotation
   if ($annot_type eq 'phi_interaction_phenotype' or $annot_type eq 'phi_pathogen_phenotype' or $annot_type eq 'phi_host_phenotype') {

      print "ANNOTATION IS GENOTYPE - INTERACTION\n";

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

   } elsif ($annot_type eq 'molecular_function' or $annot_type eq 'biological_process' or $annot_type eq 'cellular_component' ) {

      print "ANNOTATION IS GENE - GO ANNOTATION\n";

      # add the current annotation to the array of GO annotations;
      push (@go_annotations, $annotations[$annot_index]);

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
  my $uniprot_acc;
  my $pubmed_id;
  my $host_taxon;
  my $host_taxon_id;
  my $tissue;
  my $interaction_partner;
  my $curator_name;
  my $curator_email;
  my $pathogen_species;
  my $host_response_id;
  my $pathogen_taxon;
  my $phenotype_outcome;
  my $disease;
  my $go_id;
  my $go_evid_code;
  my $interaction_id;
  my $experiment_spec;
  my $inducer;
  my $anti_infective;
  my $host_interacting_protein;
  my $pathogen_interacting_protein;
  my $pathogen_strain_id;
  my $creation_date;

  my $gene_or_genotype_id;

  # declare array to hold all pathogen gene
  # involved in the interaction
  my @pathogen_genes;

  # declare array to hold all pathogen gene alleles
  # involved in the interaction
  my @pathogen_genes_alleles;

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

         my @genotype_alleles = @{ $pathogen_genotypes{$gene_or_genotype_id}{"alleles"} };

         # iterate through the array of alleles in the genotype
         foreach my $genotype_allele_index (0 .. $#genotype_alleles) {
      
	   my %allele = %{ $genotype_alleles[$genotype_allele_index] };

           

# array of pathogen_genes for the current interaction
  # each element is a hash of pathogen_gene_alleles,
  # with the key being the UniProt ID of the gene
     # the hash is an array of allele details, with each element listing the name, type, & description of each allele
# pathogen_genes [
#      QSD4F: [
#          {
#              allele_name: QSD4F+
#              allele_type: wild_type
#              allele_desc: Wild type of ...
#          },
#          {
#              allele_name: QSD4Fdelta
#              allele_type: deletion
#              allele_desc: deletion of ...
#          },   
#      ],
#      QSD4T: [
#          {
#              allele_name: QSD4T+
#              allele_type: wild_type
#              allele_desc: Wild type of ...
#          },
#          {
#              allele_name: QSD4Tdelta
#              allele_type: deletion
#              allele_desc: deletion of ...
#          },   
#      ],

#  ]


            # find the alleles relevant to the current gene
            
            # foreach allele in genotype
              # retrieve the id
              # use the id to get the allele hash from the full list of alleles
              # get the gene identifier from the allele hash
              # get the relevant gene hash from the full list of genes
              # get the uniquename for the gene
              # compare the uniquename with the uniprot identifier
                 # if uniprot acc eq current uniprot id
                   # retreive other allele details from the allele hash - name, allele_type, description
                   # use allele details + gene id to insert a new pathogen_gene_allele record.

	   my $allele_id = $allele{"id"};
	   print "Genotype allele id for annot:$allele_id\n";

	   my $allele_expression = $allele{"expression"};
	   print "Genotype allele expression for annot:$allele_expression\n";

	   my $allele_identifier = $pathogen_alleles{$allele_id}{"primary_identifier"};
	   print "Pathogen Allele Identifier for annot:$allele_identifier\n";

	   my $allele_name = $pathogen_alleles{$allele_id}{"name"};
	   print "Pathogen Allele Name for annot:$allele_name\n";

	   my $allele_gene = $pathogen_alleles{$allele_id}{"gene"};
	   print "Pathogen Allele Gene for annot:$allele_gene\n";

	   my $allele_type = $pathogen_alleles{$allele_id}{"allele_type"};
	   print "Pathogen Allele Type for annot:$allele_type\n";

	   my $allele_desc = $pathogen_alleles{$allele_id}{"description"};
	   print "Pathogen Allele Description for annot:$allele_desc\n";

	   my $path_organism = $pathogen_genes{$allele_gene}{"organism"};
	   print "Pathogen Gene organism for annot:$path_organism\n";

	   #my $path_uniprot_acc = $allele{"uniquename"};
	   my $path_uniprot_acc = $pathogen_genes{$allele_gene}{"uniquename"};
	   print "Pathogen Gene UniProt accession for annot:$path_uniprot_acc\n";

           # find if the current gene is already in the gene list,
           # if it is, then add the allele details to the existing pathogen_gene element (add as a hash)
           
           # if not, then add the allele info to a new pathogen gene

	   # add the uniprot accession for this gene to the array
	   # of genes, which includes all interaction partner genes
	   # (for multiple gene interaction) and the original gene
	   push(@pathogen_genes,$path_uniprot_acc);

	 } # end foreach genotype allele

         print "Genotype ID: $gene_or_genotype_id\n";

         my $annot_extension_list = $annotation{'annotation_extension'};
         my @annot_extensions = split(/,/,$annot_extension_list);

         # identify each annotation extension in the list
         foreach my $annot_ext (@annot_extensions) {

            # if the annotation extension begins with 'experimental_host', then assign value between brackets to host
            if ($annot_ext =~ /^experimental_host/) {
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
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $disease = $annot_ext[1];
            }

=pod
            # if the annotation extension begins with 'host_responsee', then assign value between brackets to host_response
            if ($annot_ext =~ /^host_response/) {
              my @annot_ext = split(/[\(\)]/,$annot_ext);
              $host_response_id = $annot_ext[1];
            }
=cut

	    # if the annotation extension begins with 'induced_by', then assign value between brackets to inducer
	    if ($annot_ext =~ /^induced_by/) {
	       my @annot_ext = split(/[\(\)]/,$annot_ext);
	       $inducer = $annot_ext[1];
	    }

	    # if the annotation extension begins with 'anti_infective', then assign value between brackets to anti_infective
	    if ($annot_ext =~ /^anti_infective/) {
	       my @annot_ext = split(/[\(\)]/,$annot_ext);
	       $anti_infective = $annot_ext[1];
	    }

	    # if the annotation extension begins with 'host_interacting_protein',
	    # then assign value between brackets to host_interacting_protein
	    if ($annot_ext =~ /^host_interacting_protein/) {
	       my @annot_ext = split(/[\(\)]/,$annot_ext);
	       $host_interacting_protein = $annot_ext[1];
	       # UniProt accession is the interacting protein with 'UniProt:' prefix removed
	       $host_interacting_protein = substr($host_interacting_protein,8);
	    }

	    # if the annotation extension begins with 'pathogen_interacting_protein',
	    # then assign value between brackets to pathogen_interacting_protein
	    if ($annot_ext =~ /^pathogen_interacting_protein/) {
	       my @annot_ext = split(/[\(\)]/,$annot_ext);
	       $pathogen_interacting_protein = $annot_ext[1];
	       # UniProt accession is the interacting protein with 'UniProt:' prefix removed
	       $pathogen_interacting_protein = substr($pathogen_interacting_protein,8);
	    }

	    # if the annotation extension begins with 'pathogen_strain', then assign value between brackets to pathogen_strain_id
	    if ($annot_ext =~ /^pathogen_strain/) {
	       my @annot_ext = split(/[\(\)]/,$annot_ext);
	       $pathogen_strain_id = $annot_ext[1];
	    }

=pod
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
=cut

         } # end foreach annotation extension


         foreach my $pathogen_gene (@pathogen_genes) {
           print "Pathogen Gene: $pathogen_gene\n";
         }

         print "Host Taxon: $host_taxon\n";
         print "Tissue: $tissue\n" if defined $tissue;
         print "Disease: $disease\n" if defined $disease;
#         print "Host Response: $host_response_id\n";
         print "Inducer: $inducer\n" if defined $inducer;
         print "Anti-infective: $anti_infective\n" if defined $anti_infective;
         print "Host Interaction protein: $host_interacting_protein\n" if defined $host_interacting_protein;
         print "Pathogen Interaction protein: $pathogen_interacting_protein\n" if defined $pathogen_interacting_protein;
         print "Pathogen Strain ID: $pathogen_strain_id\n" if defined $pathogen_strain_id;

         $curator_name = $annotation{'curator'}{'name'};
         $curator_email = $annotation{'curator'}{'email'};
         print "Curator: $curator_name, $curator_email\n";
         print "Approver: $approver_name, $approver_email\n";

#         $pathogen_species = $annotation{'genes'}{$gene_organism_name}{'organism'};
#         print "Pathogen Species: $pathogen_species\n";

         $creation_date = $annotation{'creation_date'};
         print "Creation Date: $creation_date\n";


         # declare array to hold the list of IDs from the pathogen_gene_mutant table 
         my @pathogen_gene_mutant_ids;

         # if a pathogen strain taxon ID has been explicitly given in the annotation extension
         # then this taxon ID should be used for all interactions, otherwise use the taxon ID
         # retrieved from the UniProt accessions supplied
         if (defined $pathogen_strain_id) {
	    # pathogen taxon id becomes the strain taxon with 'NCBItaxon:' prefix removed
            $pathogen_taxon_id = substr($pathogen_strain_id,10);
         }


	 # insert data into the pathogen_gene table, for each pathogen gene
	 # if it does not exist already (based on combination of taxon id and UniProt accession)
	 foreach my $pathogen_gene_uniprot_acc (@pathogen_genes) {

	    my $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
				    SELECT $pathogen_taxon_id,'$pathogen_gene_uniprot_acc'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM pathogen_gene
				      WHERE ncbi_taxon_id = $pathogen_taxon_id
				      AND uniprot_accession = '$pathogen_gene_uniprot_acc'
				   ));

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


            # find the alleles relevant to the current gene
            
            # foreach allele in the full list of alleles
              # retrieve the id
              # get the gene identifier from the allele hash
              # get the relevant gene hash from the full list of genes
              # get the uniquename for the gene
              # compare the uniquename with the uniprot identifier
                 # if uniprot acc eq current uniprot id
                   # retreive other allele details from the allele hash - name, allele_type, description
                   # use allele details + gene id to insert a new pathogen_gene_allele record.


            foreach my $allele_id (keys %pathogen_alleles) {

               # get the hash of allele details
               my %pathogen_allele = %{ $pathogen_alleles{$allele_id} };

               my $pathogen_allele_gene_id = $pathogen_allele{"gene"};
               print "Pathogen Allele Gene ID:$pathogen_allele_gene_id\n";

               # get the hash of gene details from the pathogen_genes hash
               my %pathogen_allele_gene = %{ $pathogen_genes{$pathogen_allele_gene_id} };

               my $pathogen_allele_gene_uniprot_acc = $pathogen_allele_gene{"uniquename"};
               print "Pathogen Allele Gene UniProt:$pathogen_allele_gene_uniprot_acc\n";

               if ($pathogen_allele_gene_uniprot_acc eq $pathogen_gene_uniprot_acc) {
                    # now we know that this allele needs to be added as a pathogen_gene_allele record
                    # for the current pathogen gene

                   # retreive other allele details from the allele hash - name, allele_type, description
                   # use allele details + gene id to insert a new pathogen_gene_allele record.

                   my $gene_allele_name = $pathogen_allele{"name"};
                   my $gene_allele_type = $pathogen_allele{"allele_type"};
                   my $gene_allele_desc = $pathogen_allele{"description"};

                   print "Allele name for current gene:$gene_allele_name\n";
                   print "Allele type for current gene:$gene_allele_type\n";
                   print "Allele desc for current gene:$gene_allele_desc\n";

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

               }

            } # end foreach pathogen allele


=pod

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
=cut
         } # end foreach pathogen gene UniProt identifier
=pod

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
            # if the host interacting protein is defined, then this is also added to the interaction_host table,
            # otherwise, it is left out of the SQL insert statement
            my $inner_sql_statement;
            if (defined $host_interacting_protein) {
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id,first_target_uniprot_accession) 
					   VALUES ($interaction_id,$host_taxon_id,'$host_interacting_protein');
					);
            } else {
	       $inner_sql_statement = qq(
					 INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
					   VALUES ($interaction_id,$host_taxon_id);
					);
            }
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


         if (defined $host_response_id) {

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


         if (defined $inducer) {

	    # insert data into the chemical_table,
	    # if it does not exist already (based on the ChEBI ID)
	    my $sql_statement = qq(INSERT INTO chemical (chebi_id) 
				      SELECT '$inducer'
				    WHERE NOT EXISTS (
				      SELECT 1 FROM chemical
				      WHERE chebi_id = '$inducer'
				    )
				  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;

	    # before we can insert the inducer chemical,
	    # we need to get the identifier for the chemical
	    $sql_statement = qq(SELECT id FROM chemical
				     WHERE chebi_id = '$inducer'
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

         } # end if inducer


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


 @genotype_alleles;
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
=cut
#      if ($annot_type eq 'pathogen_taxon') {
#         $pathogen_taxon = $annot_value;
#         print "Pathogen Taxon: $pathogen_taxon\n";
#         # the pathogen taxon is split between the name and value, separated by colon
#         # need to extract the taxon value
#         my @pathogen_taxon_parts = split(/[:]/,$pathogen_taxon);
#         $pathogen_taxon_id = $pathogen_taxon_parts[1];
#         print "Pathogen Taxon ID: $pathogen_taxon_id\n";
#      } elsif ($annot_type eq 'phenotype_outcome') {
=pod
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


my $go_mol_function_count = 0;
my $go_biol_process_count = 0;
my $go_cell_location_count = 0;

# iterate through all GO annotations
foreach my $go_annot_ref (@go_annotations) {

  # declare variables
  my $uniprot_acc;
  my $pubmed_id;
  my $pathogen_species;
  my $pathogen_taxon;
  my $go_id;
  my $go_evid_code;
  my $pathogen_strain_id;

  $go_annotation_count++;
 
  print "IN GO ANNOTATION REF:$go_annot_ref\n";

  # get the hash of GO annotation details
  my %go_annotation = %{ $go_annot_ref };

  my $annot_type = $go_annotation{'type'};
  print "GO ANNOTATION TYPE:$annot_type\n";
  if ($annot_type eq 'molecular_function') {
    $go_mol_function_count++;
  } elsif ($annot_type eq 'biological_process') {
    $go_biol_process_count++;
  } elsif ($annot_type eq 'cellular_component' ) {
    $go_cell_location_count++;
  }

  my $gene_id = $go_annotation{'gene'};

  #get the gene from the hash of pathogen genes, using the ID
  my %gene = %{ $pathogen_genes{$gene_id} };
           
  my $path_organism = $gene{"organism"};
  print "Pathogen Gene organism for annot:$path_organism\n";

  my $pathogen_gene_uniprot_acc = $gene{"uniquename"};
  print "Pathogen Gene UniProt Acc for annot:$pathogen_gene_uniprot_acc\n";

  $pubmed_id = $go_annotation{'publication'};
  print "PubMed Article: $pubmed_id\n";

  my $annot_extension_list = $go_annotation{'annotation_extension'};
  my @annot_extensions = split(/,/,$annot_extension_list);

  # identify each annotation extension in the list
  # and find out if a specific pathogen strain taxon ID has been given
  foreach my $annot_ext (@annot_extensions) {

    # if the annotation extension begins with 'pathogen_strain', then assign value between brackets to pathogen_strain_id
    if ($annot_ext =~ /^pathogen_strain/) {
       my @annot_ext = split(/[\(\)]/,$annot_ext);
       $pathogen_strain_id = $annot_ext[1];
    }

  } # end foreach annotation extension

  print "Pathogen Strain ID: $pathogen_strain_id\n" if defined $pathogen_strain_id;

  # if a pathogen strain taxon ID has been explicitly given in the annotation extension
  # then this taxon ID should be used for all interactions, otherwise use the taxon ID
  # retrieved from the UniProt accessions supplied
  if (defined $pathogen_strain_id) {
    # pathogen taxon id becomes the strain taxon with 'NCBItaxon:' prefix removed
    $pathogen_taxon_id = substr($pathogen_strain_id,10);
  }


  # insert data into the pathogen_gene table,
  # if it does not exist already (based on combination of taxon id and UniProt accession)
  my $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id, uniprot_accession) 
			  SELECT $pathogen_taxon_id,'$pathogen_gene_uniprot_acc'
			  WHERE NOT EXISTS (
			    SELECT 1 FROM pathogen_gene
			    WHERE ncbi_taxon_id = $pathogen_taxon_id
			    AND uniprot_accession = '$pathogen_gene_uniprot_acc'
			 ));

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

  print "GO Annotation: $go_id ($go_evid_code)\n";

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

} # end foreach GO annotation


print "\nTotal interactions:$interaction_count\n";
print "Total annotations :$annotation_count\n";
print "Total phenotype outcomes: $phenotype_outcome_count\n";
#print "Total diseases: $disease_count\n";
print "Total GO annotations: $go_annotation_count\n";
print "GO molecular function annotations: $go_mol_function_count\n";
print "GO biological process annotations: $go_biol_process_count\n";
print "GO cellular location annotations: $go_cell_location_count\n";
#print "Total host responses: $host_response_count\n";

#close (JSON_OUTPUT_FILE);

$db_conn->disconnect() or die "Failed to disconnect database\n";

