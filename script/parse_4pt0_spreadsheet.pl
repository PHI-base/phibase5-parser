#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module
use Array::Utils qw(:all);

use phibase_subroutines qw(connect_to_phibase query_uniprot ontology_mapping); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database


# parse text file that maps columns headings of the spreadsheet to database field names
# saving the column name and db field name as key/value pairs in a hash
open (COL_NAMES_FILE, "../input/column2accession_4pt0.txt") || die "Error opening input file\n";

# hash to map spreadsheet headings to PHI-base3 db field names
my %column_mapping;

# each row of the file contains a "spreadsheet heading" and the corresponding "db field name", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the spreadsheet heading & value becomes db field name
while (<COL_NAMES_FILE>) {
  chomp;
  my ($sheet_heading,$db_field) = split(/\t/,$_);
  $column_mapping{$sheet_heading} = $db_field; 
  # create or refresh output text files for each column
  open (COL_TEXT_FILE, "> ../output/column/spreadsheet_column_$db_field.txt") or die "Error opening output file column_$db_field\n";
  close (COL_TEXT_FILE);
}
close (COL_NAMES_FILE);

# Read in the relevant ontologies mapping files
print "Reading ontology mapping files...\n";

# parse tab-separated file that maps experimental evidence values of the spreadsheet
# to identifiers in the experiment specification ontology
# saving the value and identifiers as key/value pairs in a hash
open (EXP_SPEC_MAPPINGS_FILE, "../mapping/experiment_spec_mapping_phibase_4pt0.tsv") || die "Error opening input file\n";

# hash to map experimental evidence to ontology identifiers
my %exp_spec_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "experiment spec ontology identifiers", separated by tab
# mutliple ontology identifiers are separated by semi-colon
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes ontology identifiers
# (note that the identifiers themselves will be separated later)
while (<EXP_SPEC_MAPPINGS_FILE>) {
  chomp;
  my ($exp_spec_value,$exp_spec_ontology_id_list) = split(/\t/,$_);
  $exp_spec_mapping{$exp_spec_value} = $exp_spec_ontology_id_list;
}
close (EXP_SPEC_MAPPINGS_FILE);


# parse tab-separated file that maps host response values of the spreadsheet
# to identifiers in the PHI phenotype ontology
# saving the value and identifiers as key/value pairs in a hash
open (HOST_RES_MAPPINGS_FILE, "../mapping/host_response_mapping_phibase_4pt0.tsv") || die "Error opening input file\n";

# hash to map host response text to ontology identifiers
my %host_response_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "phi phenotype ontology identifiers", separated by tab
# mutliple ontology identifiers are separated by semi-colon
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes ontology identifiers
# (note that the identifiers themselves will be separated later)
while (<HOST_RES_MAPPINGS_FILE>) {
  chomp;
  my ($host_res_value,$host_res_ontology_id_list) = split(/\t/,$_);
  $host_response_mapping{$host_res_value} = $host_res_ontology_id_list;
}
close (HOST_RES_MAPPINGS_FILE);


# parse tab-separated file that maps phenotype outcome values of the spreadsheet
# to the identifier in the PHI phenotype ontology
# saving the value and identifier as key/value pairs in a hash
open (PHEN_OUTCOME_MAPPINGS_FILE, "../mapping/phenotype_outcome_mapping_phibase_4pt0.tsv") || die "Error opening input file\n";

# hash to map phenotype outcome text to ontology identifier
my %phenotype_outcome_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "PHI phenotype ontology identifier", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes the ontology identifier
while (<PHEN_OUTCOME_MAPPINGS_FILE>) {
  chomp;
  my ($phen_outcome_value,$phen_outcome_ontology_id_list) = split(/\t/,$_);
  $phenotype_outcome_mapping{$phen_outcome_value} = $phen_outcome_ontology_id_list;
}
close (PHEN_OUTCOME_MAPPINGS_FILE);

=pod
NO LONGER QUERYING DISEASE ONTOLOGY DIRECTLY - USING PRE-COMPOSED MAPPING FILE INSTEAD
(JUST LIKE THE OTHER ONTOLOGIES)
# disease terms for PHI-base should include the combined terms from
# both the Plant Disease ontology and the [Human] Disease Ontology,
# so a hash of key/value pairs is created where the key is the ontology identifier
# and the value is the term name
my %combined_disease_mapping = (
      ontology_mapping('../ontology/Disease/PlantDisease/plant_disease_ontology.obo'),
      ontology_mapping('../ontology/Disease/HumanDisease/doid.obo')
   );
=cut

# parse tab-separated file that maps disease values of the spreadsheet
# to identifiers in the [Human] Disease Ontology or Plant Disease Ontology
# saving the value and identifiers as key/value pairs in a hash
open (DISEASE_MAPPINGS_FILE, "../mapping/disease_mapping_phibase_4pt0.tsv") || die "Error opening input file\n";

# hash to map disease to ontology identifiers
my %disease_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "ontology identifiers", separated by tab (for either the human or plant disease ontologies)
# mutliple ontology identifiers are separated by semi-colon
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes ontology identifiers
# (note that the identifiers themselves will be separated later)
while (<DISEASE_MAPPINGS_FILE>) {
  chomp;
  my ($disease_value,$disease_ontology_id_list) = split(/\t/,$_);
  $disease_mapping{$disease_value} = $disease_ontology_id_list;
}
close (DISEASE_MAPPINGS_FILE);


# parse tab-separated file that maps tissue values of the spreadsheet
# to identifiers in the BRENDA tissue ontology
# saving the value and identifiers as key/value pairs in a hash
open (TISSUE_MAPPINGS_FILE, "../mapping/tissue_mapping_phibase_4pt0.tsv") || die "Error opening input file\n";

# hash to map tissue text to ontology identifiers
my %tissue_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "tissue ontology identifiers", separated by tab
# mutliple ontology identifiers are separated by semi-colon
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes ontology identifiers
# (note that the identifiers themselves will be separated later)
while (<TISSUE_MAPPINGS_FILE>) {
  chomp;
  my ($tissue_value,$tissue_ontology_id_list) = split(/\t/,$_);
  $tissue_mapping{$tissue_value} = $tissue_ontology_id_list;
}
close (TISSUE_MAPPINGS_FILE);


# open the tab separated values (TSV) version of the PHI-base spreadsheet
#my $phibase_tsv_filename = '../input/phi-base-1_vs36_reduced_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_to_v3_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_columns_from_phi560.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_columns_from_phi2060.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_columns_from_phi2661.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_columns_until_phi3043.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2014-12-11_reduced_cols_with_record_id.tsv';
#my $phibase_tsv_filename = '../input/phi4_2015-02-04_reduced_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2015-02-04_reduced_columns_fixed_mult_mut.tsv';
#my $phibase_tsv_filename = '../input/phi4_2015-04-18_reduced_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2015-09-24_PHI4pt0_reduced_columns.tsv';
#my $phibase_tsv_filename = '../input/phi4_2016-01-15_PHI4pt1_reduced_columns.tsv';
my $phibase_tsv_filename = '../input/phi4_2016-02-03_PHI4pt1_reduced_columns.tsv';
open (TSV_FILE, $phibase_tsv_filename) || die "Error opening input file\n";
print "Processing PHI-base data from $phibase_tsv_filename...\n";
print "Inserting data for valid annotations into PHI-base v5 database...\n";

# open output files
my $pathogen_phenotype_filename = '../output/phibase_pathogen_phenotypes.tsv';
my $invalid_pathogen_phenotype_filename = '../error/phibase_invalid_pathogen_phenotypes.tsv';
my $go_with_evid_filename = '../output/phibase_go_with_evid.tsv';
my $go_without_evid_filename = '../output/phibase_go_without_evid.tsv';
my $invalid_go_filename = '../error/phibase_invalid_go.tsv';
my $exp_spec_term_filename = '../output/phibase_exp_spec_terms.tsv';
my $invalid_exp_spec_filename = '../error/phibase_invalid_exp_specs.tsv';
my $host_res_term_filename = '../output/phibase_host_response_terms.tsv';
my $invalid_host_res_filename = '../error/phibase_invalid_host_responses.tsv';
my $phen_outcome_term_filename = '../output/phibase_phenotype_outcome_terms.tsv';
my $invalid_phen_outcome_filename = '../error/phibase_invalid_phenotype_outcomes.tsv';
my $disease_term_filename = '../output/phibase_disease_terms.tsv';
my $invalid_disease_filename = '../error/phibase_invalid_diseases.tsv';
my $without_disease_filename = '../error/phibase_without_diseases.tsv';
my $tissue_term_filename = '../output/phibase_tissue_terms.tsv';
my $invalid_tissue_filename = '../error/phibase_invalid_tissues.tsv';
my $invalid_required_data_filename = '../error/phibase_invalid_required_data.tsv';
my $invalid_path_taxon_filename = '../error/phibase_invalid_path_taxon_ids.tsv';
my $invalid_host_taxon_filename = '../error/phibase_invalid_host_taxon_ids.tsv';
my $invalid_literature_filename = '../error/phibase_invalid_literature.tsv';
my $invalid_uniprot_filename = '../error/phibase_invalid_uniprot_acc.tsv';
my $invalid_gene_name_filename = '../error/phibase_invalid_gene_names.tsv';
my $invalid_curator_filename = '../error/phibase_invalid_curators.tsv';
open (PATH_PHENOTYPE_FILE, "> $pathogen_phenotype_filename") or die "Error opening output file\n";
open (INVALID_PATH_PHENOTYPE_FILE, "> $invalid_pathogen_phenotype_filename") or die "Error opening output file\n";
open (GO_WITH_EVID_FILE, "> $go_with_evid_filename") or die "Error opening output file\n";
open (GO_WITHOUT_EVID_FILE, "> $go_without_evid_filename") or die "Error opening output file\n";
open (INVALID_GO_FILE, "> $invalid_go_filename") or die "Error opening output file\n";
open (EXP_SPEC_TERM_FILE, "> $exp_spec_term_filename") or die "Error opening output file\n";
open (INVALID_EXP_SPEC_FILE, "> $invalid_exp_spec_filename") or die "Error opening output file\n";
open (HOST_RES_TERM_FILE, "> $host_res_term_filename") or die "Error opening output file\n";
open (INVALID_HOST_RES_FILE, "> $invalid_host_res_filename") or die "Error opening output file\n";
open (PHEN_OUTCOME_TERM_FILE, "> $phen_outcome_term_filename") or die "Error opening output file\n";
open (INVALID_PHEN_OUTCOME_FILE, "> $invalid_phen_outcome_filename") or die "Error opening output file\n";
open (DISEASE_TERM_FILE, "> $disease_term_filename") or die "Error opening output file\n";
open (INVALID_DISEASE_FILE, "> $invalid_disease_filename") or die "Error opening output file\n";
open (WITHOUT_DISEASE_FILE, "> $without_disease_filename") or die "Error opening output file\n";
open (TISSUE_TERM_FILE, "> $tissue_term_filename") or die "Error opening output file\n";
open (INVALID_TISSUE_FILE, "> $invalid_tissue_filename") or die "Error opening output file\n";
open (INVALID_REQUIRED_DATA_FILE, "> $invalid_required_data_filename") or die "Error opening output file\n";
open (INVALID_PATH_TAXON_FILE, "> $invalid_path_taxon_filename") or die "Error opening output file\n";
open (INVALID_HOST_TAXON_FILE, "> $invalid_host_taxon_filename") or die "Error opening output file\n";
open (INVALID_LITERATURE_FILE, "> $invalid_literature_filename") or die "Error opening output file\n";
open (INVALID_UNIPROT_FILE, "> $invalid_uniprot_filename") or die "Error opening output file\n";
open (INVALID_GENE_NAME_FILE, "> $invalid_gene_name_filename") or die "Error opening output file\n";
open (INVALID_CURATOR_FILE, "> $invalid_curator_filename") or die "Error opening output file\n";

# first line gives the spreadsheet column headings
chomp(my $col_header_line = <TSV_FILE>);
my @col_headers = split(/\t/,$col_header_line);

# array to store invalid PHI-base accessions
my @invalid_phibase_acc;

# create another hash for the subset with only required fields
my %required_fields_data;
# create another hash for all annotation that meet the required data criteria
my %required_criteria_annotations;

# counter for interaction number (used for new PHI-base accession)
my $interaction_num = 0;

# counters to gather statistics
my $annotation_count = 0;
my $pathogen_phenotype_count = 0;
my $invalid_pathogen_phenotype_count = 0;
my $curator_count = 0;
my $species_expert_count = 0;
my $anti_infective_count = 0;
my $no_anti_infective_count = 0;
my $inducer_count = 0;
my $no_inducer_count = 0;
my $go_annotation_count = 0;
my $go_with_evid_count = 0;
my $go_without_evid_count = 0;
my $invalid_go_count = 0;
my $exp_spec_count = 0;
my $exp_spec_term_count = 0;
my $invalid_exp_spec_count = 0;
my $host_response_count = 0;
my $host_response_term_count = 0;
my $invalid_host_response_count = 0;
my $phenotype_outcome_count = 0;
my $phenotype_outcome_term_count = 0;
my $invalid_phenotype_outcome_count = 0;
my $disease_count = 0;
my $disease_term_count = 0;
my $invalid_disease_count = 0;
my $without_disease_count = 0;
my $tissue_count = 0;
my $tissue_term_count = 0;
my $invalid_tissue_count = 0;
my $invalid_required_data_count = 0;
my $invalid_path_taxon_id_count = 0;
my $invalid_host_taxon_id_count = 0;
my $invalid_literature_count = 0;
my $invalid_uniprot_acc_count = 0;
my $invalid_gene_name_count = 0;
my $invalid_curator_count = 0;

# go through each of the remaining lines of the TSV file (each representing a single annotation)
# save the values of each column to the approriate output file
while (<TSV_FILE>) {

   # increment annotation counter
   $annotation_count++;

   # remove newline from input string
   chomp($_);

   # each value is separated based on the tab, then saved as an element of the array
   my @phi_array = split(/\t/,$_);

   # initiate column iterator
   my $column_num=0;

   # hash to store all values of the current annotation
   my %phi_base_annotation;

   # iterate through each column of the annotation, saving its value to the appropriate text file
   # the name of the text file is determined by mapping the column header to db field name
   foreach my $phi_value (@phi_array) {

      #print "Col Value:$phi_value\n";
      #print "Col num:$column_num\n";
      #print "Col header:$col_headers[$column_num]\n";

      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ../output/column/spreadsheet_column_".$column_mapping{$col_headers[$column_num]}.".txt")
         or die "Error opening output file\n";
      print COLUMN_FILE "$phi_value\n";

      # add data to the annotation hash
      $phi_base_annotation{$column_mapping{$col_headers[$column_num]}} = $phi_value;

      # increment the column number
      $column_num++;

   } # end foreach column of the annotation

   # variable for PHI-base accession number
   my $phi_acc_num;

   # check if the PHI-base accession exists & if it has a valid prefix
   # then remove the PHI: prefix, leaving just the PHI-base accession number
   # and add the current annotation to the overall data hash, using the accession number as a key
   if ($phi_acc_num = $phi_base_annotation{"phi_base_acc"} and $phi_acc_num =~ /^PHI:/) {

     $phi_acc_num =~ s/PHI://;

     # get subset of the hash containing only the required fields
     my @required_fields = (
     "phi_base_acc",
     "db_type",
     "accession",
     "gene_name",
     "patho_tax",
     "host_tax",
     "literature_id",
     "literature_source",
     "entered_by"
     );
     my %required_fields_annot;
     @required_fields_annot{@required_fields} = @phi_base_annotation{@required_fields};
     $required_fields_data{$phi_acc_num} = {%required_fields_annot};


     # for UniProt IDs, remove any leading or trailing spaces
     if (defined $required_fields_annot{"accession"}) {
        $required_fields_annot{"accession"} =~ s/^\s+//; # remove blank space from start of accession
        $required_fields_annot{"accession"} =~ s/\s+$//; # remove blank space from end of accession
     }


     # get subset of these where all required fields have been defined
     # A UniProt accession and PubMed ID are also required
     if ( defined $required_fields_annot{"phi_base_acc"}
          and defined $required_fields_annot{"db_type"}
          and defined $required_fields_annot{"accession"}
          and defined $required_fields_annot{"gene_name"}
          and defined $required_fields_annot{"patho_tax"}
          and defined $required_fields_annot{"host_tax"}
          and defined $required_fields_annot{"literature_id"}
          and defined $required_fields_annot{"entered_by"}
          and $required_fields_annot{"phi_base_acc"} ne ""
          # UniProt entry must be either 6 or 10 alphanumeric characters
          and ($required_fields_annot{"accession"} =~ /^[a-zA-Z\d]{6}$/
                or $required_fields_annot{"accession"} =~ /^[a-zA-Z\d]{10}$/)
          and $required_fields_annot{"gene_name"} ne ""
          and $required_fields_annot{"host_tax"} =~ /^\d+$/  # taxon ID must be an integer
          and $required_fields_annot{"literature_id"} ne ""
          and $required_fields_annot{"literature_id"} ne "no data found"
          and $required_fields_annot{"entered_by"} ne ""
          and lc $required_fields_annot{"db_type"} eq "uniprot"
          #and lc $required_fields_annot{"literature_source"} eq "pubmed" # now we accept interactions from other sources
          and $required_fields_annot{"patho_tax"} =~ /^\d+$/  # taxon ID must be an integer
          #and $required_fields_annot{"patho_tax"} == 5518  # taxon ID for Fusarium gram
          #and $required_fields_annot{"patho_tax"} == 148305  # taxon ID for Magnaporthe oryzae
          #and $required_fields_annot{"patho_tax"} == 1307  # taxon ID for Streptococcus suis (causes human disease)
        ) {

        # get the record ID for the current annotation
        my $annot_record_id = $phi_base_annotation{"record_id"};
        #print "Record ID:$annot_record_id\tPHI-base Acc:$required_fields_annot{'phi_base_acc'}\n";

        # add the required fields of the current annotation to the required criteria annotations hash
        $required_criteria_annotations{$annot_record_id} = {%phi_base_annotation};

	#print "PHI-base Accession:$required_fields_annot{'phi_base_acc'}\n";
	#print "UniProt Accession:$required_fields_annot{'accession'}\n";
	#print "Gene Name: $required_fields_annot{'gene_name'}\n";
	#print "Pathogen Species NCBI Taxon ID:$required_fields_annot{'patho_tax'}\n";
	#print "Host Species NCBI Taxon ID:$required_fields_annot{'host_tax'}\n";
	#print "PubMed ID:$required_fields_annot{'literature_id'}\n\n";

        # if a pathogen strain taxon id is defined, then this should be used for the pathogen,
        # since it represents a more precise taxonomy. Otherwise, use the pathogen species taxon ID
        my $path_taxon_id;
        if (defined $phi_base_annotation{"patho_strain_tax"} and $phi_base_annotation{"patho_strain_tax"} ne "" and $phi_base_annotation{"patho_strain_tax"} ne "na" and $phi_base_annotation{"patho_strain_tax"} ne "no data found") {
           $path_taxon_id = $phi_base_annotation{"patho_strain_tax"};
        } else {
	   $path_taxon_id = $required_fields_annot{"patho_tax"};
        }

	# insert data into the pathogen_gene table,
        # if it does not exist already (based on combination of taxon id, gene name, and uniprot_accession)
        # if Gene Locus Identifier has been supplied, then this should also be inserted
	my $sql_statement2; 
        if (defined $phi_base_annotation{"gene_locus_id"} and $phi_base_annotation{"gene_locus_id"} ne "") {
	   $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_species_taxon_id, ncbi_taxon_id, gene_name, uniprot_accession, gene_locus_id, gene_locus_id_type) 
				SELECT $phi_base_annotation{"patho_tax"}, $path_taxon_id, '$required_fields_annot{"gene_name"}', '$required_fields_annot{"accession"}', '$phi_base_annotation{"gene_locus_id"}', '$phi_base_annotation{"gene_locus_id_type"}'
				WHERE NOT EXISTS (
				  SELECT 1 FROM pathogen_gene
				  WHERE ncbi_taxon_id = $path_taxon_id
				  AND gene_name = '$required_fields_annot{"gene_name"}'
				  AND uniprot_accession = '$required_fields_annot{"accession"}'
			       ));
        } else { # gene locus data not supplied
	   $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_species_taxon_id, ncbi_taxon_id,gene_name,uniprot_accession) 
				SELECT $phi_base_annotation{"patho_tax"}, $path_taxon_id,'$required_fields_annot{"gene_name"}','$required_fields_annot{"accession"}'
				WHERE NOT EXISTS (
				  SELECT 1 FROM pathogen_gene
				  WHERE ncbi_taxon_id = $path_taxon_id
				  AND gene_name = '$required_fields_annot{"gene_name"}'
				  AND uniprot_accession = '$required_fields_annot{"accession"}'
			       ));
        }

	my $sql_result2 = $db_conn->prepare($sql_statement2);
	$sql_result2->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene record
        my $sql_statement4 = qq(SELECT id FROM pathogen_gene
                                WHERE ncbi_taxon_id = $path_taxon_id
                                AND gene_name = '$required_fields_annot{"gene_name"}'
                                AND uniprot_accession = '$required_fields_annot{"accession"}'
                               );

	my $sql_result4 = $db_conn->prepare($sql_statement4);
	$sql_result4->execute() or die $DBI::errstr;
	my @row4 = $sql_result4->fetchrow_array();
	my $pathogen_gene_id = shift @row4;

        # insert data into pathogen_gene_allele table, including foreign key to pathogen_gene table
        # NEED TO EXTENT TO INCLUDE ALLELE-SPECIFIC DATA 
	my $sql_statement3 = qq(INSERT INTO pathogen_gene_allele (pathogen_gene_id) 
			         SELECT $pathogen_gene_id
                                 WHERE NOT EXISTS (
                                   SELECT 1 FROM pathogen_gene_allele
			           WHERE pathogen_gene_id = $pathogen_gene_id
                                 )
                               );

	my $sql_result3 = $db_conn->prepare($sql_statement3);
	$sql_result3->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene_allele record
        my $sql_statement5 = qq(SELECT id FROM pathogen_gene_allele
			        WHERE pathogen_gene_id = $pathogen_gene_id
                               );

	my $sql_result5 = $db_conn->prepare($sql_statement5);
	$sql_result5->execute() or die $DBI::errstr;
	my @row5 = $sql_result5->fetchrow_array();
	my $pathogen_gene_allele_id = shift @row5;

        # before inserting a new interaction, we need to find out if the current PHI-base accession
        # should be part of an existing interaction (i.e. in a multiple mutation)

        # mutliple mutation flag
        my $multiple_mutation = 0;
        # OLD PHI-base accession number of multiple mutation partner - MAY NEED TO CONVERT TO ARRAY TO COPE WITH MULTIPLE PARTNERS
        my $first_multi_mut_old_phi_acc_num;


        my @orig_multi_mutants_list;
        my @orig_multi_mutants_list2;
        my @record_id_list;
        my $first_mult_mut_partner_record_id;
        my $first_mult_mut_partner_old_phi_acc;
        my $first_mult_mut_partner_interaction_id;
        my $first_mult_mut_partner_new_phi_acc;
        my $first_mult_mut_partner_host_taxon_id;
        my $first_mult_mut_partner_host_strain;
        my $first_mult_mut_partner_path_strain;

        # check if the annotation is part of a "multiple mutation" interaction
        if ($phi_base_annotation{"multiple_mutation"} ne "" and $phi_base_annotation{"multiple_mutation"} ne "no" and $phi_base_annotation{"multiple_mutation"} ne "na") {

          #print $required_fields_annot{"phi_base_acc"}."\t";
          #print $required_fields_annot{"gene_name"}."\t";
          #print $required_fields_annot{"accession"}."\t";
          #print $required_fields_annot{"host_tax"}."\t";
          #print $phi_base_annotation{"multiple_mutation"}."\n";

          #print "Current Multi Mut string:$phi_base_annotation{'multiple_mutation'}\n";

 	  # get the PHI-base accessions for the multiple mutation partners based on semi-colon delimiter
          my @multi_mutants_list = split(";",$phi_base_annotation{"multiple_mutation"});
          @orig_multi_mutants_list2 = @multi_mutants_list;
          @orig_multi_mutants_list = @multi_mutants_list;

	  foreach my $orig_mm (@orig_multi_mutants_list) {
	    $orig_mm =~ s/^\s+//; # remove blank space from start of accession
	    $orig_mm =~ s/\s+$//; # remove blank space from end of accession
	  }
	  my @sorted_orig_multi_mutants_list = sort @orig_multi_mutants_list;


          foreach my $multi_mut_accession (@sorted_orig_multi_mutants_list) {

            $multi_mut_accession =~ s/^\s+//; # remove blank space from start of accession
            $multi_mut_accession =~ s/\s+$//; # remove blank space from end of accession
            $multi_mut_accession  =~ s/PHI://;

            # find the record IDs for all of the interactions with the given old PHI-base accession
            @record_id_list = grep { $required_criteria_annotations{$_}{'phi_base_acc'} eq 'PHI:'.$multi_mut_accession } keys %required_criteria_annotations;

            #print "Current RECORD ID:$annot_record_id\n";

            foreach my $record_id (@record_id_list) {

              my @local_multi_mutants_list = @orig_multi_mutants_list2;

              #print "Partner RECORD ID:$record_id\n";
              #print "Current Multi Mut string:$phi_base_annotation{'multiple_mutation'}\n";
              #print "Partner Multi Mut String:$required_criteria_annotations{$record_id}{'multiple_mutation'}\n";

              # get the PHI-base accessions for the multiple mutation partners based on semi-colon delimiter
              my @partner_multi_mutants_list = split(";",$required_criteria_annotations{$record_id}{'multiple_mutation'});

              push (@partner_multi_mutants_list, $required_criteria_annotations{$record_id}{'phi_base_acc'});
              foreach my $partner_mm (@partner_multi_mutants_list) {
                $partner_mm =~ s/^\s+//; # remove blank space from start of accession
                $partner_mm =~ s/\s+$//; # remove blank space from end of accession
              }
              my @sorted_partner_multi_mutants_list = sort @partner_multi_mutants_list;

              push (@local_multi_mutants_list, $phi_base_annotation{"phi_base_acc"});
              foreach my $mm (@local_multi_mutants_list) {
                $mm =~ s/^\s+//; # remove blank space from start of accession
                $mm =~ s/\s+$//; # remove blank space from end of accession
                #print "Current MM:$mm\n";
              }
              my @sorted_local_multi_mutants_list = sort @local_multi_mutants_list;

              #check if there is no difference between the arrays
              if ( !array_diff(@sorted_local_multi_mutants_list,@sorted_partner_multi_mutants_list) ) {
                 #print "Arrays the same\n";
                 $first_mult_mut_partner_record_id = $record_id;
                 $first_mult_mut_partner_old_phi_acc = $required_criteria_annotations{$record_id}{'phi_base_acc'};
                 $first_mult_mut_partner_interaction_id = $required_criteria_annotations{$record_id}{'interaction_id'};
                 $first_mult_mut_partner_new_phi_acc = $required_criteria_annotations{$record_id}{'new_phibase_acc'};
                 #print "Other annotation record ID: $first_mult_mut_partner_record_id\n";
                 #print "Other annotation old PHI Acc: $first_mult_mut_partner_old_phi_acc\n";
                 #print "Other annotation new PHI Acc: $first_mult_mut_partner_new_phi_acc\n";
                 #print "Other annotation interaction ID: $first_mult_mut_partner_interaction_id\n";


                 $first_mult_mut_partner_host_taxon_id = $required_criteria_annotations{$record_id}{'host_tax'};
                 $first_mult_mut_partner_host_strain = $required_criteria_annotations{$record_id}{'host_strain_tax'};
                 $first_mult_mut_partner_path_strain = $required_criteria_annotations{$record_id}{'strain_name'};
                 #print "Current Host tax ID:$required_fields_annot{'host_tax'}\n";
                 #print "Partner Host tax ID:$first_mult_mut_partner_host_taxon_id\n";
                 #print "Current Host strain:$phi_base_annotation{'host_strain_tax'}\n";
                 #print "Partner Host strain:$first_mult_mut_partner_host_strain\n";
                 #print "Current Path strain:$phi_base_annotation{'strain_name'}\n";
                 #print "Partner Path strain:$first_mult_mut_partner_path_strain\n";

                 if ( $required_fields_annot{"host_tax"} == $first_mult_mut_partner_host_taxon_id
                      and $phi_base_annotation{"host_strain_tax"} eq $first_mult_mut_partner_host_strain
                      and $phi_base_annotation{"strain_name"} eq $first_mult_mut_partner_path_strain ) {
                    #print "HOSTS & PATHS MATCH\n";
                    $multiple_mutation = 1;

                    # TODO: NEED TO ADD DETAILS OF THE ALLELE EXPRESSION, WHERE KNOWN
		    my $inner_sql_statement = qq(
						 INSERT INTO interaction_pathogen_gene_allele (interaction_id,pathogen_gene_allele_id) 
						   VALUES ($first_mult_mut_partner_interaction_id,$pathogen_gene_allele_id);
						 INSERT INTO obsolete (phi_base_accession,obsolete_accession)
						   VALUES ('$first_mult_mut_partner_new_phi_acc','$required_fields_annot{"phi_base_acc"}')
						);
		    my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

		    #print "Multiple mutation interaction_pathogen_gene_allele record inserted successfully\n";

		    # TODO: NEED TO ADDITIONALLY FIND OUT IF HOST TAXONOMY IDs MATCH (AND POSSIBLY PUBMED IDs)
		    last;

                 } else {
                    #print "HOSTS & PATHS DO NOT MATCH\n";
                 }

              } else {
                 #print "Arrays different\n";
              }

            } # end foreach record ID of the multiple mutant partner

            # if the multiple mutation partner has already been found,
            # then no need to examine other partners
            if ($multiple_mutation) {
              last;
            }

          } # end foreach multiple mutant partner

        } # end if multiple mutation


        if (not $multiple_mutation) {

          # increment interaction counter, to become new PHI-base accession
          $interaction_num++;
          my $phi_base_accession = "PHI:I".$interaction_num;

          # insert a new record into the interaction table, returning the interaction identifier
	  my $sql_statement1 = qq(INSERT INTO interaction (phi_base_accession,curation_date) 
	                            VALUES ('$phi_base_accession',current_date) RETURNING id;);
	  my $sql_result1 = $db_conn->prepare($sql_statement1);
	  $sql_result1->execute() or die $DBI::errstr;
	  my @row1 = $sql_result1->fetchrow_array();
	  my $interaction_id = shift @row1;

	  # add the interaction ID and new PHI-base accession to the required_criteria_annotations array
	  $required_criteria_annotations{$annot_record_id}{'interaction_id'} = $interaction_id;
	  $required_criteria_annotations{$annot_record_id}{'new_phibase_acc'} = $phi_base_accession;

          # insert record to reference back to old phibase accession
          my $sql_statement6 = qq(INSERT INTO obsolete (phi_base_accession,obsolete_accession)
                                    VALUES ('$phi_base_accession','$required_fields_annot{"phi_base_acc"}');
                                 );
	  my $sql_result6 = $db_conn->prepare($sql_statement6);
	  $sql_result6->execute() or die $DBI::errstr;

          # declare an array for the list of pubmed ids and digital object identifiers (DOIs)
          my @pubmed_id_list;
          my @doi_list;

          if ( $interaction_id and $pathogen_gene_allele_id ) {
             # add records for the literature and pathogen gene allele tables associated with the interaction,
             # using the interaction id as a foreign key to the interaction table
             # TODO: NEED TO ADD DETAILS OF THE ALLELE EXPRESSION, WHERE KNOWN
	     my $inner_sql_statement = qq(
				          INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
				            VALUES ($interaction_id,$required_fields_annot{"host_tax"});
				          INSERT INTO interaction_pathogen_gene_allele (interaction_id,pathogen_gene_allele_id) 
					    VALUES ($interaction_id,'$pathogen_gene_allele_id');
				         );
	     my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

             # Add the PubMed literature to the interaction_literature table,
             # inserting a separate record for each PubMed ID, which are separated by semi-colon

             if (lc $required_fields_annot{"literature_source"} eq "pubmed") {
		# separate list of PubMed IDs, based on semi-colon delimiter
		@pubmed_id_list = split (";",$phi_base_annotation{"literature_id"});
             }

             if (defined $phi_base_annotation{"doi"} and $phi_base_annotation{"doi"} ne "" and $phi_base_annotation{"doi"} ne "no data found") {
		# separate list of DOIs, based on semi-colon delimiter
		@doi_list = split (";",$phi_base_annotation{"doi"});
             }

             if (@pubmed_id_list) {
		# insert interaction_lite:serature record for each PubMed ID
		foreach my $pubmed_id (@pubmed_id_list) {

		   $pubmed_id =~ s/^\s+//; # remove blank space from start of PubMed ID
		   $pubmed_id =~ s/\s+$//; # remove blank space from end of PubMed ID

		   # insert the PubMed ID into interaction_literature table,
		   # using the interaction id and as a foreign key
                   # and include a DOI, if available
                   if ( my $doi = pop(@doi_list) ) {
		      $inner_sql_statement = qq(INSERT INTO interaction_literature (interaction_id, pubmed_id, doi)
					        VALUES ($interaction_id, '$pubmed_id', '$doi');
					       );
                   } else {
		      $inner_sql_statement = qq(INSERT INTO interaction_literature (interaction_id, pubmed_id)
					        VALUES ($interaction_id, '$pubmed_id');
					       );
                   }
		   $inner_sql_result = $db_conn->prepare($inner_sql_statement);
		   #$inner_sql_result->execute() or die $DBI::errstr;
		   $inner_sql_result->execute(); # ignore errors due to possible duplication in spreadsheet

		}
             } elsif (@doi_list) { # DOIs available, but no PubMed IDs
print "DOI provided, but no PubMed ID for $required_fields_annot{'phi_base_acc'}\n";
		foreach my $doi (@doi_list) {

		   $doi =~ s/^\s+//; # remove blank space from start of PubMed ID
		   $doi =~ s/\s+$//; # remove blank space from end of PubMed ID

		   $inner_sql_statement = qq(INSERT INTO interaction_literature (interaction_id, doi)
		   	                     VALUES ($interaction_id, '$doi');
					    );
		   $inner_sql_result = $db_conn->prepare($inner_sql_statement);
		   #$inner_sql_result->execute() or die $DBI::errstr;
		   $inner_sql_result->execute(); # ignore errors due to possible duplication in spreadsheet

		}
             } else {
#		$invalid_literature_count++;
		print STDERR "PHI-base ERROR: No PubMed ID or DOI provided for $required_fields_annot{'phi_base_acc'}\n";
#		print INVALID_LIT_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\n";
             }

	  }


          # for the curators, need to split list of curators base on semi-colon delimiter
          my $curators_string = $required_fields_annot{"entered_by"};

	  # get the unique identifier for the curator,
          # based on initials given in the entered_by field
          my @curators = split(";",$curators_string);

          # for each curator, need to get the curator identifier,
          # then insert an interaction_curator record
          foreach my $curator_init (@curators) {

            $curator_count++;

            $curator_init =~ s/^\s+//; # remove blank space from start of curator initials
            $curator_init =~ s/\s+$//; # remove blank space from end of curator initials

            my $sql_statement = qq(SELECT id FROM curator
                                     WHERE initials = '$curator_init';
                                  );

	    my $sql_result = $db_conn->prepare($sql_statement);
	    $sql_result->execute() or die $DBI::errstr;
	    my @row = $sql_result->fetchrow_array();
	    my $curator_id = shift @row;

            # insert data into interaction_curator table,
            # with foreign keys to the interaction table and the curator table 
	    $sql_statement = qq(INSERT INTO interaction_curator (interaction_id, curator_id)
                                  VALUES ($interaction_id, $curator_id);
                               );
	    $sql_result = $db_conn->prepare($sql_statement);
	    #$sql_result->execute() or die $DBI::errstr;
	    $sql_result->execute(); # error handling removed due to initials not available in spreadsheet

          } # end foreach curator

          # if a valid species expert exists,
          # need to retrieve the curator identifier, based on the initials
          # then insert the appropriate species expert record (if it does not already exist)
          # based on both the curator id and the pathogen taxon id
          my $species_expert_init = $phi_base_annotation{"species_expert"};

          $species_expert_init =~ s/^\s+//; # remove blank space from start of species expert initials
          $species_expert_init =~ s/\s+$//; # remove blank space from end of species expert initials

          my $sql_statement = qq(SELECT id FROM curator
                                   WHERE initials = '$species_expert_init';
                                );

	  my $sql_result = $db_conn->prepare($sql_statement);
	  $sql_result->execute() or die $DBI::errstr;
	  my @row = $sql_result->fetchrow_array();
	  my $curator_id = shift @row;

          # insert data into species_expert table (if it does not already exist),
          # using the pathogen taxon id and foreign key to the curator table 
	  if (defined $curator_id and $curator_id ne "" and $curator_id ne "na") {
             $species_expert_count++;
  	     $sql_statement = qq(INSERT INTO species_expert (ncbi_taxon_id, curator_id)
                                 SELECT $path_taxon_id, $curator_id
                                 WHERE NOT EXISTS (
                                     SELECT 1 FROM species_expert
                                     WHERE ncbi_taxon_id = $path_taxon_id
                                     AND curator_id = $curator_id
                                   )
                                );
	     $sql_result = $db_conn->prepare($sql_statement);
	     $sql_result->execute() or die $DBI::errstr;
          }


          # create hash of key/value pairs for attribute/value of all the pathogen_phenotypes
          my %pathogen_phenotypes = (
                           'Mating Defect'     => $phi_base_annotation{"mating_defect"},
                           'Pre-penetration'   => $phi_base_annotation{"prepenetration"},
                           'Penetration'       => $phi_base_annotation{"penetration"},
                           'Post-penetration'  => $phi_base_annotation{"postprepenetration"},
                           'Vegetative Spores' => $phi_base_annotation{"vegetative_spores"},
                           'Sexual Spores'     => $phi_base_annotation{"sexual_spores"},
                           'In Vitro Growth'   => $phi_base_annotation{"in_vitro_growth"},
                           'Spore Germination' => $phi_base_annotation{"spore_germination"}
                        );


          # for each of the pathogen_phenotypes, retrieve the id for the relevant pathogen_phenotype,
          # then retrieve the id for each pathogen_phenotype value (if available)
          # then insert a interaction_phi_pathogen_phenotype record,
          # based on combination of interaction id and phenotype id
          foreach my $pathogen_phenotype_attribute (keys %pathogen_phenotypes)
          {

             my $pathogen_phenotype_values_string = $pathogen_phenotypes{$pathogen_phenotype_attribute};

             if (defined $pathogen_phenotype_values_string and $pathogen_phenotype_values_string ne "" and lc($pathogen_phenotype_values_string) ne "none" and lc($pathogen_phenotype_values_string) ne "na" and lc($pathogen_phenotype_values_string) ne "nd" and lc($pathogen_phenotype_values_string) ne "no data found" and lc($pathogen_phenotype_values_string) ne "no data") {

                # separate list of pathogen_phenotype values, based on semi-colon delimiter
                my @pathogen_phenotype_values = split (";",$pathogen_phenotype_values_string);

                # insert interaction_phi_pathogen_phenotype record for each pathogen_phenotype value
                foreach my $pathogen_phenotype_value (@pathogen_phenotype_values) {

                   $pathogen_phenotype_value =~ s/^\s+//; # remove blank space from start of pathogen_phenotype value
                   $pathogen_phenotype_value =~ s/\s+$//; # remove blank space from end of pathogen_phenotype value

                   my @phi_phenotype_id_list;

		   # insert data into interaction_phi_pathogen_phenotype table,
		   # using the interaction id and PHI phenotype id as foreign keys 
		   if (defined $pathogen_phenotype_value) {

                      if ($pathogen_phenotype_attribute eq "Mating Defect" and lc($pathogen_phenotype_value) eq "yes") {
                        push(@phi_phenotype_id_list,'PHIPO:0000015');
                        push(@phi_phenotype_id_list,'PHIPO:0000098');
                      } 
                      elsif ($pathogen_phenotype_attribute eq "Pre-penetration" and lc($pathogen_phenotype_value) eq "yes") {
                        push(@phi_phenotype_id_list,'PHIPO:0000098');
                      }
                      elsif ($pathogen_phenotype_attribute eq "Penetration" and lc($pathogen_phenotype_value) eq "yes") {
                        push(@phi_phenotype_id_list,'PHIPO:0000099');
                      }
                      elsif ($pathogen_phenotype_attribute eq "Post-penetration" and lc($pathogen_phenotype_value) eq "yes") {
                        push(@phi_phenotype_id_list,'PHIPO:0000097');
                      }
                      elsif ($pathogen_phenotype_attribute eq "Vegetative Spores") {

                        # note that more than one pattern can be present in the value string
                        if ( lc($pathogen_phenotype_value) =~ /aberrant/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000021');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /abnormal/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000021');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /altered/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000018');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /decreased/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000020');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /defective/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000021');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /deficient/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000020');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /enhanced/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000019');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /increased/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000019');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /no conidia/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000091');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /no sporulation/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000100');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /none/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000091');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /reduced/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000020');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /unaffected/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000034');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /wild type/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000034');
                        }
                      
                      }
                      elsif ($pathogen_phenotype_attribute eq "Sexual Spores") {

                        # note that more than one pattern can be present in the value string
                        if ( lc($pathogen_phenotype_value) =~ /aberrant/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000026');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /defective/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000026');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /enhanced/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000024');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /few/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000025');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /increased/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000024');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /none/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000025');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /reduced/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000025');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /sterile/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000025');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /unaffected/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000035');
                        }
                        if ( lc($pathogen_phenotype_value) =~ /wild type/ ) {
                          push(@phi_phenotype_id_list,'PHIPO:0000035');
                        }

                      }
                      elsif ($pathogen_phenotype_attribute eq "In Vitro Growth") {

                        # NOT CURRENTLY IMPLEMENTED
                        # MAY NEED TO EXPAND THE PHI PHENOTYPE ONTOLOGY TO INCLUDE THE IN VITRO GROWTH TERMS

                      }
                      elsif ($pathogen_phenotype_attribute eq "Spore Germination") {
                       
                        # NOT CURRENTLY IMPLEMENTED
                        # AMBIGUITY IF SPORE GERMINATION REFERS TO SEXUAL SPORES OR VEGETATIVE SPORES
                        
                      }
                       
                   } # end if pathogen phenotype defined

                   if (@phi_phenotype_id_list) {

                      foreach my $phi_phenotype_id (@phi_phenotype_id_list) {

#                        print "PATH_PHENOTYPE ATT:$pathogen_phenotype_attribute\n";
#                        print "PATH_PHENOTYPE VALUE:$pathogen_phenotype_value\n";
#                        print "PHI PHENOTYPE ID:$phi_phenotype_id\n";
#                        print "INTERACTION ID:$interaction_id\n";
#                        print "PHI-BASE ACC:$required_fields_annot{'phi_base_acc'}\n";
                         $pathogen_phenotype_count++;
                         
                         # since some spreadsheet columns, such as "Mating defect prior to penetration" and "Pre-penetration defect",
                         # will match to the same PHI phenotype identifier, the SQL insert should check if the relevant values already exist
		         $sql_statement = qq(INSERT INTO interaction_phi_pathogen_phenotype (interaction_id, phi_phenotype_id)
		  			     SELECT $interaction_id, '$phi_phenotype_id'
                                             WHERE NOT EXISTS (
                                               SELECT 1 FROM interaction_phi_pathogen_phenotype
                                               WHERE interaction_id = $interaction_id
                                               AND phi_phenotype_id = '$phi_phenotype_id'
					    ));
		         $sql_result = $db_conn->prepare($sql_statement);
		         $sql_result->execute() or die $DBI::errstr;
                         print PATH_PHENOTYPE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$pathogen_phenotype_attribute\t$pathogen_phenotype_value\n";
                      }

                   } else {
                      $invalid_pathogen_phenotype_count++;
                      print INVALID_PATH_PHENOTYPE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$pathogen_phenotype_attribute\t$pathogen_phenotype_value\n";
                   }

                } # end foreach pathogen phenotype value

             } # end if pathogen phenotype values

          } # end foreach pathogen phenotype attribute

          
          # Get the GO annotations
          my $go_annot_string = $phi_base_annotation{"go_annotation"};

          $go_annot_string =~ s/^\s+//; # remove blank space from start of string
          $go_annot_string =~ s/\s+$//; # remove blank space from end of string

          # if the GO annotations string is empty, no GO annotations have been supplied
          if (defined $go_annot_string and $go_annot_string ne "" and $go_annot_string ne "data not found" and $go_annot_string ne "no data found") {

            # need to split list based on semi-colon delimiter
            my @go_annotations = split(";",$go_annot_string);

            # for each GO annotation, need to get the GO ID and the evidence code,
            # then insert an interaction_go_annotation record for this GO annotation
            foreach my $go_annotation (@go_annotations) {

               $go_annotation_count++;

               # GO annotations should have both a term identifier
               # and an evidence code, which are separated by a comma
               my @go_parts = split(",",$go_annotation);

               # first part will be the GO ID
               my $go_id = shift @go_parts;
               $go_id =~ s/^\s+//; # remove blank space from start of GO ID
               $go_id =~ s/\s+$//; # remove blank space from end of GO ID

               # second part will be the GO evidence code, if supplied
               my $go_evid_code = shift @go_parts;
               if (defined $go_evid_code) {
                 $go_evid_code =~ s/^\s+//; # remove blank space from start of GO evidence code
                 $go_evid_code =~ s/\s+$//; # remove blank space from end of GO evidence code
               }

               if (defined $go_evid_code and $go_evid_code ne "") {

                  # insert data into interaction_go_annotation table,
                  # with foreign keys to the interaction table and the go_evidence_code_table 
#	          $sql_statement = qq(INSERT INTO interaction_go_annotation (interaction_id, go_id, go_evidence_code)
#                                        VALUES ($interaction_id, '$go_id', '$go_evid_code');
#                                     );

                  # insert data into pathogen_gene_go_annotation table,
                  # with foreign keys to the pathogen_gene table, PubMed, and the go_evidence_code_table
                  # since there may be more than one PubMed ID, need to insert a record for each PubMed ID
                  foreach my $pubmed_id (@pubmed_id_list) {

                     $pubmed_id =~ s/^\s+//; # remove blank space from start of PubMed ID
                     $pubmed_id =~ s/\s+$//; # remove blank space from end of PubMed ID

                     # since we only want one GO annotation per gene, rather than per interaction,
                     # it is necessary to check if the current GO annotation has already been inserted
	             $sql_statement = qq(INSERT INTO pathogen_gene_go_annotation (pathogen_gene_id, pubmed_id, go_id, go_evidence_code)
                                         SELECT $pathogen_gene_id, '$pubmed_id', '$go_id', '$go_evid_code'
                                         WHERE NOT EXISTS (
                                           SELECT 1 FROM pathogen_gene_go_annotation
                                           WHERE pathogen_gene_id = $pathogen_gene_id
                                           AND pubmed_id = '$pubmed_id'
                                           AND go_id = '$go_id'
                                        ));
	             $sql_result = $db_conn->prepare($sql_statement);

		     # execute SQL insert and, if successful, print to file
		     if ( $sql_result->execute() ) {
		        $go_with_evid_count++;
		        print GO_WITH_EVID_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$go_id\t$go_evid_code\n";
		     } else {  # SQL insert unsuccessful, then log error
		       $invalid_go_count++;
	               print STDERR "\nPHI-base ERROR: Evidence code $go_evid_code is not valid for $required_fields_annot{'phi_base_acc'}, $go_id\n\n";
		       print INVALID_GO_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$go_id\t$go_evid_code\n";
		     }

		  } # end foreach PubMed ID

               } else { # GO evidence code not supplied

                  # insert data into pathogen_gene_go_annotation table,
                  # with foreign keys to the pathogen_gene table, PubMed, and the go_evidence_code_table
                  # since there may be more than one PubMed ID, need to insert a record for each PubMed ID
                  foreach my $pubmed_id (@pubmed_id_list) {

                     $pubmed_id =~ s/^\s+//; # remove blank space from start of PubMed ID
                     $pubmed_id =~ s/\s+$//; # remove blank space from end of PubMed ID

                     # since we only want one GO annotation per gene, rather than per interaction,
                     # it is necessary to check if the current GO annotation has already been inserted
	             $sql_statement = qq(INSERT INTO pathogen_gene_go_annotation (pathogen_gene_id, pubmed_id, go_id)
                                         SELECT $pathogen_gene_id, '$pubmed_id', '$go_id'
                                         WHERE NOT EXISTS (
                                           SELECT 1 FROM pathogen_gene_go_annotation
                                           WHERE pathogen_gene_id = $pathogen_gene_id
                                           AND pubmed_id = '$pubmed_id'
                                           AND go_id = '$go_id'
                                        ));
	             $sql_result = $db_conn->prepare($sql_statement);
	             $sql_result->execute() or die $DBI::errstr;

                     $go_without_evid_count++;
                     print GO_WITHOUT_EVID_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$go_id\n";
    
		  } # end foreach PubMed ID

=pod
                  # insert data into interaction_go_annotation table,
                  # with foreign key to the interaction table, but without GO evidence code 
	          $sql_statement = qq(INSERT INTO interaction_go_annotation (interaction_id, go_id)
                                        VALUES ($interaction_id, '$go_id');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);
	          #$sql_result->execute() or die $DBI::errstr;
	          $sql_result->execute(); # ignore error due to possible duplications in spreadsheet

                  $go_without_evid_count++;
                  print GO_WITHOUT_EVID_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$go_id\n";
=cut
    
               } # end else GO evidence code
              
             } # end foreach GO term

          } # end if GO terms supplied 

          
          # anti-infective chemicals are given as CAS Registry IDs
          my $cas_string = $phi_base_annotation{"cas"};

          # if the CAS string is empty, no anti-infectives have been supplied
          if (defined $cas_string and $cas_string ne "" and $cas_string ne "no data found") {

            # need to split list based on semi-colon delimiter
            my @cas_entries = split(";",$cas_string);

            # for each anti-infective, need to get the CAS Registry ID, then insert
            # a chemical record (if it does not already exist) and 
            # an interaction_anti-infective_chemical record for this interaction
            foreach my $cas_entry (@cas_entries) {

               $anti_infective_count++;

               # entries sometimes have a name before the ID, separated by a colon,
               # other times no name is given before the CAS Registry ID,
               # so need to extract the ID, based on colon delimiter
               # (which will always be the last part)
               my @cas_parts = split(":",$cas_entry);
               my $cas_id = pop(@cas_parts);

               $cas_id =~ s/^\s+//; # remove blank space from start of CAS ID
               $cas_id =~ s/\s+$//; # remove blank space from end of CAS ID

   	       # insert data into the chemical_table,
               # if it does not exist already (based on CAS ID)
  	       my $sql_statement = qq(INSERT INTO chemical (cas_registry) 
	 		                 SELECT '$cas_id'
                                       WHERE NOT EXISTS (
                                         SELECT 1 FROM chemical
                                         WHERE cas_registry = '$cas_id'
                                       )
                                     );

	       my $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;

	       # get the unique identifier for the inserted chemical
               $sql_statement = qq(SELECT id FROM chemical
                                     WHERE cas_registry = '$cas_id';
                                  );
	       $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;
	       @row = $sql_result->fetchrow_array();
	       my $chemical_id = shift @row;

               # insert data into interaction_anti_infective_chemical table,
               # with foreign keys to the interaction table and the chemical table 
	       $sql_statement = qq(INSERT INTO interaction_anti_infective_chemical (interaction_id, chemical_id)
                                     VALUES ($interaction_id, $chemical_id);
                                  );
	       $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;
              
             } # end foreach anti-infective chemical
 
          } else {  # no anti-infective supplied 
             $no_anti_infective_count++;
          }

          
          # get the inducer chemical names
          my $inducer_string = $phi_base_annotation{"inducer"};

          # if the inducer string is empty, no inducers have been supplied
          if (defined $inducer_string and $inducer_string ne "" and $inducer_string ne "na" and $inducer_string ne "no data found") {

            # need to split list based on semi-colon delimiter
            my @inducer_names = split(";",$inducer_string);

            # for each inducer, need to get the chemical identifier from the
            # chemical table record, then insert 
            # an interaction_inducer_chemical record for this interaction
            foreach my $inducer_name (@inducer_names) {

               $inducer_count++;

               $inducer_name =~ s/^\s+//; # remove blank space from start of string
               $inducer_name =~ s/\s+$//; # remove blank space from end of string

               # to permit single-quote in inducer name
               # need to escape it by adding a second single quote
               $inducer_name =~ s/'/''/g; 

   	       # insert data into the chemical_table,
               # if it does not exist already (based on the inducer name)
  	       my $sql_statement = qq(INSERT INTO chemical (name) 
	 		                 SELECT lower('$inducer_name')
                                       WHERE NOT EXISTS (
                                         SELECT 1 FROM chemical
                                         WHERE name = lower('$inducer_name')
                                       )
                                     );

	       my $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;

	       # get the unique identifier for the chemical
               $sql_statement = qq(SELECT id FROM chemical
                                     WHERE name = lower('$inducer_name');
                                  );
	       $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;
	       @row = $sql_result->fetchrow_array();
	       my $chemical_id = shift @row;

               # insert data into interaction_inducer_chemical table,
               # with foreign keys to the interaction table and the chemical table 
	       $sql_statement = qq(INSERT INTO interaction_inducer_chemical (interaction_id, chemical_id)
                                     VALUES ($interaction_id, $chemical_id);
                                  );
	       $sql_result = $db_conn->prepare($sql_statement);
	       $sql_result->execute() or die $DBI::errstr;
              
             } # end foreach inducer chemical

          } else {  # no inducers supplied 
             $no_inducer_count++;
          }

          
          # get the experimental evidence
          my $exp_evid_string = $phi_base_annotation{"experimental_evidence"};

          # if the experimental evidence string is empty, no evidence have been supplied
          if (defined $exp_evid_string and $exp_evid_string ne "" and $exp_evid_string ne "no data found") {

             $exp_spec_count++;

             $exp_evid_string =~ s/^\s+//; # remove blank space from start of string
             $exp_evid_string =~ s/\s+$//; # remove blank space from end of string

             # using the experiment specification mappings,
             # get the appropriate ontology identifiers associated with the exp evidence
             my $exp_spec_id_list = $exp_spec_mapping{$exp_evid_string};

             # if identifiers are present, then insert the appropriate
             # records into the interaction_experiment_spec table
             if ($exp_spec_id_list) {

                $exp_spec_id_list =~ s/^\s+//; # remove blank space from start of string
                $exp_spec_id_list =~ s/\s+$//; # remove blank space from end of string

                # need to split list based on semi-colon delimiter
                my @exp_spec_ontology_ids = split(";",$exp_spec_id_list);

                # for each exp spec id, insert data into interaction_experiment_spec table,
                # with foreign keys to the interaction_host table and the experiment spec ontology
                foreach my $exp_spec_id (@exp_spec_ontology_ids) {
                  $exp_spec_term_count++;
	          $sql_statement = qq(INSERT INTO interaction_experiment_spec (interaction_id, experiment_spec_id)
                                        VALUES ($interaction_id, '$exp_spec_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print EXP_SPEC_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$exp_evid_string\t$exp_spec_id\n";
                }

             } else { # no experiment spec identifiers
                $invalid_exp_spec_count++;
                print STDERR "ERROR:Experiment evidence $exp_evid_string given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_EXP_SPEC_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$exp_evid_string\n";
             }

          } # end if experiment evidence supplied 
         
 
          # get the host response
          my $host_response = $phi_base_annotation{"host_response"};

          # if the host response string is empty, no host response has been supplied
          if (defined $host_response and $host_response ne "" and $host_response ne "no data found") {

             $host_response_count++;

             $host_response =~ s/^\s+//; # remove blank space from start of string
             $host_response =~ s/\s+$//; # remove blank space from end of string

             # using the host response mappings,
             # get the appropriate ontology identifiers associated with the host response
             my $phi_phenotype_id_list = $host_response_mapping{$host_response};

             # if identifiers are present, then insert the appropriate
             # records into the interaction_phi_host_phenotype table
             if ($phi_phenotype_id_list) {

                my $sql_statement = qq(SELECT id FROM interaction_host
                                         WHERE interaction_id = $interaction_id
                                         AND ncbi_taxon_id = $required_fields_annot{"host_tax"}
                                      );

	        my $sql_result = $db_conn->prepare($sql_statement);
	        $sql_result->execute() or die $DBI::errstr;
	        my @row = $sql_result->fetchrow_array();
	        my $interaction_host_id = shift @row;

                $phi_phenotype_id_list =~ s/^\s+//; # remove blank space from start of string
                $phi_phenotype_id_list =~ s/\s+$//; # remove blank space from end of string

                # need to split list based on semi-colon delimiter
                my @phi_phenotype_ids = split(";",$phi_phenotype_id_list);

                # for each phi phenotype id, insert data into interaction_phi_host_phenotype table,
                # with foreign keys to the interaction_host table and the PHI phenotype ontology
                foreach my $phi_phenotype_id (@phi_phenotype_ids) {
                  $host_response_term_count++;
	          $sql_statement = qq(INSERT INTO interaction_phi_host_phenotype (interaction_host_id, phi_phenotype_id)
                                        VALUES ($interaction_host_id, '$phi_phenotype_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print HOST_RES_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$host_response\t$phi_phenotype_id\n";
                }

             } else { # no host response identifiers
                $invalid_host_response_count++;
                print STDERR "ERROR:Host response $host_response given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_HOST_RES_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$host_response\n";
             }

          } # end if host response supplied 
          
          
          # get the phenotype outcome
          my $phenotype_outcome_string = $phi_base_annotation{"phenotype"};

          # if the phenotype outcome string is empty, no phenotype has been supplied
          if (defined $phenotype_outcome_string and $phenotype_outcome_string ne "" and $phenotype_outcome_string ne "no data found") {

             $phenotype_outcome_count++;

             $phenotype_outcome_string =~ s/^\s+//; # remove blank space from start of string
             $phenotype_outcome_string =~ s/\s+$//; # remove blank space from end of string

             # using the phenotype outcome mappings,
             # get the appropriate ontology identifier associated with the phenotype outcome
             my $phi_phenotype_id = $phenotype_outcome_mapping{$phenotype_outcome_string};

             # if identifier is present, then insert the appropriate
             # record into the interaction_phi_interaction_phenotype table
             if ($phi_phenotype_id) {

                $phi_phenotype_id =~ s/^\s+//; # remove blank space from start of string
                $phi_phenotype_id =~ s/\s+$//; # remove blank space from end of string

                # insert data into the appropriate interaction_phi_interaction_phenotype table,
                # with for a foreign key to the phi_phenotype ontology
                $phenotype_outcome_term_count++;
	        $sql_statement = qq(INSERT INTO interaction_phi_interaction_phenotype
                                      (interaction_id, phi_phenotype_id)
                                      VALUES ($interaction_id, '$phi_phenotype_id');
                                   );
	        $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                print PHEN_OUTCOME_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$phenotype_outcome_string\t$phi_phenotype_id\n";

             } else { # no phi phenotype identifier
                $invalid_phenotype_outcome_count++;
                print STDERR "ERROR:Phenotype outcome $phenotype_outcome_string given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_PHEN_OUTCOME_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$phenotype_outcome_string\n";
             }

          } # end if phenotype outcome supplied 
          
          
          # get the diseases for this interaction
          my $disease_string = $phi_base_annotation{"disease_name"};

          # if the disease string is empty, no disease has been supplied
          if (defined $disease_string and $disease_string ne "" and $disease_string ne "no data found") {

             $disease_count++;

             $disease_string =~ s/^\s+//; # remove blank space from start of string
             $disease_string =~ s/\s+$//; # remove blank space from end of string

             # using the disease mappings,
             # get the appropriate ontology identifiers associated with the disease
             my $disease_id_list = $disease_mapping{lc($disease_string)};

             # if identifiers are present, then insert the appropriate
             # records into the interaction_disease table
             if ($disease_id_list) {

                $disease_id_list =~ s/^\s+//; # remove blank space from start of string
                $disease_id_list =~ s/\s+$//; # remove blank space from end of string

                # need to split list based on semi-colon delimiter
                my @disease_ontology_ids = split(";",$disease_id_list);

                # for each disease id, insert data into interaction_disease table,
                # with foreign keys to the interaction table and either the disease ontology
                # or plant disease ontology
                foreach my $disease_id (@disease_ontology_ids) {
                  $disease_term_count++;
	          $sql_statement = qq(INSERT INTO interaction_disease (interaction_id, disease_id)
                                        VALUES ($interaction_id, '$disease_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print DISEASE_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$disease_string\t$disease_id\n";
                }

             } else { # no experiment spec identifiers
                $invalid_disease_count++;
                #print STDERR "ERROR:Disease $disease_string given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_DISEASE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$disease_string\n";
             }

          } else { # no disease supplied
              $without_disease_count++;
              print WITHOUT_DISEASE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\n";
          }
          
         
=pod
NOW FINDING DISEASE BASED ON MAPPING FILE, RATHER THAN NAME SEARCH IN ONTOLOGY FILE
          # if the disease string is empty, no disease has been supplied
          if (defined $disease_string and $disease_string ne "" and $disease_string ne "no data found") {

             $disease_count++;

             # need to split list based on semi-colon delimiter
             my @disease_names = split(";",$disease_string);

             foreach my $disease_name (@disease_names) {

               $disease_name =~ s/^\s+//; # remove blank space from start of string
               $disease_name =~ s/\s+$//; # remove blank space from end of string

               # using the disease mappings,
               # get the appropriate ontology identifiers associated with the diseases
               my $disease_id = $combined_disease_mapping{ lc($disease_name) };

               # if identifier is present, then insert the appropriate
               # record into the interaction_disease table
               if ($disease_id) {

                  # insert data into interaction_disease table,
                  # with foreign keys to the interaction table and the disease ontology
                  $disease_term_count++;
  	          #$sql_statement = qq(INSERT INTO interaction_disease (interaction_id, disease_id, disease_severity_id)
	          $sql_statement = qq(INSERT INTO interaction_disease (interaction_id, disease_id)
                                        VALUES ($interaction_id, '$disease_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print DISEASE_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$disease_id\t$disease_name\n";

               } else { # no disease identifier

                  $invalid_disease_count++;
                  #print STDERR "ERROR:Disease $disease_name given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                  print INVALID_DISEASE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$disease_name\n";

               }

             } # end foreach disease name
          
          } else { # no disease supplied

              $without_disease_count++;
              print WITHOUT_DISEASE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\n";
          }
=cut          

 
          # get the tissue
          my $tissue = $phi_base_annotation{"tissue"};

          # if the tissue string is empty, no tissue has been supplied
          if (defined $tissue and $tissue ne "" and $tissue ne "no data found") {

             $tissue_count++;

             $tissue =~ s/^\s+//; # remove blank space from start of string
             $tissue =~ s/\s+$//; # remove blank space from end of string

             # using the tissue mappings,
             # get the appropriate ontology identifiers associated with the tissue
             my $tissue_id_list = $tissue_mapping{$tissue};

             # if identifiers are present, then insert the appropriate
             # records into the interaction_host_tissue table
             if ($tissue_id_list) {

                my $sql_statement = qq(SELECT id FROM interaction_host
                                         WHERE interaction_id = $interaction_id
                                         AND ncbi_taxon_id = $required_fields_annot{"host_tax"}
                                      );

	        my $sql_result = $db_conn->prepare($sql_statement);
	        $sql_result->execute() or die $DBI::errstr;
	        my @row = $sql_result->fetchrow_array();
	        my $interaction_host_id = shift @row;

                $tissue_id_list =~ s/^\s+//; # remove blank space from start of string
                $tissue_id_list =~ s/\s+$//; # remove blank space from end of string

                # need to split list based on semi-colon delimiter
                my @tissue_ids = split(";",$tissue_id_list);

                # for each tissue id, insert data into interaction_host_tissue table,
                # with foreign keys to the interaction_host table and the BRENDA tissue ontology
                foreach my $tissue_id (@tissue_ids) {
                  $tissue_term_count++;
	          $sql_statement = qq(INSERT INTO interaction_host_tissue (interaction_host_id, brenda_tissue_id)
                                        VALUES ($interaction_host_id, '$tissue_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print TISSUE_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$tissue\t$tissue_id\n";
                }

             } else { # no tissue identifiers
                $invalid_tissue_count++;
                print STDERR "ERROR:Tissue $tissue given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_TISSUE_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$tissue\n";
             }

          } # end if tissue supplied 
          
          
        } # end else multiple mutation

     } else { # required data criteria not met, so find out which data are invalid 

        # The conditional statements below identify which field(s) are invalid for this PHI-base entry
        # Note that multiple problems may be associated with an individual annotation
  
        # add to general file of PHI-base accessions with invalid data
        $invalid_required_data_count++;
        print INVALID_REQUIRED_DATA_FILE "$required_fields_annot{'phi_base_acc'}\n";

        # detect invalid pathogen taxon ID
        if (not defined $required_fields_annot{'patho_tax'}
            or $required_fields_annot{'patho_tax'} eq "") {
           $invalid_path_taxon_id_count++;
           print INVALID_PATH_TAXON_FILE "$required_fields_annot{'phi_base_acc'}\n";
        } elsif ($required_fields_annot{'patho_tax'} !~ /^\d+$/) {
           $invalid_path_taxon_id_count++;
           print INVALID_PATH_TAXON_FILE "$required_fields_annot{'phi_base_acc'}\t$required_fields_annot{'patho_tax'}\n";
        }

        # detect invalid host taxon ID
        if (not defined $required_fields_annot{'host_tax'}
            or $required_fields_annot{'host_tax'} eq "") {
           $invalid_host_taxon_id_count++;
           print INVALID_HOST_TAXON_FILE "$required_fields_annot{'phi_base_acc'}\n";
        } elsif ($required_fields_annot{'host_tax'} !~ /^\d+$/) {
           $invalid_host_taxon_id_count++;
           print INVALID_HOST_TAXON_FILE "$required_fields_annot{'phi_base_acc'}\t$required_fields_annot{'host_tax'}\n";
        }

        # detect invalid PubMed ID
        if (not defined $required_fields_annot{'literature_id'} 
            or not defined $required_fields_annot{'literature_source'}
            or $required_fields_annot{'literature_id'} eq ""
            or $required_fields_annot{'literature_source'} eq "") {
           $invalid_literature_count++;
           print INVALID_LITERATURE_FILE "$required_fields_annot{'phi_base_acc'}\n";
        } elsif (lc $required_fields_annot{'literature_source'} ne "pubmed") {
           $invalid_literature_count++;
           print INVALID_LITERATURE_FILE "$required_fields_annot{'phi_base_acc'}\t$required_fields_annot{'literature_source'}\t$required_fields_annot{'literature_id'}\n";
        }

        # detect invalid UniProt accession
        if (not defined $required_fields_annot{'db_type'} 
            or not defined $required_fields_annot{'accession'}
            or $required_fields_annot{'db_type'} eq ""
            or $required_fields_annot{'accession'} eq "") {
           $invalid_uniprot_acc_count++;
           print INVALID_UNIPROT_FILE "$required_fields_annot{'phi_base_acc'}\n";
        } elsif (lc $required_fields_annot{'db_type'} ne "uniprot") { # database must be UniProt
           $invalid_uniprot_acc_count++;
           print INVALID_UNIPROT_FILE "$required_fields_annot{'phi_base_acc'}\t$required_fields_annot{'db_type'}\t$required_fields_annot{'accession'}\n";
        } elsif ( not ($required_fields_annot{'accession'} =~ /^[a-zA-Z\d]{6}$/ # UniProt entry must be either 6 or 10 alphanumeric characters
                        or $required_fields_annot{'accession'} =~ /^[a-zA-Z\d]{10}$/) ) {
           $invalid_uniprot_acc_count++;
           print INVALID_UNIPROT_FILE "$required_fields_annot{'phi_base_acc'}\t$required_fields_annot{'db_type'}\t$required_fields_annot{'accession'}\n";
        }

        # detect if no gene name is given
        if (not defined $required_fields_annot{'gene_name'} 
            or $required_fields_annot{'gene_name'} eq "") {
           $invalid_gene_name_count++;
           print INVALID_GENE_NAME_FILE "$required_fields_annot{'phi_base_acc'}\n";
        }

        # detect if no curator is provided
        if (not defined $required_fields_annot{'entered_by'} 
            or $required_fields_annot{'entered_by'} eq "") {
           $invalid_curator_count++;
           print INVALID_CURATOR_FILE "$required_fields_annot{'phi_base_acc'}\n";
        }


     }

   } else { # else PHI-base accession does not exist, or is not valid

     # add the PHI-base accession string to the invalid accession array
     push(@invalid_phibase_acc,$phi_base_annotation{"phi_base_acc"});

   }

   # print message for every 500th PHI-base annotation processed
   print "PHI-base annotations processed:$annotation_count\n" unless ($annotation_count % 500);

} # end of file

close (TSV_FILE);
close (PATH_PHENOTYPE_FILE);
close (INVALID_PATH_PHENOTYPE_FILE);
close (GO_WITH_EVID_FILE);
close (GO_WITHOUT_EVID_FILE);
close (INVALID_GO_FILE);
close (EXP_SPEC_TERM_FILE);
close (INVALID_EXP_SPEC_FILE);
close (HOST_RES_TERM_FILE);
close (INVALID_HOST_RES_FILE);
close (PHEN_OUTCOME_TERM_FILE);
close (INVALID_PHEN_OUTCOME_FILE);
close (DISEASE_TERM_FILE);
close (INVALID_DISEASE_FILE);
close (WITHOUT_DISEASE_FILE);
close (INVALID_REQUIRED_DATA_FILE);
close (INVALID_PATH_TAXON_FILE);
close (INVALID_HOST_TAXON_FILE);
close (INVALID_LITERATURE_FILE);
close (INVALID_UNIPROT_FILE);
close (INVALID_GENE_NAME_FILE);
close (INVALID_CURATOR_FILE);

# save the invalid PHI-base accessions to a separate file
my $invalid_accessions_filename = '../error/invalid_phibase_accessions.txt';
open (INVALID_FILE, "> $invalid_accessions_filename") or die "Error opening output file\n";
foreach my $invalid_accession (@invalid_phibase_acc) {
   print INVALID_FILE "$invalid_accession\n" if defined $invalid_accession;
}
close (INVALID_FILE);

print "\nProcess completed successfully.\n\n";

print "Total PHI-base annotations processed:$annotation_count\n\n";

print "Total PHI-base annotations with an invalid accession: ".scalar(@invalid_phibase_acc)."\n\n";

print "Total annotations with invalid required data: $invalid_required_data_count\n";
print "Annotations without a valid pathogen taxon ID: $invalid_path_taxon_id_count\n";
print "Annotations without a valid host taxon ID: $invalid_host_taxon_id_count\n";
print "Annotations without a valid PubMed ID: $invalid_literature_count\n";
print "Annotations without a valid UniProt accession: $invalid_uniprot_acc_count\n";
print "Annotations without a gene name: $invalid_gene_name_count\n";
print "Annotations without a curator: $invalid_curator_count\n\n";

print "Total PHI-base accessions meeting the criteria for the required data: ".scalar(keys %required_criteria_annotations)."\n\n";

print "Total curator entries: $curator_count\n";
print "Total species expert entries: $species_expert_count\n";
print "Total valid pathogen_phenotypes: $pathogen_phenotype_count\n";
print "Total invalid pathogen_phenotypes: $invalid_pathogen_phenotype_count\n";
print "Total GO annotations: $go_annotation_count\n";
print "GO annotations with evidence code: $go_with_evid_count\n";
print "GO annotations without evidence code: $go_without_evid_count\n";
print "Invalid GO annotations: $invalid_go_count\n";
print "Total anti-infective chemicals: $anti_infective_count\n";
print "Annotations without anti-infective: $no_anti_infective_count\n";
print "Total inducer chemicals: $inducer_count\n";
print "Annotations without inducer: $no_inducer_count\n";
print "Total annotations with a experiment specification: $exp_spec_count\n";
print "Total experiment specification terms: $exp_spec_term_count\n";
print "Invalid experiment specifications: $invalid_exp_spec_count\n";
print "Total annotations with a host response: $host_response_count\n";
print "Total host response terms: $host_response_term_count\n";
print "Invalid host responses: $invalid_host_response_count\n";
print "Total annotations with a phenotype outcome: $phenotype_outcome_count\n";
print "Total phenotype outcome terms: $phenotype_outcome_term_count\n";
print "Invalid phenotype outcomes: $invalid_phenotype_outcome_count\n";
print "Total disease names: $disease_count\n";
print "Valid disease terms: $disease_term_count\n";
print "Invalid disease terms: $invalid_disease_count\n";
print "Annotations without disease terms: $without_disease_count\n";
print "Total annotations with a tissue: $tissue_count\n";
print "Total tissue terms: $tissue_term_count\n";
print "Invalid tissues: $invalid_tissue_count\n\n";

print "Output file of accession string with an invalid PHI-base accession (if available): $invalid_accessions_filename\n\n";

print "Output file of all PHI-base annotations with invalid required data: $invalid_required_data_filename\n";
print "Output file of annotations with invalid pathogen taxon ID: $invalid_path_taxon_filename\n";
print "Output file of annotations with invalid host taxon ID: $invalid_host_taxon_filename\n";
print "Output file of annotations with invalid PubMed ID: $invalid_literature_filename\n";
print "Output file of annotations with invalid UniProt accession: $invalid_uniprot_filename\n";
print "Output file of annotations without a Gene Name: $invalid_gene_name_filename\n";
print "Output file of annotations without a curator: $invalid_curator_filename\n\n";

print "Output file of valid pathogen phenotypes: $pathogen_phenotype_filename\n";
print "Output file of invalid pathogen phenotypes: $invalid_pathogen_phenotype_filename\n";
print "Output file of GO annotations with evidence code: $go_with_evid_filename\n";
print "Output file of GO annotations without evidence code: $go_without_evid_filename\n";
print "Output file of invalid GO annotations: $invalid_go_filename\n";
print "Output file of experiment spec annotations: $exp_spec_term_filename\n";
print "Output file of invalid experiment specs: $invalid_exp_spec_filename\n";
print "Output file of host response annotations: $host_res_term_filename\n";
print "Output file of invalid host responses: $invalid_host_res_filename\n";
print "Output file of phenotype outcome annotations: $phen_outcome_term_filename\n";
print "Output file of invalid phenotype outcomes: $invalid_phen_outcome_filename\n";
print "Output file of disease annotations: $disease_term_filename\n";
print "Output file of invalid diseases: $invalid_disease_filename\n";
print "Output file of annotation without a diseases: $without_disease_filename\n";
print "Output file of tissue annotations: $tissue_term_filename\n";
print "Output file of invalid tissues: $invalid_tissue_filename\n\n";

