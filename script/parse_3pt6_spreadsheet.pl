#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot ontology_mapping); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database


# parse text file that maps columns headings of the spreadsheet to database field names
# saving the column name and db field name as key/value pairs in a hash
open (COL_NAMES_FILE, "../input/column2accession.txt") || die "Error opening input file\n";

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
  open (COL_TEXT_FILE, "> ../output/column/column_$db_field.txt") or die "Error opening output file column_$db_field\n";
  close (COL_TEXT_FILE);
}
close (COL_NAMES_FILE);

# Read in the relevant ontologies mapping files
print "Reading ontology mapping files...\n";

# parse tab-separated file that maps experimental evidence values of the spreadsheet
# to identifiers in the experiment specification ontology
# saving the value and identifiers as key/value pairs in a hash
open (EXP_SPEC_MAPPINGS_FILE, "../mapping/experiment_spec_mapping_phibase_3pt6.tsv") || die "Error opening input file\n";

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
# to identifiers in the host response ontology
# saving the value and identifiers as key/value pairs in a hash
open (HOST_RES_MAPPINGS_FILE, "../mapping/host_response_mapping_phibase_3pt6.tsv") || die "Error opening input file\n";

# hash to map host response text to ontology identifiers
my %host_response_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "host response ontology identifiers", separated by tab
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
# to the identifier in the phenotype outcome ontology
# saving the value and identifier as key/value pairs in a hash
open (PHEN_OUTCOME_MAPPINGS_FILE, "../mapping/phenotype_outcome_mapping_phibase_3pt6.tsv") || die "Error opening input file\n";

# hash to map phenotype outcome text to ontology identifier
my %phenotype_outcome_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "phenotype outcome ontology identifier", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes the ontology identifier
while (<PHEN_OUTCOME_MAPPINGS_FILE>) {
  chomp;
  my ($phen_outcome_value,$phen_outcome_ontology_id_list) = split(/\t/,$_);
  $phenotype_outcome_mapping{$phen_outcome_value} = $phen_outcome_ontology_id_list;
}
close (PHEN_OUTCOME_MAPPINGS_FILE);

# disease terms for PHI-base should include the combined terms from
# both the Plant Disease ontology and the [Human] Disease Ontology,
# so a hash of key/value pairs is created where the key is the ontology identifier
# and the value is the term name
my %combined_disease_mapping = (
      ontology_mapping('../ontology/Disease/PlantDisease/plant_disease_ontology.obo'),
      ontology_mapping('../ontology/Disease/HumanDisease/doid.obo')
   );

# open the tab separated values (TSV) version of the PHI-base spreadsheet
my $phibase_tsv_filename = '../input/phi-base-1_vs36_reduced_columns.tsv';
open (TSV_FILE, $phibase_tsv_filename) || die "Error opening input file\n";
print "Processing PHI-base data from $phibase_tsv_filename...\n";
print "Inserting data for valid Fusarium graminearum annotations into PHI-base v5 database...\n";

# open output files
my $defect_filename = '../output/phibase_defects.tsv';
my $invalid_defect_filename = '../error/phibase_invalid_defects.tsv';
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
my $invalid_required_data_filename = '../error/phibase_invalid_required_data.tsv';
my $invalid_path_taxon_filename = '../error/phibase_invalid_path_taxon_ids.tsv';
my $invalid_host_taxon_filename = '../error/phibase_invalid_host_taxon_ids.tsv';
my $invalid_literature_filename = '../error/phibase_invalid_literature.tsv';
my $invalid_uniprot_filename = '../error/phibase_invalid_uniprot_acc.tsv';
my $invalid_gene_name_filename = '../error/phibase_invalid_gene_names.tsv';
my $invalid_curator_filename = '../error/phibase_invalid_curators.tsv';
open (DEFECT_FILE, "> $defect_filename") or die "Error opening output file\n";
open (INVALID_DEFECT_FILE, "> $invalid_defect_filename") or die "Error opening output file\n";
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

# remove the last item from the array (which is blank)
pop(@col_headers);

# new array that will have blank header names removed
my @new_col_headers;

# iterate through the column names,
# saving only those that are not empty to the new array
foreach my $header (@col_headers) {
   chomp($header);
   push(@new_col_headers,$header) if $header;
}

# hash to store all annotations with a valid PHI-base accession
my %valid_phibase_data;
# array to store invalid PHI-base accessions
my @invalid_phibase_acc;

# create another hash for the subset with only required fields
my %required_fields_data;
# create another hash for fusarium graminearum data
my %fusarium_gram_data;

# counter for interaction number (used for new PHI-base accession)
my $interaction_num = 0;

# counters to gather statistics
my $annotation_count = 0;
my $defect_count = 0;
my $invalid_defect_count = 0;
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

   # each value is separated based on the tab, then saved as an element of the array
   my @phi_array = split(/\t/,$_);

   # initiate column iterator
   my $column_num=0;

   # hash to store all values of the current annotation
   my %phi_base_annotation;

   # iterate through each column of the annotation, saving its value to the appropriate text file
   # the name of the text file is determined by mapping the column header to db field name
   foreach my $phi_value (@phi_array) {

      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ../output/column/spreadsheet_column_".$column_mapping{$new_col_headers[$column_num]}.".txt")
         or die "Error opening output file\n";
      print COLUMN_FILE "$phi_value\n";

      # add data to the annotation hash
      $phi_base_annotation{$column_mapping{$new_col_headers[$column_num]}} = $phi_value;

      # increment the column number
      $column_num++;

      last if ($column_num == @new_col_headers);  # values after this are just blank, so exit out of foreach loop

   } # end foreach column of the annotation

   # variable for PHI-base accession number
   my $phi_acc_num;

   # check if the PHI-base accession exists & if it has a valid prefix
   # then remove the PHI: prefix, leaving just the PHI-base accession number
   # and add the current annotation to the overall data hash, using the accession number as a key
   if ($phi_acc_num = $phi_base_annotation{"phi_base_acc"} and $phi_acc_num =~ /^PHI:/) {

     $phi_acc_num =~ s/PHI://;
     $valid_phibase_data{$phi_acc_num} = {%phi_base_annotation};

     # get subset of the hash containing only the required fields
     my @required_fields = (
     "phi_base_acc",
     "db_type",
     "accession",
     "gene_name",
     "multiple_mutation",
     "patho_tax",
     "host_tax",
     "literature_id",
     "literature_source",
     "entered_by"
     );
     my %required_fields_annot;
     @required_fields_annot{@required_fields} = @phi_base_annotation{@required_fields};
     $required_fields_data{$phi_acc_num} = {%required_fields_annot};

     # get subset of these where pathogen taxonomy ID = 5518 (Fusarium graminearum),
     # with all required fields defined (except multiple mutation)
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
          and $required_fields_annot{"entered_by"} ne ""
          and lc $required_fields_annot{"db_type"} eq "uniprot"
          and lc $required_fields_annot{"literature_source"} eq "pubmed"
          and $required_fields_annot{"patho_tax"} =~ /^\d+$/  # taxon ID must be an integer
          #and $required_fields_annot{"patho_tax"} == 5518  # taxon ID for Fusarium gram
          #and $required_fields_annot{"patho_tax"} == 148305  # taxon ID for Magnaporthe oryzae
          #and $required_fields_annot{"patho_tax"} == 1307  # taxon ID for Streptococcus suis (causes human disease)
        ) {

        # add the required fields of the current annotation to the fusarium hash
        $fusarium_gram_data{$phi_acc_num} = {%required_fields_annot};

	#print "PHI-base Accession:$required_fields_annot{'phi_base_acc'}\n";
	#print "UniProt Accession:$required_fields_annot{'accession'}\n";
	#print "Gene Name: $required_fields_annot{'gene_name'}\n";
	#print "Pathogen Species NCBI Taxon ID:$required_fields_annot{'patho_tax'}\n";
	#print "Host Species NCBI Taxon ID:$required_fields_annot{'host_tax'}\n";
	#print "PubMed ID:$required_fields_annot{'literature_id'}\n\n";

	# insert data into the pathogen_gene table,
        # if it does not exist already (based on combination of taxon id and gene name)
	my $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id,gene_name) 
			        SELECT $required_fields_annot{"patho_tax"},'$required_fields_annot{"gene_name"}'
                                WHERE NOT EXISTS (
                                  SELECT 1 FROM pathogen_gene
                                  WHERE ncbi_taxon_id = $required_fields_annot{"patho_tax"}
                                  AND gene_name = '$required_fields_annot{"gene_name"}'
                               ));

	my $sql_result2 = $db_conn->prepare($sql_statement2);
	$sql_result2->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene record
        my $sql_statement4 = qq(SELECT id FROM pathogen_gene
                                WHERE ncbi_taxon_id = $required_fields_annot{"patho_tax"}
                                AND gene_name = '$required_fields_annot{"gene_name"}');

	my $sql_result4 = $db_conn->prepare($sql_statement4);
	$sql_result4->execute() or die $DBI::errstr;
	my @row4 = $sql_result4->fetchrow_array();
	my $pathogen_gene_id = shift @row4;

        # insert data into pathogen_gene_mutant table, including foreign key to pathogen_gene table 
	my $sql_statement3 = qq(INSERT INTO pathogen_gene_mutant (pathogen_gene_id,ncbi_taxon_id,uniprot_accession) 
			         SELECT $pathogen_gene_id,$required_fields_annot{"patho_tax"},
                                        '$required_fields_annot{"accession"}'
                                 WHERE NOT EXISTS (
                                   SELECT 1 FROM pathogen_gene_mutant
			           WHERE pathogen_gene_id = $pathogen_gene_id
                                   AND ncbi_taxon_id = $required_fields_annot{"patho_tax"}
                                   AND uniprot_accession = '$required_fields_annot{"accession"}'
                                 )
                               );

	my $sql_result3 = $db_conn->prepare($sql_statement3);
	$sql_result3->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene_mutant record
        my $sql_statement5 = qq(SELECT id FROM pathogen_gene_mutant
			        WHERE pathogen_gene_id = $pathogen_gene_id
                                AND ncbi_taxon_id = $required_fields_annot{"patho_tax"}
                                AND uniprot_accession = '$required_fields_annot{"accession"}');

	my $sql_result5 = $db_conn->prepare($sql_statement5);
	$sql_result5->execute() or die $DBI::errstr;
	my @row5 = $sql_result5->fetchrow_array();
	my $pathogen_gene_mutant_id = shift @row5;

        # before inserting a new interaction, we need to find out if the current PHI-base accession
        # should be part of an existing interaction (i.e. in a multiple mutation)

        # mutliple mutation flag
        my $multiple_mutation = 0;
        # PHI-base accession number of multiple mutation partner - MAY NEED TO CONVERT TO ARRAY TO COPE WITH MULTIPLE PARTNERS
        my $multi_mut_phi_acc_num;

        # check if the annotation part of a is a "multiple mutation" interaction
        if ($required_fields_annot{"multiple_mutation"} ne "") {

          print $required_fields_annot{"phi_base_acc"}."\t";
          print $required_fields_annot{"gene_name"}."\t";
          print $required_fields_annot{"accession"}."\t";
          print $required_fields_annot{"host_tax"}."\t";
          print $required_fields_annot{"multiple_mutation"}."\n";

          $multi_mut_phi_acc_num  = $required_fields_annot{"multiple_mutation"}; # MAY NEED TO SPLIT BASED ON SEMI-COLON
          $multi_mut_phi_acc_num  =~ s/PHI://;
	  # confirm if the multiple mutation partner gene already exists
          # only an annotation where the partner already exists needs to be treated differently from other annotations
	  if (exists $fusarium_gram_data{$multi_mut_phi_acc_num}) {
            $multiple_mutation = 1;
	    #print "Other annotation exists: $fusarium_gram_data{$multi_mut_phi_acc_num}{'phi_base_acc'}\n";
	  }

        } # end if multiple mutation


        if ($multiple_mutation) {
          print "In multiple mutation for: $required_fields_annot{'phi_base_acc'}, linking to existing $fusarium_gram_data{$multi_mut_phi_acc_num}{'phi_base_acc'}\n";
          # need to find the correct interaction_id for the corresponding multiple mutant gene
          # there could be several interactions for this gene, so needs to be based on a combination
          # of phi_base_acc + host_tax
          # then insert new interaction_pathogen_gene_mutant record, based on the returned interaction_id

          my $sql_query = qq(SELECT interaction.id 
                   FROM interaction, interaction_host
                   WHERE interaction.phi_base_accession = '$fusarium_gram_data{$multi_mut_phi_acc_num}{"phi_base_acc"}'
                   AND interaction.id = interaction_host.interaction_id
                   AND interaction_host.ncbi_taxon_id = $required_fields_annot{"host_tax"}
                 ;);
          # print "\n$sql_query\n";
          my $sql_stmt = $db_conn->prepare($sql_query);
          $sql_stmt->execute() or die $DBI::errstr;
	  my @mult_mut_row = $sql_stmt->fetchrow_array();
	  my $mult_mut_interaction_id = shift @mult_mut_row;

          if ( $mult_mut_interaction_id and $pathogen_gene_mutant_id ) {
	     print "Mult Mutant Partner Interaction ID: ".$mult_mut_interaction_id."\n";
	     print "Pathogen_gene_mutant ID: ".$pathogen_gene_mutant_id."\n";
	  
	     my $inner_sql_statement = qq(
				          INSERT INTO interaction_pathogen_gene_mutant (interaction_id,pathogen_gene_mutant_id) 
					    VALUES ($mult_mut_interaction_id,$pathogen_gene_mutant_id);
				         );
	     my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

	     print "Multiple mutation interaction_pathogen_gene_mutant record inserted successfully\n";
	  }


        } else {  # annotation is not a multiple mutant, so insert new interaction records

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

          # insert record to reference back to old phibase accession
          my $sql_statement6 = qq(INSERT INTO obsolete (phi_base_accession,obsolete_accession)
                                    VALUES ('$phi_base_accession','$required_fields_annot{"phi_base_acc"}');
                                 );
	  my $sql_result6 = $db_conn->prepare($sql_statement6);
	  $sql_result6->execute() or die $DBI::errstr;

          if ( $interaction_id and $pathogen_gene_mutant_id ) {
             # add records for the other tables associated with the interaction,
             # using the interaction id as a foreign key to the interaction table
	     my $inner_sql_statement = qq(
		  		          INSERT INTO interaction_literature (interaction_id,pubmed_id) 
					    VALUES ($interaction_id,'$required_fields_annot{"literature_id"}');
				          INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
				            VALUES ($interaction_id,$required_fields_annot{"host_tax"});
				          INSERT INTO interaction_pathogen_gene_mutant (interaction_id,pathogen_gene_mutant_id) 
					    VALUES ($interaction_id,'$pathogen_gene_mutant_id');
				         );
	     my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
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
	    $sql_result->execute() or die $DBI::errstr;

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
                                 SELECT $required_fields_annot{"patho_tax"}, $curator_id
                                 WHERE NOT EXISTS (
                                     SELECT 1 FROM species_expert
                                     WHERE ncbi_taxon_id = $required_fields_annot{"patho_tax"}
                                     AND curator_id = $curator_id
                                   )
                                );
	     $sql_result = $db_conn->prepare($sql_statement);
	     $sql_result->execute() or die $DBI::errstr;
          }


          # create hash of key/value pairs for attribute/value of all the defects
          my %defects = (
                           'Mating Defect'     => $phi_base_annotation{"mating_defect"},
                           'Pre-penetration'   => $phi_base_annotation{"prepenetration"},
                           'Penetration'       => $phi_base_annotation{"penetration"},
                           'Post-penetration'  => $phi_base_annotation{"postprepenetration"},
                           'Vegetative Spores' => $phi_base_annotation{"vegetative_spores"},
                           'Sexual Spores'     => $phi_base_annotation{"sexual_spores"},
                           'In Vitro Growth'   => $phi_base_annotation{"in_vitro_growth"},
                           'Spore Germination' => $phi_base_annotation{"spore_germination"}
                        );


          # for each of the defects, retrieve the id for the relevant defect,
          # then retrieve the id for each defect value (if available)
          # then insert a interaction_defect record,
          # based on combination of interaction id, defect attribute id and defect value id
          foreach my $defect_attribute (keys %defects)
          {

             my $defect_values_string = $defects{$defect_attribute};

             if (defined $defect_values_string and $defect_values_string ne "" and lc($defect_values_string) ne "none") {

		# get the id for the current defect attribute
		my $sql_statement = qq(SELECT id FROM defect_attribute
		               	         WHERE attribute = '$defect_attribute';
				      );
		my $sql_result = $db_conn->prepare($sql_statement);
		$sql_result->execute() or die $DBI::errstr;
		my @row = $sql_result->fetchrow_array();
		my $defect_attr_id = shift @row;

                # separate list of defect values, based on semi-colon delimiter
                my @defect_values = split (";",$defect_values_string);

                # insert interaction_defect record for each defect value
                foreach my $defect_value (@defect_values) {

                   $defect_value =~ s/^\s+//; # remove blank space from start of defect value
                   $defect_value =~ s/\s+$//; # remove blank space from end of defect value

		   # get the id for the current defect value, if it is a valid value
		   $sql_statement = qq(SELECT id FROM defect_value
		   	               WHERE lower(value) = lower('$defect_value');
				      );
		   $sql_result = $db_conn->prepare($sql_statement);
		   $sql_result->execute() or die $DBI::errstr;
		   @row = $sql_result->fetchrow_array();
		   my $defect_value_id = shift @row;

		   # insert data into interaction_defect table,
		   # using the interaction id, defect attribute id and defect value id as foreign keys 
		   if (defined $defect_value_id) {
                      $defect_count++;
		      $sql_statement = qq(INSERT INTO interaction_defect (interaction_id, defect_attribute_id, defect_value_id)
		  			    VALUES ($interaction_id, $defect_attr_id, $defect_value_id);
					 );
		      $sql_result = $db_conn->prepare($sql_statement);
		      $sql_result->execute() or die $DBI::errstr;
                      print DEFECT_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$defect_attribute\t$defect_value\n";
                   } else {
                      $invalid_defect_count++;
                      print INVALID_DEFECT_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$defect_attribute\t$defect_value\n";
                   }

                } # end foreach defect value

             } # end if defect values

          } # end foreach defect attribute

          
          # Get the GO annotations
          my $go_annot_string = $phi_base_annotation{"go_annotation"};

          # if the GO annotations string is empty, no GO annotations have been supplied
          if (defined $go_annot_string and $go_annot_string ne "") {

            # need to split list based on semi-colon delimiter
            my @go_annotations = split(";",$go_annot_string);

            # for each GO annotation, need to get the GO ID and the evidence code,
            # then insert an interaction_go_annotation record for this GO annotation
            foreach my $go_annotation (@go_annotations) {

               $go_annotation_count++;

               # entries sometimes have a name before the ID, separated by a colon,
               # other times no name is given before the CAS Registry ID,
               # so need to extract the ID, based on colon delimiter
               # (which will always be the last part)
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
	          $sql_statement = qq(INSERT INTO interaction_go_annotation (interaction_id, go_id, go_evidence_code)
                                        VALUES ($interaction_id, '$go_id', '$go_evid_code');
                                     );
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

               } else { # GO evidence code not supplied

                  # insert data into interaction_go_annotation table,
                  # with foreign key to the interaction table, but without GO evidence code 
	          $sql_statement = qq(INSERT INTO interaction_go_annotation (interaction_id, go_id)
                                        VALUES ($interaction_id, '$go_id');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);
	          $sql_result->execute() or die $DBI::errstr;

                  $go_without_evid_count++;
                  print GO_WITHOUT_EVID_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$go_id\n";
    
               } # else GO evidence code
              
             } # end foreach GO term

          } # if GO terms supplied 

          
          # anti-infective chemicals are given as CAS Registry IDs
          my $cas_string = $phi_base_annotation{"cas"};

          # if the CAS string is empty, no anti-infectives have been supplied
          if (defined $cas_string and $cas_string ne "") {

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
          if (defined $inducer_string and $inducer_string ne "") {

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
          if (defined $exp_evid_string and $exp_evid_string ne "") {

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
          if (defined $host_response and $host_response ne "") {

             $host_response_count++;

             $host_response =~ s/^\s+//; # remove blank space from start of string
             $host_response =~ s/\s+$//; # remove blank space from end of string

             # using the host response mappings,
             # get the appropriate ontology identifiers associated with the host response
             my $host_response_id_list = $host_response_mapping{$host_response};

             # if identifiers are present, then insert the appropriate
             # records into the interaction_host_response table
             if ($host_response_id_list) {

                my $sql_statement = qq(SELECT id FROM interaction_host
                                         WHERE interaction_id = $interaction_id
                                         AND ncbi_taxon_id = $required_fields_annot{"host_tax"}
                                      );

	        my $sql_result = $db_conn->prepare($sql_statement);
	        $sql_result->execute() or die $DBI::errstr;
	        my @row = $sql_result->fetchrow_array();
	        my $interaction_host_id = shift @row;

                $host_response_id_list =~ s/^\s+//; # remove blank space from start of string
                $host_response_id_list =~ s/\s+$//; # remove blank space from end of string

                # need to split list based on semi-colon delimiter
                my @host_res_ontology_ids = split(";",$host_response_id_list);

                # for each host response id, insert data into interaction_host_response table,
                # with foreign keys to the interaction_host table and the host response ontology
                foreach my $host_response_id (@host_res_ontology_ids) {
                  $host_response_term_count++;
	          $sql_statement = qq(INSERT INTO interaction_host_response (interaction_host_id, host_response_id)
                                        VALUES ($interaction_host_id, '$host_response_id');
                                     );
	          $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                  print HOST_RES_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$host_response\t$host_response_id\n";
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
          if (defined $phenotype_outcome_string and $phenotype_outcome_string ne "") {

             $phenotype_outcome_count++;

             $phenotype_outcome_string =~ s/^\s+//; # remove blank space from start of string
             $phenotype_outcome_string =~ s/\s+$//; # remove blank space from end of string

             # using the phenotype outcome mappings,
             # get the appropriate ontology identifier associated with the phenotype outcome
             my $phenotype_outcome_id = $phenotype_outcome_mapping{$phenotype_outcome_string};

             # if identifier is present, then insert the appropriate
             # record into the interaction_phenotype_outcome table
             if ($phenotype_outcome_id) {

                $phenotype_outcome_id =~ s/^\s+//; # remove blank space from start of string
                $phenotype_outcome_id =~ s/\s+$//; # remove blank space from end of string

                # insert data into the appropriate interaction_phenotype_outcome table,
                # with for a foreign key to the phenotype_outcome ontology
                $phenotype_outcome_term_count++;
	        $sql_statement = qq(INSERT INTO interaction_phenotype_outcome
                                      (interaction_id, phenotype_outcome_id)
                                      VALUES ($interaction_id, '$phenotype_outcome_id');
                                   );
	        $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
                print PHEN_OUTCOME_TERM_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$phenotype_outcome_string\t$phenotype_outcome_id\n";

             } else { # no phenotype outcome identifier
                $invalid_phenotype_outcome_count++;
                print STDERR "ERROR:Phenotype outcome $phenotype_outcome_string given for $required_fields_annot{'phi_base_acc'} is not valid\n";
                print INVALID_PHEN_OUTCOME_FILE "$phi_base_accession\t$required_fields_annot{'phi_base_acc'}\t$phenotype_outcome_string\n";
             }

          } # end if phenotype outcome supplied 
          
          
          # get the diseases for this interaction
          my $disease_string = $phi_base_annotation{"disease_name"};

          # if the disease string is empty, no disease has been supplied
          if (defined $disease_string and $disease_string ne "") {

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
close (DEFECT_FILE);
close (INVALID_DEFECT_FILE);

# save all of the valid phibase data to file (sorted by PHI-base accession)
# formatted with the PHI-base accession as the first row
# then a separate row for each column heading and value, separated by tab
# with a blank line between each PHI-base annotationi
my $all_data_filename = '../output/all_phibase_data.txt';
open (ALL_DATA_FILE, "> $all_data_filename") or die "Error opening output file\n";
foreach my $phi_base_ann (sort {$a<=>$b} keys %valid_phibase_data) {
   print ALL_DATA_FILE "PHI:$phi_base_ann\n";
   foreach my $col_name (sort keys %{ $valid_phibase_data{$phi_base_ann} }) {
     print ALL_DATA_FILE "$col_name\t$valid_phibase_data{$phi_base_ann}{$col_name}\n";
   }
   print ALL_DATA_FILE "\n";
}
close (ALL_DATA_FILE);

# save the invalid PHI-base accessions to a separate file
my $invalid_accessions_filename = '../error/invalid_phibase_accessions.txt';
open (INVALID_FILE, "> $invalid_accessions_filename") or die "Error opening output file\n";
foreach my $invalid_accession (@invalid_phibase_acc) {
   print INVALID_FILE "$invalid_accession\n" if defined $invalid_accession;
}
close (INVALID_FILE);

# save required data to a separate file
my $required_data_filename = '../output/required_phibase_data.txt';
open (REQUIRED_FIELDS_FILE, "> $required_data_filename") or die "Error opening output file\n";
foreach my $phi_base_ann (sort {$a<=>$b} keys %required_fields_data) {
   print REQUIRED_FIELDS_FILE "PHI:$phi_base_ann\n";
   foreach my $col_name (sort keys %{ $required_fields_data{$phi_base_ann} }) {
     # check that the value is defined before attempting to display it
     if (defined $required_fields_data{$phi_base_ann}{$col_name}) {
       print REQUIRED_FIELDS_FILE "$col_name\t$required_fields_data{$phi_base_ann}{$col_name}\n";
     }
   }
   print REQUIRED_FIELDS_FILE "\n";
}
close (REQUIRED_FIELDS_FILE);

# save all the fusarium gram data to file, using same format as above
my $fus_gram_filename = '../output/species_fusarium_gram_data.txt';
open (SPECIES_FILE, "> $fus_gram_filename") or die "Error opening output file\n";
foreach my $phi_base_ann (sort {$a<=>$b} keys %fusarium_gram_data) {
   print SPECIES_FILE "PHI:$phi_base_ann\n";
   foreach my $col_name (sort keys %{ $fusarium_gram_data{$phi_base_ann} }) {
     # check that the value is defined before attempting to display it
     if (defined $fusarium_gram_data{$phi_base_ann}{$col_name}) {
       print SPECIES_FILE "$col_name\t$fusarium_gram_data{$phi_base_ann}{$col_name}\n";
     }
   }
   print SPECIES_FILE "\n";
}
close (SPECIES_FILE);


print "\nProcess completed successfully.\n\n";
print "Total PHI-base annotations with valid data: ".scalar(keys %valid_phibase_data)."\n";
print "Total PHI-base annotations with an invalid accession: ".scalar(@invalid_phibase_acc)."\n";
print "Total valid interactions for Fusarium graminearum: ".scalar(keys %fusarium_gram_data)."\n";
print "Total curator entries for F gram: $curator_count\n";
print "Total species expert entries for F gram: $species_expert_count\n";
print "Total valid defects for F gram: $defect_count\n";
print "Total invalid defects for F gram: $invalid_defect_count\n";
print "Total GO annotations for F gram: $go_annotation_count\n";
print "GO annotations with evidence code for F gram: $go_with_evid_count\n";
print "GO annotations without evidence code for F gram: $go_without_evid_count\n";
print "Invalid GO annotations for F gram: $invalid_go_count\n";
print "Total anti-infective chemicals for F gram: $anti_infective_count\n";
print "Annotations without anti-infective for F gram: $no_anti_infective_count\n";
print "Total inducer chemicals for F gram: $inducer_count\n";
print "Annotations without inducer for F gram: $no_inducer_count\n";
print "Total annotations with a experiment specification for F gram: $exp_spec_count\n";
print "Total experiment specification terms for F gram: $exp_spec_term_count\n";
print "Invalid experiment specifications for F gram: $invalid_exp_spec_count\n";
print "Total annotations with a host response for F gram: $host_response_count\n";
print "Total host response terms for F gram: $host_response_term_count\n";
print "Invalid host responses for F gram: $invalid_host_response_count\n";
print "Total annotations with a phenotype outcome for F gram: $phenotype_outcome_count\n";
print "Total phenotype outcome terms for F gram: $phenotype_outcome_term_count\n";
print "Invalid phenotype outcomes for F gram: $invalid_phenotype_outcome_count\n";
print "Total disease terms for F gram: $disease_term_count\n";
print "Invalid disease terms for F gram: $invalid_disease_count\n";
print "Annotations without disease terms for F gram: $without_disease_count\n\n";

print "Total annotations with invalid required data: $invalid_required_data_count\n";
print "Annotations without a valid pathogen taxon ID: $invalid_path_taxon_id_count\n";
print "Annotations without a valid host taxon ID: $invalid_host_taxon_id_count\n";
print "Annotations without a valid PubMed ID: $invalid_literature_count\n";
print "Annotations without a valid UniProt accession: $invalid_uniprot_acc_count\n";
print "Annotations without a gene name: $invalid_gene_name_count\n";
print "Annotations without a curator: $invalid_curator_count\n\n";

print "Total annotations retrieved from database: $annotation_count\n\n";

print "Output file of all PHI-base annotations with valid data: $all_data_filename\n";
print "Output file of accession string with an invalid PHI-base accession (if available): $invalid_accessions_filename\n";
print "Output file of only the required data of valid PHI-base annotations: $required_data_filename\n";
print "Output file of valid data from fusarium graminearum: $fus_gram_filename\n";
print "Output file of valid defects from fusarium graminearum: $defect_filename\n";
print "Output file of invalid defects from fusarium graminearum: $invalid_defect_filename\n";
print "Output file of GO annotations with evidence code from fusarium graminearum: $go_with_evid_filename\n";
print "Output file of GO annotations without evidence code from fusarium graminearum: $go_without_evid_filename\n";
print "Output file of invalid GO annotations from fusarium graminearum: $invalid_go_filename\n";
print "Output file of experiment spec annotations from fusarium graminearum: $exp_spec_term_filename\n";
print "Output file of invalid experiment specs from fusarium graminearum: $invalid_exp_spec_filename\n";
print "Output file of host response annotations from fusarium graminearum: $host_res_term_filename\n";
print "Output file of invalid host responses from fusarium graminearum: $invalid_host_res_filename\n";
print "Output file of phenotype outcome annotations from fusarium graminearum: $phen_outcome_term_filename\n";
print "Output file of invalid phenotype outcomes from fusarium graminearum: $invalid_phen_outcome_filename\n";
print "Output file of disease annotations from fusarium graminearum: $disease_term_filename\n";
print "Output file of invalid diseases from fusarium graminearum: $invalid_disease_filename\n";
print "Output file of annotation without a diseases from fusarium graminearum: $without_disease_filename\n\n";

print "Output file of all PHI-base annotations with invalid required data: $invalid_required_data_filename\n";
print "Output file of annotations with invalid pathogen taxon ID: $invalid_path_taxon_filename\n";
print "Output file of annotations with invalid host taxon ID: $invalid_host_taxon_filename\n";
print "Output file of annotations with invalid PubMed ID: $invalid_literature_filename\n";
print "Output file of annotations with invalid UniProt accession: $invalid_uniprot_filename\n";
print "Output file of annotations without a Gene Name: $invalid_gene_name_filename\n";
print "Output file of annotations without a curator: $invalid_gene_name_filename\n\n";


