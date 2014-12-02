#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# parse text file that maps columns headings of the spreadsheet to database field names
# saving the column name and db field name as key/value pairs in a hash
open (COL_NAMES_FILE, "column2accession.txt") || die "Error opening input file\n";

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
  open (COL_TEXT_FILE, "> ./output/column_$db_field.txt") or die "Error opening output file column_$db_field\n";
  close (COL_TEXT_FILE);
}
close (COL_NAMES_FILE);

# open the tab separated values (TSV) version of the PHI-base spreadsheet
my $phibase_tsv_filename = 'phi-base-1_vs36_reduced_columns.tsv';
open (TSV_FILE, $phibase_tsv_filename) || die "Error opening input file\n";
print "Processing PHI-base data from $phibase_tsv_filename...\n";
print "Inserting data for valid Fusarium graminearum annotations into PHI-base v5 database...\n";

# open output files
my $defect_filename = './output/phibase_defects.tsv';
my $invalid_defect_filename = './output/phibase_invalid_defects.tsv';
my $go_with_evid_filename = './output/phibase_go_with_evid.tsv';
my $go_without_evid_filename = './output/phibase_go_without_evid.tsv';
my $invalid_go_filename = './output/phibase_invalid_go.tsv';
open (DEFECT_FILE, "> $defect_filename") or die "Error opening output file\n";
open (INVALID_DEFECT_FILE, "> $invalid_defect_filename") or die "Error opening output file\n";
open (GO_WITH_EVID_FILE, "> $go_with_evid_filename") or die "Error opening output file\n";
open (GO_WITHOUT_EVID_FILE, "> $go_without_evid_filename") or die "Error opening output file\n";
open (INVALID_GO_FILE, "> $invalid_go_filename") or die "Error opening output file\n";

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
my $go_annotation_count = 0;
my $go_with_evid_count = 0;
my $go_without_evid_count = 0;
my $invalid_go_count = 0;
my $exp_spec_count = 0;

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
      open (COLUMN_FILE, ">> ./output/column_".$column_mapping{$new_col_headers[$column_num]}.".txt")
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
          and $required_fields_annot{"accession"} ne ""
          and $required_fields_annot{"gene_name"} ne ""
          and $required_fields_annot{"host_tax"} ne ""
          and $required_fields_annot{"literature_id"} ne ""
          and $required_fields_annot{"entered_by"} ne ""
          and $required_fields_annot{"db_type"} eq "Uniprot"
          and lc $required_fields_annot{"literature_source"} eq "pubmed"
          and $required_fields_annot{"patho_tax"} == 5518 ) {

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
			        SELECT '$required_fields_annot{"patho_tax"}','$required_fields_annot{"gene_name"}'
                                WHERE NOT EXISTS (
                                  SELECT 1 FROM pathogen_gene
                                  WHERE ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                  AND gene_name = '$required_fields_annot{"gene_name"}'
                               ));

	my $sql_result2 = $db_conn->prepare($sql_statement2);
	$sql_result2->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene record
        my $sql_statement4 = qq(SELECT id FROM pathogen_gene
                                WHERE ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                AND gene_name = '$required_fields_annot{"gene_name"}');

	my $sql_result4 = $db_conn->prepare($sql_statement4);
	$sql_result4->execute() or die $DBI::errstr;
	my @row4 = $sql_result4->fetchrow_array();
	my $pathogen_gene_id = $row4[0];

        # insert data into pathogen_gene_mutant table, including foreign key to pathogen_gene table 
	my $sql_statement3 = qq(INSERT INTO pathogen_gene_mutant (pathogen_gene_id,ncbi_taxon_id,uniprot_accession) 
			         SELECT $pathogen_gene_id,'$required_fields_annot{"patho_tax"}',
                                        '$required_fields_annot{"accession"}'
                                 WHERE NOT EXISTS (
                                   SELECT 1 FROM pathogen_gene_mutant
			           WHERE pathogen_gene_id = $pathogen_gene_id
                                   AND ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                   AND uniprot_accession = '$required_fields_annot{"accession"}'
                                 )
                               );

	my $sql_result3 = $db_conn->prepare($sql_statement3);
	$sql_result3->execute() or die $DBI::errstr;

	# get the unique identifier for the inserted pathogen_gene_mutant record
        my $sql_statement5 = qq(SELECT id FROM pathogen_gene_mutant
			        WHERE pathogen_gene_id = $pathogen_gene_id
                                AND ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                AND uniprot_accession = '$required_fields_annot{"accession"}');

	my $sql_result5 = $db_conn->prepare($sql_statement5);
	$sql_result5->execute() or die $DBI::errstr;
	my @row5 = $sql_result5->fetchrow_array();
	my $pathogen_gene_mutant_id = $row5[0];

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
                   AND interaction_host.ncbi_taxon_id = '$required_fields_annot{"host_tax"}'
                 ;);
          # print "\n$sql_query\n";
          my $sql_stmt = $db_conn->prepare($sql_query);
          $sql_stmt->execute() or die $DBI::errstr;
	  my @mult_mut_row = $sql_stmt->fetchrow_array();
	  my $mult_mut_interaction_id = $mult_mut_row[0];

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
	  my $interaction_id = $row1[0];

          # insert record to reference back to old phibase accession
          my $sql_statement6 = qq(INSERT INTO obsolete_reference (phi_base_accession,obsolete_accession)
                                    VALUES ('$phi_base_accession','$required_fields_annot{"phi_base_acc"}');
                                 );
	  my $sql_result6 = $db_conn->prepare($sql_statement6);
	  $sql_result6->execute() or die $DBI::errstr;

          if ( $interaction_id and $pathogen_gene_mutant_id ) {
	     #print "Interaction ID: ".$interaction_id."\n";
	     #print "Pathogen_gene_mutant ID: ".$pathogen_gene_mutant_id."\n";
	  
             # add records for the other tables associated with the interaction,
             # using the interaction id as a foreign key to the interaction table
	     my $inner_sql_statement = qq(
		  		          INSERT INTO interaction_literature (interaction_id,pubmed_id) 
					    VALUES ($interaction_id,'$required_fields_annot{"literature_id"}');
				          INSERT INTO interaction_host (interaction_id,ncbi_taxon_id) 
				            VALUES ($interaction_id,'$required_fields_annot{"host_tax"}');
				          INSERT INTO interaction_pathogen_gene_mutant (interaction_id,pathogen_gene_mutant_id) 
					    VALUES ($interaction_id,'$pathogen_gene_mutant_id');
				         );
	     my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;
	     #print "Interaction_literature, interaction_host, interaction_pathogen_gene_mutant records inserted successfully\n";
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
                                 SELECT '$required_fields_annot{"patho_tax"}', $curator_id
                                 WHERE NOT EXISTS (
                                     SELECT 1 FROM species_expert
                                     WHERE ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
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
                      print DEFECT_FILE "$required_fields_annot{'phi_base_acc'}\t$defect_attribute\t$defect_value\n";
                   } else {
                      $invalid_defect_count++;
                      print INVALID_DEFECT_FILE "$required_fields_annot{'phi_base_acc'}\t$defect_attribute\t$defect_value\n";
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
            # then insert an interaction_go_term record for this GO annotation
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

                  # insert data into interaction_go_term table,
                  # with foreign keys to the interaction table and the go_evidence_code_table 
	          $sql_statement = qq(INSERT INTO interaction_go_term (interaction_id, go_id, go_evidence_code)
                                        VALUES ($interaction_id, '$go_id', '$go_evid_code');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);

                  # execute SQL insert and, if successful, print to file
	          if ( $sql_result->execute() ) {
                     $go_with_evid_count++;
                     print GO_WITH_EVID_FILE "$required_fields_annot{'phi_base_acc'}\t$go_id\t$go_evid_code\n";
                  } else {  # SQL insert unsuccessful, then log error
                     $invalid_go_count++;
                     print STDERR "\nPHI-base ERROR: Evidence code $go_evid_code is not valid for $required_fields_annot{'phi_base_acc'}, $go_id\n\n";
                     print INVALID_GO_FILE "$required_fields_annot{'phi_base_acc'}\t$go_id\t$go_evid_code\n";
                  }

               } else { # GO evidence code not supplied

                  # insert data into interaction_go_term table,
                  # with foreign key to the interaction table, but without GO evidence code 
	          $sql_statement = qq(INSERT INTO interaction_go_term (interaction_id, go_id)
                                        VALUES ($interaction_id, '$go_id');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);
	          $sql_result->execute() or die $DBI::errstr;

                  $go_without_evid_count++;
                  print GO_WITHOUT_EVID_FILE "$required_fields_annot{'phi_base_acc'}\t$go_id\n";
    
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

          } # if anti-infectives supplied 

          
          # get the experimental evidence
          my $exp_evid_string = $phi_base_annotation{"experimental_evidence"};

          # if the CAS string is empty, no anti-infectives have been supplied
          if (defined $exp_evid_string and $exp_evid_string ne "") {

            # need to split list based on semi-colon delimiter
            my @exp_evid_entries = split(";",$exp_evid_string);

            # for each anti-infective, need to get the CAS Registry ID, then insert
            # a chemical record (if it does not already exist) and 
            # an interaction_anti-infective_chemical record for this interaction
            foreach my $exp_evid_entry (@exp_evid_entries) {

               $exp_spec_count++;

               # entries sometimes have a name before the ID, separated by a colon,
               # other times no name is given before the CAS Registry ID,
               # so need to extract the ID, based on colon delimiter
               # (which will always be the last part)
               #my @exp_evid_parts = split(":",$exp_evid_entry);
               #my $exp_evid_id = pop(@exp_evid_parts);

               $exp_evid_entry =~ s/^\s+//; # remove blank space from start of CAS ID
               $exp_evid_entry =~ s/\s+$//; # remove blank space from end of CAS ID

               my $exp_spec_id;
               my $exp_spec_id2;

               # swap the old experimental evidence string for an
               # equivalent experiment specification ontology identifier
               #for ($exp_evid_entry) {
               if ($exp_evid_entry eq 'gene disruption') { $exp_spec_id = "ESO:0000001" }
               elsif ($exp_evid_entry eq 'gene deletion') { $exp_spec_id = "ESO:0000002" }
               elsif ($exp_evid_entry eq 'altered gene expression / gene regulation') { $exp_spec_id = "ESO:0000024" }
               elsif ($exp_evid_entry eq 'altered gene expression / gene regulation: overexpression') { $exp_spec_id = "ESO:0000011" }
               elsif ($exp_evid_entry eq 'altered gene expression / gene regulation: downregulation') { $exp_spec_id = "ESO:0000012" }
               elsif ($exp_evid_entry eq 'altered gene expression / gene regulation: down- and upregulation') 
                       { $exp_spec_id = "ESO:0000011"; $exp_spec_id2 = "ESO:0000012" }
               elsif ($exp_evid_entry eq 'altered gene expression / gene regulation: silencing') { $exp_spec_id = "ESO:0000013" }
               elsif ($exp_evid_entry eq 'biochemical analysis') { $exp_spec_id = "ESO:0000005" }
               elsif ($exp_evid_entry eq 'functional test in host') { $exp_spec_id = "ESO:0000006" }
               elsif ($exp_evid_entry eq 'functional test in host: direct injection') { $exp_spec_id = "ESO:0000014" }
               elsif ($exp_evid_entry eq 'functional test in host: transient expression') { $exp_spec_id = "ESO:0000015" }
               elsif ($exp_evid_entry eq 'mutation' or 'gene mutation') { $exp_spec_id = "ESO:00000019" }
               elsif ($exp_evid_entry eq 'mutation: characterised' or 'gene mutation: characterised') { $exp_spec_id = "ESO:00000016" }
               elsif ($exp_evid_entry eq 'complementation') { $exp_spec_id = "ESO:0000007" }
               elsif ($exp_evid_entry eq 'sequence analysis of sensitive and resistant strains') { $exp_spec_id = "ESO:0000017" }
               elsif ($exp_evid_entry eq 'sexual cross, sequencing of resistance conferring allele') { $exp_spec_id = "ESO:0000018" }
               elsif ($exp_evid_entry eq 'other evidence') { $exp_spec_id = "ESO:0000010" }
               #}
 
               # insert data into interaction_experiment_spec table,
               # with foreign keys to the interaction table and the experiment spec ontology
               # finding out if two, one, or no experiment specifications are to be entered
               if (defined $exp_spec_id2) {
                  # need two insert statements, one for each experiment spec id
	          $sql_statement = qq(INSERT INTO interaction_experiment_spec (interaction_id, experiment_spec_id)
                                        VALUES ($interaction_id, '$exp_spec_id');
	                              INSERT INTO interaction_experiment_spec (interaction_id, experiment_spec_id)
                                        VALUES ($interaction_id, '$exp_spec_id2');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);
	          $sql_result->execute() or die $DBI::errstr;
               } elsif (defined $exp_spec_id) {
                  # need only one insert statement
	          $sql_statement = qq(INSERT INTO interaction_experiment_spec (interaction_id, experiment_spec_id)
                                        VALUES ($interaction_id, '$exp_spec_id');
                                     );
	          $sql_result = $db_conn->prepare($sql_statement);
	          $sql_result->execute() or die $DBI::errstr;
               } else {
                  # the experimental evidence is not valid
                  print STDERR "ERROR: Experimental evidence $exp_evid_entry given for $required_fields_annot{'phi_base_acc'} is not valid\n";
               }
              
             } # end foreach experimental evidence

          } # if experimental evidence supplied 

        } # end else multiple mutation

     } # end if required criteria met and pathogen id = fusarium gram 

     # TODO: ELSE (REQURIED FIELDS CRITERIA NOT MET - OUTPUT PHI-BASE ACCESSIONS TO SEPARATE FILE)
     # ALSO POSSIBLE ADD ADDITIONAL IF STATEMENTS TO IDENTIFY EXACTLY WHAT THE PROLEM(S) MAY BE
     # SEPERATE IFs FOR EACH CRITERIA (RATHER THAN ONE IF ELSE), SO THAT ANY PHI-BASE ENTRY
     # WITH MULTIPLE PROBLEMS CAN HAVE ALL OF THOSE PROBLEMS IDENTIFIED (BY BEING PRESENT IN EACH OF THE ERROR FILES)

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
my $all_data_filename = './output/all_phibase_data.txt';
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
my $invalid_accessions_filename = './output/invalid_phibase_accessions.txt';
open (INVALID_FILE, "> $invalid_accessions_filename") or die "Error opening output file\n";
foreach my $invalid_accession (@invalid_phibase_acc) {
   print INVALID_FILE "$invalid_accession\n" if defined $invalid_accession;
}
close (INVALID_FILE);

# save required data to a separate file
my $required_data_filename = './output/required_phibase_data.txt';
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
my $fus_gram_filename = './output/species_fusarium_gram_data.txt';
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

# retrieve all of the data inserted into the relevant tables of phibase
# and save the file in a tab-delimited file, with appropriate headings
my $sql_query = qq(SELECT * FROM interaction, obsolete_reference, interaction_pathogen_gene_mutant, 
                     pathogen_gene_mutant, pathogen_gene, interaction_literature, interaction_host,
                     interaction_curator, curator, curation_organisation, species_expert
                   WHERE interaction.phi_base_accession = obsolete_reference.phi_base_accession 
                     AND interaction.id = interaction_pathogen_gene_mutant.interaction_id
                     AND pathogen_gene_mutant.id = interaction_pathogen_gene_mutant.pathogen_gene_mutant_id
                     AND pathogen_gene.id = pathogen_gene_mutant.pathogen_gene_id
                     AND interaction.id = interaction_literature.interaction_id
                     AND interaction.id = interaction_host.interaction_id
                     AND interaction.id = interaction_curator.interaction_id
                     AND interaction_curator.curator_id = curator.id
                     AND curation_organisation.id = curator.curation_organisation_id
                 ;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

my $db_data_filename = './output/database_data.tsv';  
open (DATABASE_DATA_FILE, "> $db_data_filename") or die "Error opening output file\n";

print DATABASE_DATA_FILE "int id\tPHI-base acc\tcuration date\tobsolete ref id\tobsolete ref new accession\tobsolete ref old accession\tint_path_gene_mut int_id\tint_path_gene_mut path_gene_id\tpath_gene_mutant_id\tpath_gene_mutant path_gene_id\tpath_gene_mutant ncbi_taxon_id\tuniprot_accession\tpath_gene id\tpath_gene ncbi_taxon_id\tpath_gene gene_name\tint_lit int_id\tint_lit PubMed ID\tint_host id\tint_host int_id\tint_host ncbi_taxon_id\tint_cur int_id\tint_cur cur_id\tcurator id\tcurator init\tcurator name\tcurator org_id\tcur_org id\tcur_org name\tsp_expert taxon_id\tsp_expert curator_id\n";

while (my @row = $sql_stmt->fetchrow_array()) {
  foreach my $column (@row) {
     print DATABASE_DATA_FILE "$column\t" if defined $column; 
  }
  print DATABASE_DATA_FILE "\n";
}
close (DATABASE_DATA_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

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
print "Total experiment specifications for F gram: $exp_spec_count\n";
print "Total annotations retrieved from database: $annotation_count\n\n";

print "Output file of all PHI-base annotations with valid data: $all_data_filename\n";
print "Output file of accession string with an invalid PHI-base accession (if available): $invalid_accessions_filename\n";
print "Output file of only the required data of valid PHI-base annotations: $required_data_filename\n";
print "Output file of valid data from fusarium graminearum: $fus_gram_filename\n";
print "Output file of valid defects from fusarium graminearum: $defect_filename\n";
print "Output file of invalid defects from fusarium graminearum: $invalid_defect_filename\n";
print "Output file of GO annotations with evidence code from fusarium graminearum: $go_with_evid_filename\n";
print "Output file of GO annotations without evidence code from fusarium graminearum: $go_without_evid_filename\n";
print "Output file of Invalid GO annotations from fusarium graminearum: $invalid_go_filename\n";
print "Tab-separated file of all PHI-base data inserted into relevant tables: $db_data_filename\n\n";

