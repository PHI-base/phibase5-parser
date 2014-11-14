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
my $phibase_tsv_filename = 'phi-base-1_vs36.txt';
open (TSV_FILE, $phibase_tsv_filename) || die "Error opening input file\n";
print "Processing PHI-base data from $phibase_tsv_filename...\n";
print "Inserting data for valid Fusarium graminearum annotations into PHI-base v5 database...\n";

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

# go through each of the remaining lines of the TSV file (each representing a single annotation)
# save the values of each column to the approriate output file
while (<TSV_FILE>) {

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
     "literature_source"
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
          and $required_fields_annot{"phi_base_acc"} ne ""
          and $required_fields_annot{"accession"} ne ""
          and $required_fields_annot{"gene_name"} ne ""
          and $required_fields_annot{"host_tax"} ne ""
          and $required_fields_annot{"literature_id"} ne ""
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

	# insert data into the pathogen_gene table
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

	my $sql_statement3;

        my $pathogen_gene_id;

	while ( my @row4 = $sql_result4->fetchrow_array() ) {
	  $pathogen_gene_id = $row4[0];
	  #print "Pathogen_gene ID: ".$pathogen_gene_id."\n";
	  
	  $sql_statement3 = qq(INSERT INTO pathogen_gene_mutant (pathogen_gene_id,ncbi_taxon_id,uniprot_accession) 
			       SELECT $pathogen_gene_id,'$required_fields_annot{"patho_tax"}',
                                      '$required_fields_annot{"accession"}'
                               WHERE NOT EXISTS (
                                 SELECT 1 FROM pathogen_gene_mutant
			         WHERE pathogen_gene_id = $pathogen_gene_id
                                 AND ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                 AND uniprot_accession = '$required_fields_annot{"accession"}'
                               ));
	}

	my $sql_result3 = $db_conn->prepare($sql_statement3);
	$sql_result3->execute() or die $DBI::errstr;
	#print "pathogen_gene_mutant record inserted successfully\n";

        my $sql_statement5 = qq(SELECT id FROM pathogen_gene_mutant
			        WHERE pathogen_gene_id = $pathogen_gene_id
                                AND ncbi_taxon_id = '$required_fields_annot{"patho_tax"}'
                                AND uniprot_accession = '$required_fields_annot{"accession"}');

	my $sql_result5 = $db_conn->prepare($sql_statement5);
	$sql_result5->execute() or die $DBI::errstr;

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

          while ( my @mult_mut_row = $sql_stmt->fetchrow_array() and my @row5 = $sql_result5->fetchrow_array() ) {
	     my $mult_mut_interaction_id = $mult_mut_row[0];
	     my $pathogen_gene_mutant_id = $row5[0];
	     print "Mult Mutant Partner Interaction ID: ".$mult_mut_interaction_id."\n";
	     print "Pathogen_gene_mutant ID: ".$pathogen_gene_mutant_id."\n";
	  
	     my $inner_sql_statement = qq(
				          INSERT INTO interaction_pathogen_gene_mutant (interaction_id,pathogen_gene_mutant_id) 
					    VALUES ($mult_mut_interaction_id,$pathogen_gene_mutant_id);
				         );
	     my $inner_sql_result = $db_conn->do($inner_sql_statement) or die $DBI::errstr;

	     print "Multiple mutation interaction_pathogen_gene_mutant record inserted successfully\n";
	  } # end while


        } else {  # annotation is not a multiple mutant, so insert new interaction records

#	  my $sql_statement1 = qq(INSERT INTO interaction (phi_base_accession,curation_date) 
#	                          VALUES ('$required_fields_annot{"phi_base_acc"}',current_date) RETURNING id;);
          # increment interaction counter, to become new PHI-base accession
          $interaction_num++;
          my $phi_base_accession = "PHI:I".$interaction_num;
	  my $sql_statement1 = qq(INSERT INTO interaction (phi_base_accession,curation_date) 
	                            VALUES ('$phi_base_accession',current_date) RETURNING id;);
	  my $sql_result1 = $db_conn->prepare($sql_statement1);
	  $sql_result1->execute() or die $DBI::errstr;

          # insert record to reference back to old phibase accession
          my $sql_statement6 = qq(INSERT INTO obsolete_reference (phi_base_accession,obsolete_accession)
                                    VALUES ('$phi_base_accession','$required_fields_annot{"phi_base_acc"}');
                                 );
	  my $sql_result6 = $db_conn->prepare($sql_statement6);
	  $sql_result6->execute() or die $DBI::errstr;

          while ( my @row1 = $sql_result1->fetchrow_array() and my @row5 = $sql_result5->fetchrow_array() ) {
	     my $interaction_id = $row1[0];
	     my $pathogen_gene_mutant_id = $row5[0];
	     #print "Interaction ID: ".$interaction_id."\n";
	     #print "Pathogen_gene_mutant ID: ".$pathogen_gene_mutant_id."\n";
	  
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
	  } # end while

        } # end else multiple mutation

     } # end if pathogen id = fusarium gram 

   } else { # else PHI-base accession does not exist, or is not valid
     # add the PHI-base accession string to the invalid accession array
     push(@invalid_phibase_acc,$phi_base_annotation{"phi_base_acc"});
   }

} # end of file
close (TSV_FILE);

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
                     pathogen_gene_mutant, pathogen_gene, interaction_literature, interaction_host
                   WHERE interaction.phi_base_accession = obsolete_reference.phi_base_accession 
                     AND interaction.id = interaction_pathogen_gene_mutant.interaction_id
                     AND pathogen_gene_mutant.id = interaction_pathogen_gene_mutant.pathogen_gene_mutant_id
                     AND pathogen_gene.id = pathogen_gene_mutant.pathogen_gene_id
                     AND interaction.id = interaction_literature.interaction_id
                     AND interaction.id = interaction_host.interaction_id
                 ;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

my $annotation_count = 0; 

my $db_data_filename = './output/database_data.tsv';  
open (DATABASE_DATA_FILE, "> $db_data_filename") or die "Error opening output file\n";

print DATABASE_DATA_FILE "int id\tPHI-base acc\tcuration date\tobsolete ref id\tobsolete ref new accession\tobsolete ref old accession\tint_path_gene_mut int_id\tint_path_gene_mut path_gene_id\tpath_gene_mutant_id\tpath_gene_mutant path_gene_id\tpath_gene_mutant ncbi_taxon_id\tuniprot_accession\tpath_gene id\tpath_gene ncbi_taxon_id\tpath_gene gene_name\tint_lit int_id\tint_lit PubMed ID\tint_host id\tint_host int_id\tint_host ncbi_taxon_id\n";

while (my @row = $sql_stmt->fetchrow_array()) {
  $annotation_count++;
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
print "Total annotations retrieved from database: $annotation_count\n\n";
print "Output file of all PHI-base annotations with valid data: $all_data_filename\n";
print "Output file of accession string with an invalid PHI-base accession (if available): $invalid_accessions_filename\n";
print "Output file of only the required data of valid PHI-base annotations: $required_data_filename\n";
print "Output file of valid data from fusarium graminearum: $fus_gram_filename\n";
print "Tab-separated file of all PHI-base data inserted into relevant tables: $db_data_filename\n\n";

