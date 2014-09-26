#!/usr/bin/perl
use strict;
use warnings;

# load perl postgresql module
use DBI;

# define credentials for phibase database
my $db_name = "phibase";
my $db_user = "postgres";
my $db_host = "localhost";
my $db_pw = "Jake0001"; 

# connect to database
my $db_conn = DBI->connect("DBI:Pg:dbname=$db_name;host=$db_host","$db_user","$db_pw");

# parse text file that maps columns headings of the spreadsheet to database field names
# saving the column name and db field name as key/value pairs in a hash
open (COL_NAMES_FILE, "column2accession.txt") || die "Error opening input file\n";

# hash to map spreadsheet headings to PHI-base3 db field names
my %column_mapping;

# each row of the file contains a spreadsheet heading the corresponding db field name, separated by tab
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
open (TSV_FILE, "phi-base-1_vs36.txt") || die "Error opening input file\n";

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

# create hash to store all data
my %phi_base_data;
# create another hash for the subset with only required fields
my %required_fields_data;
# create another hash for fusarium graminearum data
my %fusarium_gram_data;

# go through each of the remaining lines of the TSV file (each representing a single annotation)
# save the values of each column to the approriate output file
while (<TSV_FILE>) {
   # each value is separated based on the tab, then saved as an element of the array
   my @phi_array = split(/\t/,$_);
   my $i=0;
   # iterate through each column of the annotation, saving it to the appropriate text file

   #create separate hash to store all values
   my %phi_base_annotation;

   # the name of the text file is determined by mapping the column header to db field name, using the hash created earlier
   foreach my $phi_value (@phi_array) {
      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ./output/column_".$column_mapping{$new_col_headers[$i]}.".txt") or die "Error opening output file\n";
      print COLUMN_FILE "$phi_value\n";

      # add data to the annotation hash
      $phi_base_annotation{$column_mapping{$new_col_headers[$i]}} = $phi_value;

      $i++;
      last if ($i == @new_col_headers);  # values after this are just blank, so exit out of foreach loop
   } # end foreach annotation

   # variable for PHI-base accession number
   my $phi_acc_num;

   # check if the PHI-base accession exists & if it has a valid prefix
   # then remove the PHI: prefix, leaving just the PHI-base accession number
   # and add the current annotation to the overall data hash, using the accession number as a key
   if ($phi_acc_num = $phi_base_annotation{"phi_base_acc"} and $phi_acc_num =~ /^PHI:/) {

     $phi_acc_num =~ s/PHI://;
     $phi_base_data{$phi_acc_num} = {%phi_base_annotation};

     # get subset of the hash containing only the required fields
     my @required_fields = (
     "phi_base_acc",
     "db_type",
     "accession",
     "gene_name",
     "patho_tax",
     "host_tax",
     "literature_id"
     );
     my %required_fields_annot;
     @required_fields_annot{@required_fields} = @phi_base_annotation{@required_fields};
     $required_fields_data{$phi_acc_num} = {%required_fields_annot};

     # get subset of these where pathogen taxonomy ID = 5518 (Fusarium graminearum),
     # with all required fields defined and a UniProt accession
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
          and $required_fields_annot{"patho_tax"} == 5518 ) {

        $fusarium_gram_data{$phi_acc_num} = {%required_fields_annot};

	#print "PHI-base Accession:$required_fields_annot{'phi_base_acc'}\n";
	#print "UniProt Accession:$required_fields_annot{'accession'}\n";
	#print "Gene Name: $required_fields_annot{'gene_name'}\n";
	#print "Pathogen Species NCBI Taxon ID:$required_fields_annot{'patho_tax'}\n";
	#print "Host Species NCBI Taxon ID:$required_fields_annot{'host_tax'}\n";
	#print "PubMed ID:$required_fields_annot{'literature_id'}\n\n";

=pod
	my $sql_query = qq(SELECT * FROM interaction WHERE phi_base_accession = '$required_fields_annot{"phi_base_acc"}';);
	my $sql_stmt = $db_conn->prepare($sql_query);
	my $sql_result = $sql_stmt->execute() or die $DBI::errstr; 

	while (my @row = $sql_stmt->fetchrow_array()) {
	  print "Interaction ID: ".$row[0]."\n";
	#  print "PHI-base Acc: ".$row[1]."\n";
	#  print "Date: ".$row[2]."\n\n" if defined $row[2];
	}
=cut

	# test insert statement
	my $sql_statement1 = qq(INSERT INTO interaction (phi_base_accession,curation_date) 
			       VALUES ('$required_fields_annot{"phi_base_acc"}',current_date) RETURNING id;);
	my $sql_statement2 = qq(INSERT INTO pathogen_gene (ncbi_taxon_id,gene_name) 
			       VALUES ('$required_fields_annot{"patho_tax"}','$required_fields_annot{"gene_name"}') RETURNING id;);

	my $sql_result1 = $db_conn->prepare($sql_statement1);
	$sql_result1->execute() or die $DBI::errstr;

	my $sql_result2 = $db_conn->prepare($sql_statement2);
	$sql_result2->execute() or die $DBI::errstr;

	my $sql_statement3;

	while ( my @row2 = $sql_result2->fetchrow_array() ) {
	  my $pathogen_gene_id = $row2[0];
	  #print "Pathogen_gene ID: ".$pathogen_gene_id."\n";
	  
	  $sql_statement3 = qq(INSERT INTO pathogen_gene_mutant (pathogen_gene_id,ncbi_taxon_id,uniprot_accession) 
			       VALUES ($pathogen_gene_id,'$required_fields_annot{"patho_tax"}',
                                      '$required_fields_annot{"accession"}') RETURNING id;);
	}

	my $sql_result3 = $db_conn->prepare($sql_statement3);
	$sql_result3->execute() or die $DBI::errstr;
	#print "pathogen_gene_mutant record inserted successfully\n";

	while ( my @row1 = $sql_result1->fetchrow_array() and my @row3 = $sql_result3->fetchrow_array() ) {
	  my $interaction_id = $row1[0];
	  my $pathogen_gene_mutant_id = $row3[0];
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
	} 

     } # end if pathogen id = fusarium gram 

   }  # end if PHI-base accession exists

} # end of file
close (TSV_FILE);

# save all of the valid phibase data to file (sorted by PHI-base accession)
# formatted with the PHI-base accession as the first row
# then a separate row for each column heading and value, separated by tab
# with a blank line between each PHI-base annotation
open (ALL_DATA_FILE, "> ./output/all_phibase_data.tsv") or die "Error opening output file\n";
foreach my $phi_base_ann (sort {$a<=>$b} keys %phi_base_data) {
   print ALL_DATA_FILE "PHI:$phi_base_ann\n";
   foreach my $col_name (sort keys %{ $phi_base_data{$phi_base_ann} }) {
     print ALL_DATA_FILE "$col_name\t$phi_base_data{$phi_base_ann}{$col_name}\n";
   }
   print ALL_DATA_FILE "\n";
}
close (ALL_DATA_FILE);

# save required data to a separate file, using same format as above
open (REQUIRED_FIELDS_FILE, "> ./output/required_phibase_data.tsv") or die "Error opening output file\n";
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

# save required data to a separate file, using same format as above
open (SPECIES_FILE, "> ./output/species_fusarium_gram_data.tsv") or die "Error opening output file\n";
my $interaction_counter = 0;
foreach my $phi_base_ann (sort {$a<=>$b} keys %fusarium_gram_data) {
   print SPECIES_FILE "PHI:$phi_base_ann\n";
   foreach my $col_name (sort keys %{ $fusarium_gram_data{$phi_base_ann} }) {
     # check that the value is defined before attempting to display it
     if (defined $fusarium_gram_data{$phi_base_ann}{$col_name}) {
       print SPECIES_FILE "$col_name\t$fusarium_gram_data{$phi_base_ann}{$col_name}\n";
     }
   }
   print SPECIES_FILE "\n";
   $interaction_counter++;
}
close (SPECIES_FILE);
print "Total valid interactions for Fusarium gram: $interaction_counter\n";


my $sql_query = qq(SELECT * FROM interaction, interaction_pathogen_gene_mutant, pathogen_gene_mutant, pathogen_gene
                   WHERE interaction.id = interaction_pathogen_gene_mutant.interaction_id
                   AND pathogen_gene_mutant.id = interaction_pathogen_gene_mutant.pathogen_gene_mutant_id
                   AND pathogen_gene.id = pathogen_gene_mutant.pathogen_gene_id
                 ;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr; 
while (my @row = $sql_stmt->fetchrow_array()) {
  my $col_count = 0;
  foreach my $column (@row) {
     $col_count++;
     print "Col $col_count: $column\n" if defined $column; 
  }
  print "\n";
}

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

