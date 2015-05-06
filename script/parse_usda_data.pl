#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open the tab separated values (TSV) version of the USDA table
my $usda_tsv_filename = '../input/usda_data.tsv';
open (TSV_FILE, $usda_tsv_filename) || die "Error opening input file\n";
print "Processing USDA data from $usda_tsv_filename...\n";

# first line gives the spreadsheet column headings
chomp(my $col_header_line = <TSV_FILE>);
my @col_headers = split(/\t/,$col_header_line);


# parse tab-separated file that maps pathogen lifestyle values of the spreadsheet
# to the identifier in the pathogen lifestyle ontology
# saving the value and identifier as key/value pairs in a hash
open (PATH_LIFESTYLE_MAPPINGS_FILE, "../mapping/pathogen_lifestyle_mapping.tsv") || die "Error opening input file\n";

# hash to map pathogen lifestyle text to ontology identifier
my %pathogen_lifestyle_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "pathogen lifestyle ontology identifier", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes the ontology identifier
while (<PATH_LIFESTYLE_MAPPINGS_FILE>) {
  chomp;
  my ($path_lifestyle_value,$path_lifestyle_ontology_id_list) = split(/\t/,$_);
  $pathogen_lifestyle_mapping{$path_lifestyle_value} = $path_lifestyle_ontology_id_list;
}
close (PATH_LIFESTYLE_MAPPINGS_FILE);


# parse tab-separated file that maps natural host values of the spreadsheet
# to the identifier in the natural host ontology
# saving the value and identifier as key/value pairs in a hash
open (NATURAL_HOST_MAPPINGS_FILE, "../mapping/natural_host_mapping.tsv") || die "Error opening input file\n";

# hash to map natural host text to ontology identifier
my %natural_host_mapping;

# each row of the file contains a "valid spreadsheet value"
# and corresponding "natural host ontology identifier", separated by tab
# separate these fields and save as key/value pairs in a hash
# where key becomes the valid value & value becomes the ontology identifier
while (<NATURAL_HOST_MAPPINGS_FILE>) {
  chomp;
  my ($natural_host_value,$natural_host_ontology_id_list) = split(/\t/,$_);
  $natural_host_mapping{$natural_host_value} = $natural_host_ontology_id_list;
}
close (NATURAL_HOST_MAPPINGS_FILE);

# create hash to store all USDA data
my %all_usda_data;

# counters to gather statistics
my $pathogen_lifestyle_term_count = 0;
my $natural_host_term_count = 0;

# go through each of the remaining lines of the tab-delimited file (each representing a single entry)
# save the values of each column to the approriate output file, then insert into the appropriate database table
while (<TSV_FILE>) {

   # each value is separated based on the tab, then saved as an element of the array
   my @entry_data = split(/\t/,$_);

   # initialise column iterator
   my $column_num=0;

   # hash to store values for the current USDA entry
   my %usda_entry_hash;

   # iterate through each column of the current USDA entry, saving it to the appropriate text file
   foreach my $entry_value (@entry_data) {

      $entry_value =~ s/^\s+//; # remove blank space from start
      $entry_value =~ s/\s+$//; # remove blank space from end

      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ../output/column/usda_column_".$col_headers[$column_num].".txt") or die "Error opening output file\n";
      print COLUMN_FILE "$entry_value\n";

      # add data to the annotation hash
      $usda_entry_hash{$col_headers[$column_num]} = $entry_value;

      # increment col
      $column_num++;

   } # end foreach column of current USDA entry

   # add the current USDA entry to hash of all entries
   $all_usda_data{$usda_entry_hash{"Pathogen Taxonomy NCBI Identifier"}} = {%usda_entry_hash};


   # replace string literal '' or 'none' with integer 0
   my $usda_interactions = $usda_entry_hash{"No of host interactions USDA Fungal database"};
   if ($usda_interactions eq "" or $usda_interactions eq "none") {
     $usda_interactions = 0;
   }
   # replace string literal 'no entry' with integer 0
   if ($usda_interactions eq "no entry") {
     $usda_interactions = 0;
   }

   # replace string literal '' or 'none' with integer 0
   my $plantwise_interactions = $usda_entry_hash{"No of host interactions Plantwise Knowledge Bank"};
   if ($plantwise_interactions eq "" or $plantwise_interactions eq "none") {
     $plantwise_interactions = 0;
   }
   # replace string literal 'multiple' with integer 0 - VALUE 'MULTIPLE IS NOT VALID FOR AN INTEGER FIELD'
   if ($plantwise_interactions eq "multiple") {
     $plantwise_interactions = 0;
   }


   my $num_plant_hosts = $usda_entry_hash{"No of plant hosts"};

   # prepare insert statement
   my $sql_statement;
   # if plant host data is not available (e.g. for pathogen infecting only animal hosts)
   # then leave out the number of plant hosts from the insertion statement
   # otherwise, include the plant host data
   if (not defined $num_plant_hosts or $num_plant_hosts eq "") {
      $sql_statement = qq(INSERT INTO usda (pathogen_taxon_id,host_kingdom,
                                            num_interactions_usda,num_interactions_plantwise)
                          VALUES ($usda_entry_hash{"Pathogen Taxonomy NCBI Identifier"},
                                  '$usda_entry_hash{"Host Kingdom"}',
                                  $usda_interactions,
                                  $plantwise_interactions)
                          RETURNING id;
                         );
   } else {
      $sql_statement = qq(INSERT INTO usda (pathogen_taxon_id,host_kingdom,num_plant_hosts,
                                            num_interactions_usda,num_interactions_plantwise)
                          VALUES ($usda_entry_hash{"Pathogen Taxonomy NCBI Identifier"},
                                  '$usda_entry_hash{"Host Kingdom"}',
                                  '$usda_entry_hash{"No of plant hosts"}',
                                  $usda_interactions,
                                  $plantwise_interactions)
                          RETURNING id;
                         );
   }

   my $sql_result = $db_conn->prepare($sql_statement);
   $sql_result->execute() or die $DBI::errstr;

   # retrieve the USDA id from the returned value of the insert
   my @row = $sql_result->fetchrow_array();
   my $usda_id = shift @row;


   # get the pathogen lifestyle
   my $pathogen_lifestyle_string = $usda_entry_hash{"Pathogen Lifestyle"};

   # if the pathogen lifestyle string is empty, no pathogen lifestyle has been supplied
   if (defined $pathogen_lifestyle_string and $pathogen_lifestyle_string ne "") {

      $pathogen_lifestyle_string =~ s/^\s+//; # remove blank space from start of string
      $pathogen_lifestyle_string =~ s/\s+$//; # remove blank space from end of string

      # using the pathogen lifestyle mappings,
      # get the appropriate ontology identifiers associated with the pathogen lifestyle
      my $pathogen_lifestyle_id_list = $pathogen_lifestyle_mapping{$pathogen_lifestyle_string};

      # if identifiers are present, then insert the appropriate
      # records into the usda_pathogen_lifestyle table
      if ($pathogen_lifestyle_id_list) {

	 $pathogen_lifestyle_id_list =~ s/^\s+//; # remove blank space from start of string
	 $pathogen_lifestyle_id_list =~ s/\s+$//; # remove blank space from end of string

	 # need to split list based on semi-colon delimiter
	 my @pathogen_lifestyle_ontology_ids = split(";",$pathogen_lifestyle_id_list);

	 # for each pathogen lifestyle id, insert data into usda_pathogen_lifestyle table,
	 # with foreign keys to the interaction_host table and the pathogen lifestyle ontology
	 foreach my $pathogen_lifestyle_id (@pathogen_lifestyle_ontology_ids) {
	   $pathogen_lifestyle_term_count++;
	   $sql_statement = qq(INSERT INTO usda_pathogen_lifestyle (usda_id, pathogen_lifestyle_id)
				 VALUES ($usda_id, '$pathogen_lifestyle_id');
			      );
	   $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
	 }

      } else { # no pathogen lifestyle identifiers
	 print STDERR "ERROR:Pathogen lifestyle $pathogen_lifestyle_string is not valid\n";
      }

   } # end if pathogen lifestyle supplied           


   # get the natural host
   my $natural_host_string = $usda_entry_hash{"Natural Host Type"};

   # if the natural host string is empty, no natural host has been supplied
   if (defined $natural_host_string and $natural_host_string ne "") {

      $natural_host_string =~ s/^\s+//; # remove blank space from start of string
      $natural_host_string =~ s/\s+$//; # remove blank space from end of string

      # using the natural host mappings,
      # get the appropriate ontology identifieris associated with the natural host
      my $natural_host_id_list = $natural_host_mapping{$natural_host_string};

      # if identifiers are present, then insert the appropriate
      # records into the usda_natural_host table
      if ($natural_host_id_list) {

	 $natural_host_id_list =~ s/^\s+//; # remove blank space from start of string
	 $natural_host_id_list =~ s/\s+$//; # remove blank space from end of string

	 # need to split list based on semi-colon delimiter
	 my @natural_host_ontology_ids = split(";",$natural_host_id_list);

	 # for each natural host id, insert data into usda_natural_host table,
	 # with foreign keys to the interaction_host table and the natural host ontology
	 foreach my $natural_host_id (@natural_host_ontology_ids) {
	   $natural_host_term_count++;
	   $sql_statement = qq(INSERT INTO usda_natural_host (usda_id, natural_host_id)
				 VALUES ($usda_id, '$natural_host_id');
			      );
	   $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;
	 }

      } else { # no natural host identifiers
	 print STDERR "ERROR:Natural host $natural_host_string is not valid\n";
      }

   } # end if natural host supplied           


} # end of USDA data file

close (TSV_FILE);

# save all of the USDA entry data to file
# formatted with the pathogen taxon ID name as the first row
# then a separate row for each column heading and value, separated by tab
# with a blank line between each USDA entry
my $all_data_filename = '../output/all_usda_data.txt';
open (ALL_DATA_FILE, "> $all_data_filename") or die "Error opening output file\n";
foreach my $usda_entry (keys %all_usda_data) {
   print ALL_DATA_FILE "USDA entry\t$usda_entry\n";
   foreach my $col_name (keys %{ $all_usda_data{$usda_entry} }) {
     print ALL_DATA_FILE "$col_name\t$all_usda_data{$usda_entry}{$col_name}\n";
   }
   print ALL_DATA_FILE "\n";
}
close (ALL_DATA_FILE);


# retrieve all of the data inserted into the usda table of phibase
# and save the file in a tab-delimited file, with appropriate headings
my $sql_query = qq(SELECT * FROM usda;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

my $usda_entry_count = 0; 

# open outuput file
my $db_usda_filename = '../output/database_usda_data.tsv'; 
open (DATABASE_USDA_FILE, "> $db_usda_filename") or die "Error opening output file\n";

print DATABASE_USDA_FILE "usda_id\tpathogen_taxon_id\thost_kingdom\tnum_plant_hosts\tnum_interactions_usda\tnum_interactions_plantwise\n";

while (my @row = $sql_stmt->fetchrow_array()) {
  $usda_entry_count++;
  foreach my $column (@row) {
     if (defined $column) {
       print DATABASE_USDA_FILE "$column\t";
     } else {
       print DATABASE_USDA_FILE "\t";
     }
  }
  print DATABASE_USDA_FILE "\n";
}

close (DATABASE_USDA_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Total USDA entries:$usda_entry_count\n";
print "Total Pathogen Lifestyle terms:$pathogen_lifestyle_term_count\n";
print "Total Natural Host terms:$natural_host_term_count\n";
print "Output file of all USDA data: $all_data_filename\n";
print "Tab-separated file of all entry data inserted into usda table: $db_usda_filename\n";

