#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open the tab separated values (TSV) version of the Curator table
my $curator_tsv_filename = '../input/phibase_v3pt6_curators.tsv';
open (TSV_FILE, $curator_tsv_filename) || die "Error opening input file\n";
print "Processing curator data from $curator_tsv_filename...\n";

# first line gives the spreadsheet column headings
chomp(my $col_header_line = <TSV_FILE>);
my @col_headers = split(/\t/,$col_header_line);

# create hash to store all curator data
my %all_curator_data;

# counters for statistics
my $curator_count = 0;
my $curator_with_org_count = 0;
my $curator_without_org_count = 0;

# go through each of the remaining lines of the tab-delimited file (each representing a single curator)
# save the values of each column to the approriate output file, then insert into the appropriate database table
while (<TSV_FILE>) {

   # each value is separated based on the tab, then saved as an element of the array
   my @curator_data = split(/\t/,$_);

   # initialise column iterator
   my $column_num = 0;

   # hash to store values for the current curator
   my %curator_hash;

   # iterate through each column of the current curator, saving it to the appropriate text file
   foreach my $curator_details (@curator_data) {

      $curator_details =~ s/^\s+//; # remove blank space from start
      $curator_details =~ s/\s+$//; # remove blank space from end

      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ../output/column/curator_column_".$col_headers[$column_num].".txt") or die "Error opening output file\n";
      print COLUMN_FILE "$curator_details\n";

      # add data to the annotation hash
      $curator_hash{$col_headers[$column_num]} = $curator_details;

      # increment col
      $column_num++;

   } # end foreach column of current curator

   # add the current curator details to hash of all curators
   $all_curator_data{$curator_hash{"ID"}} = {%curator_hash};

   # declare SQL variables
   my $sql_statement;
   my $sql_result;

   if ($curator_hash{"Organisation"} eq "") {
      $curator_without_org_count++;
      print STDERR "No organisation available for curator $curator_hash{'Name'}\n";
      # prepare SQL statement to insert curator without organisation
      $sql_statement = qq(INSERT INTO curator (initials,name)
                          VALUES ('$curator_hash{"ID"}','$curator_hash{"Name"}');
                         );
   } else {

      $curator_with_org_count++;
      # statement to insert new organisation, if it does not already exist
      $sql_statement = qq(INSERT INTO curation_organisation (name) 
                             SELECT '$curator_hash{"Organisation"}'
                             WHERE NOT EXISTS (
                               SELECT 1 FROM curation_organisation
                               WHERE name = '$curator_hash{"Organisation"}'
                            ));

      $sql_result = $db_conn->prepare($sql_statement);
      $sql_result->execute() or die $DBI::errstr;

      # get the unique identifier for the curation_organisation record
      $sql_statement = qq(SELECT id FROM curation_organisation
                            WHERE name = '$curator_hash{"Organisation"}';
                         );
      $sql_result = $db_conn->prepare($sql_statement);
      $sql_result->execute() or die $DBI::errstr;
      my @row = $sql_result->fetchrow_array();
      my $curation_org_id = shift @row;

      # prepare SQL statement to insert curator with organisation id
      $sql_statement = qq(INSERT INTO curator (initials, name, curation_organisation_id)
                          VALUES ('$curator_hash{"ID"}','$curator_hash{"Name"}',$curation_org_id);
                         );

   } # end else organisation available

   # execute SQL statement to insert curator
   $sql_result = $db_conn->prepare($sql_statement);
   $sql_result->execute() or die $DBI::errstr;
   $curator_count++;

} # end of curator data file

close (TSV_FILE);


# save all of the curator data to file
# formatted with the common name as the first row
# then a separate row for each column heading and value, separated by tab
# with a blank line between each curator
my $all_data_filename = '../output/all_curator_data.txt';
open (ALL_DATA_FILE, "> $all_data_filename") or die "Error opening output file\n";
foreach my $curator_id (keys %all_curator_data) {
   print ALL_DATA_FILE "Curator ID\t$curator_id\n";
   foreach my $col_name (keys %{ $all_curator_data{$curator_id} }) {
     print ALL_DATA_FILE "$col_name\t$all_curator_data{$curator_id}{$col_name}\n";
   }
   print ALL_DATA_FILE "\n";
}
close (ALL_DATA_FILE);


# retrieve all of the data inserted into the curator and curator_organisation tables of phibase
# and save the values in a tab-delimited file, with appropriate headings
#my $sql_query = qq(SELECT * FROM curator, curation_organisation
#                   WHERE curator.curation_organisation_id = curation_organisation.id;
#                  );
my $sql_query = qq(SELECT * FROM curator;);
my $sql_stmt = $db_conn->prepare($sql_query);
my $sql_result = $sql_stmt->execute() or die $DBI::errstr;

# open outuput file
my $db_curator_filename = '../output/database_curator_data.tsv'; 
open (DATABASE_FRAC_FILE, "> $db_curator_filename") or die "Error opening output file\n";

#print DATABASE_FRAC_FILE "curator_id\tcurator_name\tcuration_org_id(FK)\tcuration_org_id\tcuration_org_name\n";
print DATABASE_FRAC_FILE "curator_id\tcurator_initials\tcurator_name\tcuration_org_id(FK)\n";

while (my @row = $sql_stmt->fetchrow_array()) {
  foreach my $column (@row) {
     print DATABASE_FRAC_FILE "$column\t" if defined $column; 
  }
  print DATABASE_FRAC_FILE "\n";
}

close (DATABASE_FRAC_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Curators entered into database:$curator_count\n";
print "Curators with organisataion:$curator_with_org_count\n";
print "Curators without organisation:$curator_without_org_count\n";
print "Output file of all curator data: $all_data_filename\n";
print "Tab-separated file of all curator data inserted into curator and curation_organisation tables: $db_curator_filename\n";

