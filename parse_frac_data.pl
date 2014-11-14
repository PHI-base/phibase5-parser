#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open the tab separated values (TSV) version of the FRAC table
my $frac_tsv_filename = 'frac_data.tsv';
open (TSV_FILE, $frac_tsv_filename) || die "Error opening input file\n";
print "Processing FRAC data from $frac_tsv_filename...\n";

# first line gives the spreadsheet column headings
chomp(my $col_header_line = <TSV_FILE>);
my @col_headers = split(/\t/,$col_header_line);

# specify array for columns that must have a unique value
my @unique_columns = ("common_name","cas_registry","chebi_id");

# create hash to store all FRAC data
my %all_frac_data;

# holder for previous chemical in list
my %previous_chem_hash;

# go through each of the remaining lines of the tab-delimited file (each representing a single chemical)
# save the values of each column to the approriate output file
while (<TSV_FILE>) {

   # each value is separated based on the tab, then saved as an element of the array
   my @chemical_data = split(/\t/,$_);

   # initialise column iterator
   my $column_num=0;

   # hash to store values for the current FRAC chemical
   my %frac_chem_hash;

   # iterate through each column of the current FRAC chemical, saving it to the appropriate text file
   foreach my $chem_value (@chemical_data) {

      $chem_value =~ s/^\s+//; # remove blank space from start
      $chem_value =~ s/\s+$//; # remove blank space from end

      # if the current value is empty (as will often be the case)
      # this is because the value should be the same as the previous chemical,
      # so the corresponding value from the previous FRAC chemical is copied to the current one
      # (the exceptions are common name, CAS registry, and ChEBI ID, which must be unique and should be left empty)
      if ( $chem_value eq "" and not $col_headers[$column_num] ~~ @unique_columns ) {
         $chem_value = $previous_chem_hash{$col_headers[$column_num]};
      }

      # add data to output file for the individual column
      open (COLUMN_FILE, ">> ./output/frac_column_".$col_headers[$column_num].".txt") or die "Error opening output file\n";
      print COLUMN_FILE "$chem_value\n";

      # add data to the annotation hash
      $frac_chem_hash{$col_headers[$column_num]} = $chem_value;

      # increment col
      $column_num++;

   } # end foreach column of current FRAC chem

   # add the current FRAC chem to hash of all chemicals
   $all_frac_data{$frac_chem_hash{"common_name"}} = {%frac_chem_hash};

   # before next FRAC chemical, assign current chemical to previous chemical hash
   %previous_chem_hash = %frac_chem_hash;

   # prepare insert statement
   my $sql_statement = qq(INSERT INTO frac (frac_code,moa_code,moa_name,target_code,target_site,
                                      group_name,chemical_group,common_name,resistance_risk,comments)
                          VALUES ('$frac_chem_hash{"frac_code"}','$frac_chem_hash{"moa_code"}',
                                  '$frac_chem_hash{"moa_name"}','$frac_chem_hash{"target_code"}',
                                  '$frac_chem_hash{"target_site"}','$frac_chem_hash{"group_name"}',
                                  '$frac_chem_hash{"chemical_group"}','$frac_chem_hash{"common_name"}',
                                  '$frac_chem_hash{"resistance_risk"}','$frac_chem_hash{"comment"}')
                          RETURNING id;
                         );

   my $sql_result = $db_conn->prepare($sql_statement);
   $sql_result->execute() or die $DBI::errstr;

   # hash of known mode in planta, manually curated for phibase_v3pt6,
   # only available for a few previously annotated FRAC chemicals
   # for convenience, these are hard-coded in this hash
   my %mode_in_planta_hash = 
      (
        "3405" => "Systemic",
        "3015" => "Systemic",
        "4520" => "Systemic",
        "3392" => "Systemic",
        "35014" => "Systemic",
        "83295" => "Systemic",
        "81960" => "Systemic",
        "50145" => "Systemic",
        "9700" => "Systemic",
        "81917" => "Systemic",
        "9242" => "Systemic",
        "81853" => "Locally systemic",
        "81736" => "Systemic",
        "7451" => "Various",
        "9448" => "Various",
        "40909" => "Preventive, curative, and systemic",
        "28909" => "Preventive, early curative, and systemic"
        # "????",  #This chemical still needs to be curated into the ChEBI database (as of 13/11/2014, CAS reg: 94361-06-5)
      );

   # retrieve the FRAC id from the returned value of the insert
   my @row = $sql_result->fetchrow_array();
   my $frac_id = $row[0];

   # retreive the mode_of_planta, if availble for the current ChEBI ID
   my $mode_in_planta  = $mode_in_planta_hash{ $frac_chem_hash{"chebi_id"} };

   # if the FRAC chemcial has a defined mode_in_planta value
   # then insert should include mode_in_planta column,
   # otherwise use insert statement without this optional column
   if ( defined $mode_in_planta ) {
       #print "ChEBI ID: $frac_chem_hash{'chebi_id'}, Mode in planta: $mode_in_planta\n"; 
       $sql_statement = qq(INSERT INTO chemical (chebi_id, cas_registry, frac_id, mode_in_planta) 
                            VALUES ('$frac_chem_hash{"chebi_id"}', '$frac_chem_hash{"cas_registry"}',
                                     $frac_id, '$mode_in_planta');
                          );
   } else {
       # create insert statement without mode_in_planta
       $sql_statement = qq(INSERT INTO chemical (chebi_id,cas_registry,frac_id) 
                            VALUES ('$frac_chem_hash{"chebi_id"}','$frac_chem_hash{"cas_registry"}',$frac_id);
                          );
   }

   $sql_result = $db_conn->prepare($sql_statement);
   $sql_result->execute() or die $DBI::errstr;

} # end of FRAC data file

close (TSV_FILE);

# insert additional chemicals that are not part of FRAC fungicide data
# (only 2 available from PHI-base v3pt6, which are hard-coded here)
# obviously, no FRAC ID is included (these chemicals are not in the FRAC table)
my $sql_statement = qq(
                       INSERT INTO chemical (chebi_id, cas_registry, mode_in_planta) 
                          VALUES ('47519', '65277-42-1','Preventive and curative');
                       INSERT INTO chemical (chebi_id, cas_registry, mode_in_planta) 
                          VALUES ('46081', '86386-73-4','Primarily fungistatic');
                      );
my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

# save all of the FRAC chemical data to file
# formatted with the common name as the first row
# then a separate row for each column heading and value, separated by tab
# with a blank line between each FRAC chemical
my $all_data_filename = './output/all_frac_data.txt';
open (ALL_DATA_FILE, "> $all_data_filename") or die "Error opening output file\n";
foreach my $frac_chem (keys %all_frac_data) {
   print ALL_DATA_FILE "FRAC chemical\t$frac_chem\n";
   foreach my $col_name (keys %{ $all_frac_data{$frac_chem} }) {
     print ALL_DATA_FILE "$col_name\t$all_frac_data{$frac_chem}{$col_name}\n";
   }
   print ALL_DATA_FILE "\n";
}
close (ALL_DATA_FILE);


# retrieve all of the data inserted into the chemical and frac tables of phibase
# and save the file in a tab-delimited file, with appropriate headings
my $sql_query = qq(SELECT * FROM chemical, frac
                   WHERE chemical.frac_id = frac.id 
                 ;);
my $sql_stmt = $db_conn->prepare($sql_query);
$sql_result = $sql_stmt->execute() or die $DBI::errstr;

my $frac_chemical_count = 0; 

# open outuput file
my $db_frac_filename = './output/database_frac_data.tsv'; 
open (DATABASE_FRAC_FILE, "> $db_frac_filename") or die "Error opening output file\n";

print DATABASE_FRAC_FILE "chem_id\tchebi_id\tcas_registry\tfrac_id(FK)\tfrac_id\tfrac_code\tmoa_code\tmoa_name\ttarget_code\ttarget_site\tgroup_name\tchemical_group\tcommon_name\tresistance_risk\tcomments\n";

while (my @row = $sql_stmt->fetchrow_array()) {
  $frac_chemical_count++;
  my $col_count = 0;
  foreach my $column (@row) {
     $col_count++;
     print DATABASE_FRAC_FILE "$column\t" if defined $column; 
  }
  print DATABASE_FRAC_FILE "\n";
}

close (DATABASE_FRAC_FILE);

$sql_stmt->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Total FRAC chemicals:$frac_chemical_count\n";
print "Output file of all FRAC data: $all_data_filename\n";
print "Tab-separated file of all chemical data inserted into chemical and frac tables: $db_frac_filename\n";


