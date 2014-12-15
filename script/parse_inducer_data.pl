#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open the tab separated values (TSV) version of the inducer data with chemical identifiers
my $inducer_tsv_filename = '../mapping/inducers_4pt0_with_chebi.tsv';
open (TSV_FILE, $inducer_tsv_filename) or die "Error opening input file\n";
print "Processing inducer data from $inducer_tsv_filename...\n";

# counters for statistics
my $inducer_count = 0;
my $chebi_and_cas_count = 0;
my $chebi_only_count = 0;
my $cas_only_count = 0;
my $no_chebi_or_cas_count = 0;

# go through each of the remaining lines of the tab-delimited file (each representing a single chemical)
# save the values of the approriate output file, then insert into the chemical database table
while (<TSV_FILE>) {

   chomp;
   $inducer_count++;

   # each value is separated based on the tab, then saved as an element of the array
   my @chemical_data = split(/\t/,$_);

   my $inducer_name = shift @chemical_data;
   my $chebi_id = shift @chemical_data;
   my $chebi_name = shift @chemical_data;
   my $cas_registry = shift @chemical_data;

   # declare SQL statement
   my $sql_stmt;

   # to permit single-quote in inducer name
   # need to escape it by adding a second single quote
   $inducer_name =~ s/'/''/g; 

   # increment counters for identifiers that are available
   if (defined $chebi_id and defined $cas_registry) {
     $chebi_and_cas_count++;
     # create insert statement to insert chemical
     # with both ChEBI ID and CAS registry
     # unless this chemical already exists
     $sql_stmt = qq(INSERT INTO chemical (name,chebi_id,cas_registry) 
                      SELECT '$inducer_name','$chebi_id','$cas_registry'
                      WHERE NOT EXISTS (
                        SELECT 1 FROM chemical
                        WHERE chebi_id = '$chebi_id'
                      )
                      AND NOT EXISTS (
                        SELECT 1 FROM chemical
                        WHERE cas_registry = '$cas_registry'
                      );
                   );
   } elsif (defined $chebi_id) {
     $chebi_only_count++;
     # create insert statement to insert chemical
     # with only the ChEBI ID
     # unless this chemical already exists
     $sql_stmt = qq(INSERT INTO chemical (name,chebi_id) 
                      SELECT '$inducer_name','$chebi_id'
                      WHERE NOT EXISTS (
                        SELECT 1 FROM chemical
                        WHERE chebi_id = '$chebi_id'
                      );
                   );
   } elsif (defined $cas_registry) {
     $cas_only_count++;
     # create insert statement to insert chemical
     # with only the CAS registry number
     # unless this chemical already exists
     $sql_stmt = qq(INSERT INTO chemical (name,cas_registry) 
                      SELECT '$inducer_name','$cas_registry'
                      WHERE NOT EXISTS (
                        SELECT 1 FROM chemical
                        WHERE cas_registry = '$cas_registry'
                      );
                   );
   } else {  # neither ChEBI nor CAS available
     $no_chebi_or_cas_count++;
     # create insert statement to insert chemical
     # without either the ChEBI ID or the CAS registry number
     $sql_stmt = qq(INSERT INTO chemical (name) 
                      VALUES ('$inducer_name');
                   );
   }
 
   my $sql_result = $db_conn->prepare($sql_stmt);
   $sql_result->execute() or die $DBI::errstr;

} # end of inducer data file

close (TSV_FILE);

$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Total inducer chemicals:$inducer_count\n";
print "Inducer chemicals with both ChEBI ID and CAS registry number:$chebi_and_cas_count\n";
print "Inducer chemicals with a ChEBI ID but no CAS registry number:$chebi_only_count\n";
print "Inducer chemicals with a CAS registry number but no ChEBI ID:$cas_only_count\n";
print "Inducer chemicals with neither a ChEBI ID nor a CAS registry number:$no_chebi_or_cas_count\n\n";

