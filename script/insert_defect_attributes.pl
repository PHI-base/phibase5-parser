#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# hard-code insert of valid defect attributes and their associated valid defect values
my $sql_statement = qq(INSERT INTO defect_attribute (attribute) VALUES ('Mating Defect'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Pre-penetration'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Penetration'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Post-penetration'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Vegetative Spores'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Sexual Spores'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('In Vitro Growth'); 
		     INSERT INTO defect_attribute (attribute) VALUES ('Spore Germination'); 
		     INSERT INTO defect_value (value) VALUES ('Yes');
		     INSERT INTO defect_value (value) VALUES ('No'); 
		     INSERT INTO defect_value (value) VALUES ('Increased'); 
		     INSERT INTO defect_value (value) VALUES ('Wild Type'); 
		     INSERT INTO defect_value (value) VALUES ('Reduced'); 
		     INSERT INTO defect_value (value) VALUES ('Defective'); 
		     INSERT INTO defect_value (value) VALUES ('Aberrant'); 
		     INSERT INTO defect_value (value) VALUES ('Unaffected'); 
		     INSERT INTO defect_value (value) VALUES ('No Data'); 
		    );

my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Defect attributes entered into database:8\n";
print "Defect values entered into database:9\n";

