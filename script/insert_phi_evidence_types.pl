#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# hard-code insert of valid defect attributes and their associated valid defect values
my $sql_statement = qq(INSERT INTO phi_evidence (phi_evidence) VALUES ('cell growth assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('chromatin immunoprecipitation'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('co-immunoprecipitation'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('count assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('direct injection into host organism'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('electrophoretic mobility shift assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('enzyme assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('epitope-tagged protein immunolocalization'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('functional test in host organism'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('gel electrophoresis'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('immunolocalization'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('microscopy'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('particle size'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('plasmid maintenance assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('quantitative PCR'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('reporter gene assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('sequence analysis'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('substance quantification'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('transcript expression level'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('transient expression in host organism'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('western blot assay'); 
                       INSERT INTO phi_evidence (phi_evidence) VALUES ('other'); 
		    );

my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

$db_conn->disconnect() or die "Failed to disconnect database\n";

print "PHI evidence types successfully entered into database\n";

