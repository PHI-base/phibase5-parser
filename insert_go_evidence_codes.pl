#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# hard-code insert of valid Gene Ontology Evidence Codes
# for more information, see the GO Evidence Code Decision Tree
# at http://geneontology.org/page/evidence-code-decision-tree
my $sql_statement = qq(INSERT INTO go_evidence_code (code) VALUES ('IMP'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IGI'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IPI'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IDA'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IEP'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IEA'); 
                       INSERT INTO go_evidence_code (code) VALUES ('ISS'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IGC'); 
                       INSERT INTO go_evidence_code (code) VALUES ('ICA'); 
                       INSERT INTO go_evidence_code (code) VALUES ('TAS'); 
                       INSERT INTO go_evidence_code (code) VALUES ('NAS'); 
                       INSERT INTO go_evidence_code (code) VALUES ('IC'); 
                       INSERT INTO go_evidence_code (code) VALUES ('ND'); 
		    );

my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Gene Ontology Evidence Codes entered into database:13\n";

