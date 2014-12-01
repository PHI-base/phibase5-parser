#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# hard-code insert of valid Gene Ontology Evidence Codes
#
# for more information, see 
# http://geneontology.org/page/guide-go-evidence-codes
# and the GO Evidence Code Decision Tree at 
# http://geneontology.org/page/evidence-code-decision-tree
#
# the list below are separated based on the following categories:
# 1. Experimental evidence codes
# 2. Computational evidece codes
# 3. Author statement evidence codes
# 4. Curational statement evidence codes
# 5. Automatically-assigned evidence codes (may not be used in PHI-base)
#
my $sql_statement = qq(INSERT INTO go_evidence (code,category,description) VALUES ('EXP','Experimental','Inferred from experiment');
                       INSERT INTO go_evidence (code,category,description) VALUES ('IMP','Experimental','Inferred from mutant phenotype'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IGI','Experimental','Inferred from genetic interaction'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IPI','Experimental','Inferred from physical interaction'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IDA','Experimental','Inferred from direct assay'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IEP','Experimental','Inferred from expression pattern');
 
                       INSERT INTO go_evidence (code,category,description) VALUES ('ISS','Computational','Inferred from sequence or structure similarity'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('ISA','Computational','Inferred from sequence alignment'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('ISO','Computational','Inferred from sequence orthology'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('ISM','Computational','Inferred from sequence model'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IGC','Computational','Inferred from genomic context'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IBA','Computational','Inferred from biological aspect of ancestor'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IBD','Computational','Inferred from biological aspect of descendant'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IKR','Computational','Inferred from key residues'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IRD','Computational','Inferred from rapid divergence'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('RCA','Computational','Inferred from reviewed computational analysis');

                       INSERT INTO go_evidence (code,category,description) VALUES ('TAS','Author statement','Traceable author statement'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('NAS','Author statement','Non-traceable author statement');

                       INSERT INTO go_evidence (code,category,description) VALUES ('IC','Curatorial statement','Inferred by curator'); 
                       INSERT INTO go_evidence (code,category,description) VALUES ('ND','Curatorial statement','No biological data available');
 
                       INSERT INTO go_evidence (code,category,description) VALUES ('IEA','Automatically-assigned','Inferred from electronic annotation'); 
		    );

my $sql_result = $db_conn->do($sql_statement) or die $DBI::errstr;

$db_conn->disconnect() or die "Failed to disconnect database\n";

print "Process completed successfully.\n";
print "Gene Ontology Evidence Codes entered into database:13\n";

