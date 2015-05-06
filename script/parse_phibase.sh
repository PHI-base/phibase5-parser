#!/bin/sh

#############################################################
#
# Parent shell script to re-create PHI-base 5 database
# and insert the required data from a previous version,
# based on data from the PHI-base spreadsheet
#
# WARNING: All existing data for PHI-base 5 will be deleted!
#
##############################################################

echo "\nStarted process to Re-create PHI-base 5 database..."

# drop the current database, then recreate it
dropdb -U postgres phibase
createdb -U postgres phibase

# insert tables from SQL schema, with appropriate restraints 
psql -U postgres phibase < ../schema/phibase5_db_schema.sql
echo "\nGenerated SQL schema successfully"

echo "\nInserting valid defect attributes and values..."
perl insert_defect_attributes.pl

echo "\nInserting Gene Ontology Evidence Codes..."
perl insert_go_evidence_codes.pl

echo "\nParsing Curator data..."
perl parse_curator_data.pl

echo "\nParsing FRAC fungicide chemical data..."
perl parse_frac_data.pl

echo "\nParsing Inducer chemical data..."
perl parse_inducer_data.pl

echo "\nParsing USDA data..."
perl parse_usda_data.pl

echo "\nPHI-base 5 database populated with all pre-requisite data.\n"
