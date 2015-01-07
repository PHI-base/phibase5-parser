#!/bin/sh

################################################################
#
# Shell script to get external data related to PHI-base 5 data
#
################################################################

echo "\nStarted process to retrieve PHI-base-related data from external databases..."

echo "\nGetting ChEBI chemical data..."
perl chebi_get_chems.pl

echo "\nGetting all GO annotations..."
perl go_get_all_annotations.pl

echo "\nGetting GO cellular component annotations..."
perl go_get_component_annotations.pl

echo "\nGetting GO molecular function annotations..."
perl go_get_function_annotations.pl

echo "\nGetting GO biological process annotations..."
perl go_get_process_annotations.pl

echo "\nGetting GO term names..."
perl go_get_term_names.pl

echo "\nGetting PubMed citation details..."
perl pubmed_get_details.pl

echo "\nGetting host taxonomy details..."
perl taxon_get_cotyledons.pl

echo "\nGetting UniProt EMBL IDs..."
perl uniprot_get_embl_ids.pl

echo "\nGetting UniProt Gene Names..."
perl uniprot_get_gene_names.pl

echo "\nGetting UniProt Gene Ontology annotations..."
perl uniprot_get_go.pl

echo "\nGetting UniProt Protein Names..."
perl uniprot_get_protein_names.pl

echo "\nGetting UniProt amino acid sequences..."
perl uniprot_seq_fasta.pl

echo "\nGetting UniProt Pathogen Taxa..."
perl uniprot_taxon.pl

