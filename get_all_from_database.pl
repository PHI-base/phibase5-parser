#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module

use phibase_subroutines qw(connect_to_phibase query_uniprot ontology_mapping); # load PHI-base functions

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open output file
my $db_data_filename = './output/all_database_data.tsv';  
open (DATABASE_DATA_FILE, "> $db_data_filename") or die "Error opening output file\n";

# print the headers for the output file
print DATABASE_DATA_FILE 
"New PHI-base Acc\tOld PHI-base Acc\tUniProt Acc\tGene Name\tPathogen Taxon ID\tDisease ID\tHost Taxon ID\tGO Annotations\tPhenotype Outcome ID\tDefects\tInducers\tCAS Registry IDs\tExperiment Spec ID\tCurators\tSpecies Experts\tPubMed IDs\tCuration Date\n";

# retrieve all of the data from the relevant tables of phibase
# and save in a tab-delimited file, with appropriate headings

# hash to store all annotations with a valid PHI-base accession
#my %valid_phibase_data;

# hash to store all values of the current annotation
#my %phi_base_annotation;

# first, get details of all interactions from the interaction table
my $sql_stmt = qq(SELECT id,phi_base_accession,curation_date FROM interaction);

my $sql_result = $db_conn->prepare($sql_stmt);
$sql_result->execute() or die $DBI::errstr;

my $interaction_count = 0;

while (my @row = $sql_result->fetchrow_array()) {

  $interaction_count++;

  my $interaction_id = shift @row;
  my $phibase_accession = shift @row;
  my $curation_date = shift @row;

  print DATABASE_DATA_FILE "$phibase_accession\t";

  # get the obsolete PHI-base accession
  my $sql_stmt2 = qq(SELECT obsolete_accession FROM obsolete
                       WHERE phi_base_accession = '$phibase_accession';
                    );

  my $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;
  my @row2 = $sql_result2->fetchrow_array();
  my $obsolete_accession = shift @row2;

  print DATABASE_DATA_FILE "$obsolete_accession\t";

  # get the pathogen gene related fields 
  $sql_stmt2 = qq(SELECT uniprot_accession,
                         gene_name,
                         pathogen_gene.ncbi_taxon_id,
                         phenotype_outcome_id
                    FROM interaction,
                         interaction_pathogen_gene_mutant, 
                         pathogen_gene_mutant,
                         pathogen_gene
                    WHERE interaction.id = $interaction_id
                      AND interaction.id = interaction_pathogen_gene_mutant.interaction_id
                      AND pathogen_gene_mutant.id = interaction_pathogen_gene_mutant.pathogen_gene_mutant_id
                      AND pathogen_gene.id = pathogen_gene_mutant.pathogen_gene_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;
  @row2 = $sql_result2->fetchrow_array();

  my $uniprot_accession = shift @row2;
  my $gene_name = shift @row2;
  my $path_taxon_id = shift @row2;
  my $phenotype_outcome = shift @row2;

  print DATABASE_DATA_FILE "$uniprot_accession\t$gene_name\t$path_taxon_id\t";


  # get the disease related fields 
  $sql_stmt2 = qq(SELECT disease_id
                    FROM interaction,
                         interaction_disease
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_disease.interaction_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;
  @row2 = $sql_result2->fetchrow_array();

  my $disease_id = shift @row2;

  # since this field is not mandatory
  # need to check if a result exists
  if (defined $disease_id) {
    print DATABASE_DATA_FILE "$disease_id\t";
  } else {
    print DATABASE_DATA_FILE "\t";
  }


  # get the host related fields 
  $sql_stmt2 = qq(SELECT interaction_host.ncbi_taxon_id
                    FROM interaction,
                         interaction_host
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_host.interaction_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;
  @row2 = $sql_result2->fetchrow_array();

  my $host_taxon_id = shift @row2;

  print DATABASE_DATA_FILE "$host_taxon_id\t";


  # get the Gene Ontology annotation fields 
  $sql_stmt2 = qq(SELECT go_id,
                         go_evidence_code
                    FROM interaction,
                         interaction_go_annotation
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_go_annotation.interaction_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for GO terms
  my $go_output_string = "";

  # since there may be multiple GO terms,
  # need to retrieve all of them and construct output string 
  # based on comma and semi-colon delimiters
  while (@row2 = $sql_result2->fetchrow_array()) {

    my $go_id = shift @row2;
    my $go_evid_code = shift @row2;

    if (defined $go_evid_code) {  # GO term with evid code
      $go_output_string .= "$go_id,$go_evid_code;";
    } else {  # GO term without evid code
      $go_output_string .= "$go_id;";
    }

  }

  # remove the final semi-colon from end of the string
  $go_output_string =~ s/;$//;
  # print the list of GO terms to file
  print DATABASE_DATA_FILE "$go_output_string\t";


  # output the Phenotype Outcome, which has already been retrieved
  # for comparison with PHI-base 3.6, it is output here 
  print DATABASE_DATA_FILE "$phenotype_outcome\t";


  # get the Defect fields 
  $sql_stmt2 = qq(SELECT defect_attribute.attribute,
                         defect_value.value
                    FROM interaction,
                         interaction_defect,
                         defect_attribute,
                         defect_value
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_defect.interaction_id
                     AND defect_attribute.id = interaction_defect.defect_attribute_id
                     AND defect_value.id = interaction_defect.defect_value_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for defects
  my $defect_output_string = "";

  # since there may be multiple defects,
  # need to retrieve all of them and construct output string 
  # based on colon and semi-colon delimiters
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $attribute = shift @row2;
    my $value = shift @row2;
    $defect_output_string .= "$attribute:$value;";
  }

  # remove the final semi-colon from end of the string
  $defect_output_string =~ s/;$//;
  # print the list of GO terms to file
  print DATABASE_DATA_FILE "$defect_output_string\t";


  # get the Inducer fields 
  $sql_stmt2 = qq(SELECT chemical.name
                    FROM interaction,
                         interaction_inducer_chemical,
                         chemical
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_inducer_chemical.interaction_id
                     AND chemical.id = interaction_inducer_chemical.chemical_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for Inducers
  my $inducer_output_string = "";

  # since there may be multiple inducers,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $chemical = shift @row2;
    $inducer_output_string .= "$chemical;";
  }

  # remove the final semi-colon from end of the string
  $inducer_output_string =~ s/;$//;
  # print the list of inducers to file
  print DATABASE_DATA_FILE "$inducer_output_string\t";


  # get the Anti-infective fields 
  $sql_stmt2 = qq(SELECT cas_registry
                    FROM interaction,
                         interaction_anti_infective_chemical,
                         chemical
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_anti_infective_chemical.interaction_id
                     AND chemical.id = interaction_anti_infective_chemical.chemical_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for Anti-infectives
  my $anti_infective_output_string = "";

  # since there may be multiple anti-infectives,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $cas_registry = shift @row2;
    $anti_infective_output_string .= "$cas_registry;";
  }

  # remove the final semi-colon from end of the string
  $anti_infective_output_string =~ s/;$//;
  # print the list of inducers to file
  print DATABASE_DATA_FILE "$anti_infective_output_string\t";


  # get the Experiment Specification fields 
  $sql_stmt2 = qq(SELECT experiment_spec_id
                    FROM interaction,
                         interaction_experiment_spec
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_experiment_spec.interaction_id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for exp spec
  my $exp_spec_output_string = "";

  # since there may be multiple experiment specifications
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $exp_spec = shift @row2;
    $exp_spec_output_string .= "$exp_spec;";
  }

  # remove the final semi-colon from end of the string
  $exp_spec_output_string =~ s/;$//;
  # print the list of inducers to file
  print DATABASE_DATA_FILE "$exp_spec_output_string\t";


  # get the Curator fields 
  $sql_stmt2 = qq(SELECT curator.id,
                         curator.name
                    FROM interaction,
                         interaction_curator,
                         curator
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_curator.interaction_id
                     AND interaction_curator.curator_id = curator.id
                 ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for both Curators and Species Experts
  my $curator_output_string = "";
  my $species_experts_string = "";

  # since there may be multiple curators,
  # need to retrieve all of them and construct output string 
  # based on comma and semi-colon delimiters
  while (@row2 = $sql_result2->fetchrow_array()) {

    my $curator_id = shift @row2;
    my $curator_name = shift @row2;

    # need to determine if the curator belongs to
    # a known organisation
    my $sql_stmt3 = qq(SELECT curation_organisation.name
                         FROM curator,
                              curation_organisation
                        WHERE curator.id = $curator_id
                          AND curation_organisation.id = curator.curation_organisation_id
                     ;);

    my $sql_result3 = $db_conn->prepare($sql_stmt3);
    $sql_result3->execute() or die $DBI::errstr;
    my @row3 = $sql_result3->fetchrow_array();
    my $organisation = shift @row3;

    if (defined $organisation) {  # curator with organisation
      $curator_output_string .= "$curator_name,$organisation;";
    } else {  # curator without organisation
      $curator_output_string .= "$curator_name;";
    }

    # need to determine if the curator is a species expert
    # based on the taxon id of the pathogen
    my $sql_stmt4 = qq(SELECT curator_id
                         FROM species_expert,
                              curator
                        WHERE species_expert.curator_id = $curator_id
                          AND species_expert.ncbi_taxon_id = $path_taxon_id
                     ;);

    my $sql_result4 = $db_conn->prepare($sql_stmt4);
    $sql_result4->execute() or die $DBI::errstr;
    my @row4 = $sql_result4->fetchrow_array();
    my $expert_curator_id = shift @row4;

    # if a curator id was returned,
    # then the curator is a species expert
    # so add their name to the experts list
    if (defined $expert_curator_id) {
      $species_experts_string .= "$curator_name;";
    }

  }

  # remove the final semi-colon from end of the strings
  $curator_output_string =~ s/;$//;
  $species_experts_string =~ s/;$//;
  # print the list of curators and species experts to file
  print DATABASE_DATA_FILE "$curator_output_string\t$species_experts_string\t";


  # get the literature fields 
  $sql_stmt2 = qq(SELECT interaction_literature.pubmed_id
                    FROM interaction,
                         interaction_literature
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_literature.interaction_id
                ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for literature
  my $pubmed_output_string = "";

  # since there may be multiple PubMed articles,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $pubmed_id = shift @row2;
    $pubmed_output_string .= "$pubmed_id;";
  }

  # remove the final semi-colon from end of the string
  $pubmed_output_string =~ s/;$//;
  # print the list of pubmed ids to file
  print DATABASE_DATA_FILE "$pubmed_output_string\t";

  # finally, output curation date
  print DATABASE_DATA_FILE "$curation_date\n";

}

close (DATABASE_DATA_FILE);

$sql_result->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "\nProcess completed successfully.\n";
print "Total interactions:$interaction_count\n";
print "Tab-separated file of all PHI-base data: $db_data_filename\n\n";

