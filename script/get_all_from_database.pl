#!/usr/bin/perl
use strict;
use warnings;
use DBI; # load perl postgresql module
use LWP::Simple;
use XML::Twig;

# load PHI-base functions
use phibase_subroutines 
  qw(connect_to_phibase 
     query_uniprot 
     ontology_mapping
    );

my $db_conn = connect_to_phibase(); # connect to PHI-base database

# open output file
my $db_data_filename = '../output/all_database_data.tsv';  
open (DATABASE_DATA_FILE, "> $db_data_filename") or die "Error opening output file\n";

# print the headers for the output file
print DATABASE_DATA_FILE 
"New PHI-base Acc\tOld PHI-base Acc\tUniProt Acc\tGene Name (PHI-base)\tGene Names (UniProt)\tPathogen Taxon\tDisease\tHost Taxon\tGO Annotations\tPhenotype Outcome\tDefects\tInducers\tCAS Registry IDs\tChEBI IDs\tFRAC Codes\tExperiment Specifications\tCurators\tSpecies Experts\tPubMed IDs\tCuration Date\n";

# first, get details of all interactions from the interaction table
my $sql_stmt = qq(SELECT id,phi_base_accession,curation_date FROM interaction);

my $sql_result = $db_conn->prepare($sql_stmt);
$sql_result->execute() or die $DBI::errstr;

my $interaction_count = 0;

# Read in the relevant ontologies
print "Reading ontology files...\n";
my $obo_parser = OBO::Parser::OBOParser->new;
my $exp_spec_ontology = $obo_parser->work("../ontology/phibase/experiment_specification.obo");
my $phen_outcome_ontology = $obo_parser->work("../ontology/phibase/phenotype_outcome.obo");
my $human_disease_ontology = $obo_parser->work("../ontology/Disease/HumanDisease/doid.obo");
my $plant_disease_ontology = $obo_parser->work("../ontology/Disease/PlantDisease/plant_disease_ontology.obo");

print "Parsing PHI-base data...\n";

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

  my $uniprot_acc = shift @row2;
  my $phibase_gene_name = shift @row2;
  my $path_taxon_id = shift @row2;
  my $phenotype_outcome_id = shift @row2;

  # print the UniProt accession and PHI-base gene name
  print DATABASE_DATA_FILE "$uniprot_acc\t$phibase_gene_name\t";

  # get corresponding data from UniProt

  # declare variable to store UniProt fields
  my $uniprot_gene_names = "";

  # RESTful URL query to get gene names for the current UniProt accession
  my $query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:$uniprot_acc&columns=genes";

  # execute query and process response
  my $gene_names_response = query_uniprot($query);
  my @gene_names_plus_header = split ("\n",$gene_names_response); # split into header & gene names
  my $gene_names_string = $gene_names_plus_header[1]; # the gene names string is second element, after the header
  # check if any gene names have been defined
  if (defined $gene_names_string) {
    my @gene_names = split (" ",$gene_names_string); # split into array of individual gene names
    foreach my $gene_name (@gene_names) {
      $uniprot_gene_names .= "$gene_name;";
    }
  }
  # remove the final semi-colon from end of the string
  $uniprot_gene_names =~ s/;$//;
  print DATABASE_DATA_FILE "$uniprot_gene_names\t";


  # get the pathogen taxon details from the ENA web service
  $query = "http://www.ebi.ac.uk/ena/data/view/Taxon:$path_taxon_id&display=xml";
  my $xml_response = get $query or die "Error getting $query";

  # use XML twig to parse the XML data
  my $xml_twig = XML::Twig->new();
  $xml_twig->parse($xml_response);

  # parse the XML data to get the relevant pathogen taxon info
  my $path_taxon = $xml_twig->root->first_child('taxon');
  my $path_taxon_name = $path_taxon->{'att'}->{'scientificName'};

  # print the pathogen taxon details
  print DATABASE_DATA_FILE "$path_taxon_id:$path_taxon_name\t";


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
  # need to check if a disease exists
  if (defined $disease_id) {

    # use the disease ontologies to retrieve the term name, based on the identifier
    # however, we need to find out which ontology it belongs to (human disease or plant disease)
    my $disease_term;
    my $disease_name;

    # first try to get name from plant disease ontology
    $disease_term = $plant_disease_ontology->get_term_by_id($disease_id);
 
    if (defined $disease_term) {
      # if found, then look up name
      $disease_name = $disease_term->name;
    } else { # if not defined, then look up human disease ontology
      $disease_term = $human_disease_ontology->get_term_by_id($disease_id);
      $disease_name = $disease_term->name;
    }

    # now print the disease id and name 
    print DATABASE_DATA_FILE "$disease_id:$disease_name\t";

  } else { # no disease found
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

  # get the host taxon details from the ENA web service
  $query = "http://www.ebi.ac.uk/ena/data/view/Taxon:$host_taxon_id&display=xml";
  $xml_response = get $query or die "Error getting $query";

  # use XML twig to parse the XML data
  $xml_twig = XML::Twig->new();
  $xml_twig->parse($xml_response);

  # parse the XML data to get the relevant host taxon info
  my $host_taxon = $xml_twig->root->first_child('taxon');
  my $host_taxon_name = $host_taxon->{'att'}->{'scientificName'};

  # print the host taxon details
  print DATABASE_DATA_FILE "$host_taxon_id:$host_taxon_name\t";


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
    my $go_term = "";

    # retrieve the name of the GO term, using the Quick REST web service
    my $query = "http://www.ebi.ac.uk/QuickGO/GTerm?id=$go_id&format=oboxml";
    my $xml_response = get $query;

    # use XML twig to parse the XML data
    my $xml_twig = XML::Twig->new();

    if (defined $xml_response) {
       # parse the XML data to get the GO term name
       $xml_twig->parse($xml_response);
       $go_term = $xml_twig->root->first_child('term')->field('name');
    } else {
       print STDERR "ERROR: Gene Ontology term not found for $go_id\n";
    }


    if (defined $go_evid_code) {  # GO term with evid code
      $go_output_string .= "$go_id($go_evid_code):$go_term;";
    } else {  # GO term without evid code
      $go_output_string .= "$go_id:$go_term;";
    }

  }

  # remove the final semi-colon from end of the string
  $go_output_string =~ s/;$//;
  # print the list of GO terms to file
  print DATABASE_DATA_FILE "$go_output_string\t";


  # output the Phenotype Outcome, which has already been retrieved
  # for comparison with PHI-base 3.6 spreadsheet, it is output in this position 
  # use the ontology to retrieve the term name, based on the identifier
  my $phen_outcome_name = $phen_outcome_ontology->get_term_by_id($phenotype_outcome_id)->name;
  print DATABASE_DATA_FILE "$phenotype_outcome_id:$phen_outcome_name\t";


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


  # get the CAS IDs
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

  # initalise output string for CAS IDs
  my $cas_output_string = "";

  # since there may be multiple anti-infectives,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $cas_registry = shift @row2; 
    $cas_output_string .= "$cas_registry;" if defined $cas_registry;
  }

  # remove the final semi-colon from end of the strings
  $cas_output_string =~ s/;$//;
  # print the list of CAS IDs
  print DATABASE_DATA_FILE "$cas_output_string\t";


  # get the ChEBI ID IDs
  $sql_stmt2 = qq(SELECT chebi_id
                    FROM interaction,
                         interaction_inducer_chemical,
                         chemical
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_inducer_chemical.interaction_id
                     AND chemical.id = interaction_inducer_chemical.chemical_id
                ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for ChEBI IDs
  my $chebi_output_string = "";

  # since there may be multiple ChEBI IDs,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $chebi_id = shift @row2;
    $chebi_output_string .= "$chebi_id;" if defined $chebi_id;
  }

  # remove the final semi-colon from end of the strings
  $chebi_output_string =~ s/;$//;
  # print the list of inducers to file
  print DATABASE_DATA_FILE "$chebi_output_string\t";


  # get the FRAC related fields
  $sql_stmt2 = qq(SELECT frac_code
                    FROM interaction,
                         interaction_anti_infective_chemical,
                         chemical,
                         frac
                   WHERE interaction.id = $interaction_id
                     AND interaction.id = interaction_anti_infective_chemical.interaction_id
                     AND chemical.id = interaction_anti_infective_chemical.chemical_id
                     AND frac.id = chemical.frac_id
                ;);

  $sql_result2 = $db_conn->prepare($sql_stmt2);
  $sql_result2->execute() or die $DBI::errstr;

  # initalise output string for FRAC codes
  my $frac_output_string = "";

  # since there may be multiple FRAC codes,
  # need to retrieve all of them and construct output string 
  # based on semi-colon delimiter
  while (@row2 = $sql_result2->fetchrow_array()) {
    my $frac_code = shift @row2;
    $frac_output_string .= "$frac_code;" if defined $frac_code;
  }

  # remove the final semi-colon from end of the strings
  $frac_output_string =~ s/;$//;
  # print the list of inducers to file
  print DATABASE_DATA_FILE "$frac_output_string\t";


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
    my $exp_spec_id = shift @row2;
    
    # get the term name from the ontology, based on the identifier
    my $exp_spec_name = $exp_spec_ontology->get_term_by_id($exp_spec_id)->name;

    $exp_spec_output_string .= "$exp_spec_id:$exp_spec_name;";
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

  # print message for every 50th PHI-base interaction processed
  print "PHI-base annotations processed:$interaction_count\n" unless ($interaction_count % 50);

}

close (DATABASE_DATA_FILE);

$sql_result->finish() or die "Failed to finish SQL statement\n";
$db_conn->disconnect() or die "Failed to disconnect database\n";

print "\nProcess completed successfully.\n";
print "Total interactions:$interaction_count\n";
print "Tab-separated file of all PHI-base data: $db_data_filename\n\n";

