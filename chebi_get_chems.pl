#!/usr/bin/perl -w
# SOAP::Lite version 0.67
# Please note: ChEBI webservices uses document/literal binding

#use SOAP::Lite + trace => qw(debug);
use SOAP::Lite;

# open output file
my $chebi_filename = './output/chebi_chem_details.txt';
open (CHEBI_FILE, "> $chebi_filename") or die "Error opening output file\n";

print "Printing to ChEBI output file $chebi_filename...\n";

# define the ChEBI ID (temperarily hard-code an example ChEBI ID)
# TODO: ChEBI ID will be a value returned from database query of chemical table
my $chebi_id = 'CHEBI:40909';

# Setup SOAP web service, using the ChEBI WSDL XML as description of available fields
my $WSDL = 'http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl';
my $nameSpace = 'http://www.ebi.ac.uk/webservices/chebi';
my $soap = SOAP::Lite
-> uri($nameSpace)
-> proxy($WSDL);

# Setup method and parameters (in this case, ChEBI ID is the parameter)
my $method = SOAP::Data->name('getCompleteEntity')
             ->attr({xmlns => $nameSpace});
#my $method = SOAP::Data->name('getLiteEntity')
#              ->attr({xmlns => $nameSpace});
my @params = ( SOAP::Data->name(chebiId => $chebi_id));
#my @params = ( SOAP::Data->name(search => '40909'),
#               SOAP::Data->name(searchCategory => 'CHEBI ID'),
#               SOAP::Data->name(maximumResults => '200'),
#               SOAP::Data->name(stars => 'ALL'));


# Call method
my $som = $soap->call($method => @params);

# Retrieve chemical data, based on metadata fields available from the WSDL XML definitions
#@stuff = $som->valueof('//OntologyParents//chebiId');
#@stuff = $som->valueof('//ListElement//chebiAsciiName');
@names = $som->valueof('//chebiAsciiName');
@synonyms = $som->valueof('//Synonyms//data');
@registry = $som->valueof('//RegistryNumbers//data');
@registry_types = $som->valueof('//RegistryNumbers//type');

# print the ChEBI ID, names & synonyms
print CHEBI_FILE "\n$chebi_id\n";
print CHEBI_FILE "Name: $_\n" foreach @names; 
print CHEBI_FILE "Synonym: $_\n" foreach @synonyms;

# only interested in the CAS registry, not other registry numbers
# so test the registry type, then print the corresponding CAS registry number
# (which will have an identical index to the registry type)
my $registry_index = 0;
foreach $registry_type (@registry_types) {
  if ($registry_type eq 'CAS Registry Number') {
    print CHEBI_FILE "CAS Registry: $registry[$registry_index]\n";
  }
  $registry_index++;
}

close (CHEBI_FILE);

