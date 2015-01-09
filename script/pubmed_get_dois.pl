#!/usr/bin/perl
use strict;
use warnings;

use LWP::Simple;
use JSON;
use DBI;

# Attempt to get PubMed IDs based on Document Object Identifiers (DOIs) from the PHI-base spreadsheet

# open the the DOI file from the output columns data
my $doi_filename = '../output/column/spreadsheet_column_doi.txt';
open (DOI_FILE, $doi_filename) || die "Error opening input file\n";
print "Processing DOI data from $doi_filename...\n";

# counters to gather statistics
my $doi_count = 0;
my $doi_found_count = 0;
my $doi_not_found_count = 0;

# open output files
my $pubmed_filename = '../output/pubmed_from_dois.tsv';
my $no_pubmed_filename = '../error/pubmed_from_dois_not_found.tsv';
open (PUBMED_FILE,"> $pubmed_filename") or die "Error opening output file\n";
open (NO_PUBMED_FILE,"> $no_pubmed_filename") or die "Error opening output file\n";
print "Printing PubMed IDs from DOIs to output file $pubmed_filename...\n";

# iterate through each DOI to get the corresponding PubMed ID
while (<DOI_FILE>) {

  # each DOI is separated based on the semi-colon, then saved as an element of the array
  my @doi_data = split(/;/,$_);

  # iterate through each DOI of the current annotation
  foreach my $doi (@doi_data) {

     $doi =~ s/^\s+//; # remove blank space from start of DOI
     $doi =~ s/\s+$//; # remove blank space from end of DOI

     # if DOI is blank, then ignore and go to next DOI in file
     if ($doi eq "" or $doi eq "no data found") {
       next;
     }

     # for URL query, other characters need to be removed
     my $url_doi = $doi;
     $url_doi =~ s/^doi//ig; # remove all prefix 'doi' (ignoring case)
     $url_doi =~ s/://g; # remove all colons from DOI
     $url_doi =~ s/\s+//g; # remove all blank spaces from DOI (including internal spaces)

     # increment count of DOIs
     $doi_count++;

=pod
     # before using doi in URL query it needs to be URL compatible
     # so certain characters, such as round brackets, need to be substituted
     my $url_doi = $doi;
     # substitue all open brackets for %28
     $url_doi =~ s/\(/%28/;
     # substitue all open brackets for %28
     $url_doi =~ s/\)/%29/;
     my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=DOI:$url_doi&format=json";
=cut
     # run REST query and get JSON response
     my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=DOI:$url_doi&format=json";
     my $json_response = get $url;
     my $text_response = decode_json($json_response);

     # parse the PubMed ID from the JSON text
     my $pubmed_id = $text_response->{'resultList'}{'result'}[0]{'id'};

     # if the pubmed ID is found, then print to file
     if (defined $pubmed_id and $pubmed_id ne "") {
        $doi_found_count++;
        print PUBMED_FILE "$doi\t$pubmed_id\n";
     } else { # article not found
        $doi_not_found_count++;
        print STDERR "ERROR:DOI $doi not found\n";
        print NO_PUBMED_FILE "$doi\n";
     }

  } # end foreach DOI in annotation

} # end while DOIs in file

close (PUBMED_FILE);
close (NO_PUBMED_FILE);

print "Total DOIs:$doi_count\n";
print "PubMed IDs found using DOI: $doi_found_count, output file $pubmed_filename\n";
print "PubMed IDs not found using DOI: $doi_not_found_count, output file $no_pubmed_filename\n";

