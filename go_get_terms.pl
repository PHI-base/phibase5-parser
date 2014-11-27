#!/usr/bin/perl
use strict;
use warnings;

use LWP::UserAgent;
use XML::XPath;
use XML::XPath::XMLParser;

my $ua = LWP::UserAgent->new;
my $req = HTTP::Request->new(GET => 'http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0003824&format=oboxml');
my $res = $ua->request($req, "term.xml");

my $xpath = XML::XPath->new(filename => "term.xml");

my $nodeset = $xpath->find("/obo/term/name");
foreach my $node ($nodeset->get_nodelist) {
  print XML::XPath::XMLParser::as_string($node) . "\n";
}

exit;
