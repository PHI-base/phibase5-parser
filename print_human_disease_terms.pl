#!/usr/bin/perl
use strict;
use warnings;

use phibase_subroutines qw(print_ontology_terms); # load PHI-base function

print_ontology_terms('human_disease','ontologies/Disease/HumanDisease/doid.obo');

