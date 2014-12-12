#!/usr/bin/perl
use strict;
use warnings;

use phibase_subroutines qw(print_ontology_terms); # load PHI-base function

print_ontology_terms('phenotype_outcome','../ontology/phibase/phenotype_outcome.obo');

