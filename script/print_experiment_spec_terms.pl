#!/usr/bin/perl
use strict;
use warnings;

use phibase_subroutines qw(print_ontology_terms); # load PHI-base function

print_ontology_terms('experiment_spec','../ontology/phibase/experiment_specification.obo');

