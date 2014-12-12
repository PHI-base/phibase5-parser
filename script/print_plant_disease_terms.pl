#!/usr/bin/perl
use strict;
use warnings;

use phibase_subroutines qw(print_ontology_terms); # load PHI-base function

print_ontology_terms('plant_disease','../ontology/Disease/PlantDisease/plant_disease_ontology.obo');

