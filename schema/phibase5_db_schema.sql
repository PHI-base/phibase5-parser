CREATE TABLE interaction (
    id SERIAL PRIMARY KEY,
    phi_base_accession varchar(50),
    curation_date date
);

CREATE TABLE interaction_host (
    id SERIAL PRIMARY KEY,
    interaction_id integer REFERENCES interaction,
    ncbi_taxon_id integer,
    host_strain_name varchar(50),
    first_target_uniprot_accession varchar(50),
    genbank_locus_id varchar(50)
);

CREATE TABLE interaction_phi_host_phenotype (
    interaction_host_id integer REFERENCES interaction_host,
    phi_phenotype_id varchar(50),
    phi_evidence varchar(50),
    PRIMARY KEY (interaction_host_id, phi_phenotype_id)
);

CREATE TABLE interaction_host_tissue (
    interaction_host_id integer REFERENCES interaction_host,
    brenda_tissue_id varchar(50),
    PRIMARY KEY (interaction_host_id, brenda_tissue_id)
);

CREATE TABLE usda (
    id SERIAL PRIMARY KEY,
    pathogen_taxon_id integer,
    host_kingdom varchar(50),
    num_plant_hosts varchar(50),
    num_interactions_usda integer,
    num_interactions_plantwise integer
);

CREATE TABLE usda_pathogen_lifestyle (
    usda_id integer REFERENCES usda,
    pathogen_lifestyle_id varchar(50),
    PRIMARY KEY (usda_id, pathogen_lifestyle_id)
);

CREATE TABLE usda_natural_host (
    usda_id integer REFERENCES usda,
    natural_host_id varchar(50),
    PRIMARY KEY (usda_id, natural_host_id)
);

CREATE TABLE pathogen_gene (
    id SERIAL PRIMARY KEY,
    ncbi_species_taxon_id integer,
    ncbi_taxon_id integer,
    pathogen_strain_name varchar(50),
    gene_name varchar(100),
    uniprot_accession varchar(50),
    uniparc_id varchar(50),
    gene_locus_id varchar(50),
    gene_locus_id_type varchar(130)
);

CREATE TABLE go_evidence (
    code varchar(50) PRIMARY KEY,
    category varchar(50),
    description varchar (50)
);

CREATE TABLE pathogen_gene_go_annotation (
    id SERIAL PRIMARY KEY,
    pathogen_gene_id integer REFERENCES pathogen_gene,
    pubmed_id varchar(50),
    go_id varchar(50),
    go_evidence_code varchar(50) REFERENCES go_evidence
    --PRIMARY KEY (pathogen_gene_id, pubmed_id, go_id)
);

CREATE TABLE pathogen_gene_go_annot_ext (
    id SERIAL PRIMARY KEY,
    pathogen_gene_go_annotation_id integer REFERENCES pathogen_gene_go_annotation,
    go_annot_ext_relation varchar(50),
    go_annot_ext_value varchar(50)
);

CREATE TABLE effector_gene (
    id SERIAL PRIMARY KEY,
    pathogen_gene_id integer REFERENCES pathogen_gene,
    pubmed_id varchar(50),
    phi_effector_id varchar(50),
    phi_effector_evidence_code varchar(50),
    location_in_host_go_id varchar(50),
    host_target_uniprot_acc varchar(50)
);

CREATE TABLE pathogen_gene_allele (
    id SERIAL PRIMARY KEY,
    pathogen_gene_id integer REFERENCES pathogen_gene,
    allele_name varchar(100),
    allele_type varchar(50),
    allele_description varchar(50)
);

CREATE TABLE interaction_pathogen_gene_allele (
    interaction_id integer REFERENCES interaction,
    pathogen_gene_allele_id integer REFERENCES pathogen_gene_allele,
    allele_expression varchar(50),
    PRIMARY KEY (interaction_id, pathogen_gene_allele_id)
);

CREATE TABLE interaction_phi_interaction_phenotype (
    interaction_id integer REFERENCES interaction,
    phi_phenotype_id varchar(50),
    phi_evidence varchar(50),
    PRIMARY KEY (interaction_id, phi_phenotype_id)
);

CREATE TABLE obsolete (
    id SERIAL PRIMARY KEY,
    phi_base_accession varchar(50),
    obsolete_accession varchar(50)
);

CREATE TABLE interaction_literature (
    interaction_id integer REFERENCES interaction,
    pubmed_id varchar(50),
    doi varchar(50)
);

CREATE TABLE pathogen_interacting_protein (
    interaction_id integer REFERENCES interaction,
    uniprot_accession varchar(50),
    PRIMARY KEY (interaction_id, uniprot_accession)
);

CREATE TABLE gene_post_trans_mod (
    id SERIAL PRIMARY KEY,
    pathogen_gene_id integer REFERENCES pathogen_gene,
    pubmed_id varchar(50),
    psi_mod_id varchar(50),
    psi_mod_evid_code varchar(50)
);

CREATE TABLE post_trans_mod_annot_ext (
    id SERIAL PRIMARY KEY,
    gene_post_trans_mod_id integer REFERENCES gene_post_trans_mod,
    post_trans_mod_annot_ext_relation varchar(50),
    post_trans_mod_annot_ext_value varchar(50)
);

CREATE TABLE interaction_transient_assay (
    interaction_id integer REFERENCES interaction,
    bioassay_ontology_id integer,
    PRIMARY KEY (interaction_id, bioassay_ontology_id)    
);

CREATE TABLE interaction_experiment_spec (
    interaction_id integer REFERENCES interaction,
    experiment_spec_id varchar(50),
    PRIMARY KEY (interaction_id, experiment_spec_id)    
);

CREATE TABLE curation_organisation (
    id SERIAL PRIMARY KEY,
    name varchar(50)    
);

CREATE TABLE curator (
    id SERIAL PRIMARY KEY,
    initials varchar(50),
    name varchar(50),
    email varchar(50),
    curation_organisation_id integer REFERENCES curation_organisation
);

CREATE TABLE interaction_curator (
    interaction_id integer REFERENCES interaction,
    curator_id integer REFERENCES curator,
    PRIMARY KEY (interaction_id, curator_id)    
);

CREATE TABLE interaction_approver (
    interaction_id integer REFERENCES interaction,
    curator_id integer REFERENCES curator,
    PRIMARY KEY (interaction_id, curator_id)    
);

CREATE TABLE species_expert (
    ncbi_taxon_id integer,
    curator_id integer REFERENCES curator,
    PRIMARY KEY (ncbi_taxon_id, curator_id)
);

CREATE TABLE frac (
    id SERIAL PRIMARY KEY,
    frac_code varchar(50),    
    moa_code varchar(50),
    moa_name varchar(100),
    target_code varchar(50),    
    target_site varchar(100),    
    group_name varchar(100),    
    chemical_group varchar(100),
    common_name varchar(100),
    resistance_risk varchar(50),
    comments varchar(1000)
);

CREATE TABLE chemical (
    id SERIAL PRIMARY KEY,
    name varchar(100),
    chebi_id varchar(50),
    cas_registry varchar(50),    
    frac_id integer REFERENCES frac,
    mode_in_planta varchar(50)
);

CREATE TABLE interaction_anti_infective_chemical (
    interaction_id integer REFERENCES interaction,
    chemical_id integer REFERENCES chemical,
    PRIMARY KEY (interaction_id, chemical_id)    
);

CREATE TABLE interaction_inducer_chemical (
    interaction_id integer REFERENCES interaction,
    chemical_id integer REFERENCES chemical,
    PRIMARY KEY (interaction_id, chemical_id)    
);

CREATE TABLE interaction_inducer_gene (
    interaction_id integer REFERENCES interaction,
    uniprot_accession varchar(50),
    PRIMARY KEY (interaction_id, uniprot_accession)    
);

CREATE TABLE disease_severity (
    id SERIAL PRIMARY KEY,
    severity varchar(50)    
);

CREATE TABLE interaction_disease (
    interaction_id integer REFERENCES interaction,
    disease_id varchar(50),
    disease_severity_id integer REFERENCES disease_severity,
    PRIMARY KEY (interaction_id, disease_id)       
);

CREATE TABLE interaction_phi_pathogen_phenotype (
    interaction_id integer REFERENCES interaction,
    phi_phenotype_id varchar(50),
    phi_evidence varchar(50),
    PRIMARY KEY (interaction_id, phi_phenotype_id)
);

CREATE TABLE interaction_disease_formation (
    interaction_id integer REFERENCES interaction,
    disease_formation_id varchar(50),
    PRIMARY KEY (interaction_id, disease_formation_id)
);

