-- Database: "genomeAnnos"

-- DROP DATABASE "genomeAnnos";

CREATE DATABASE "genomeAnnos"
  WITH OWNER = postgres
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'en_US.UTF-8'
       LC_CTYPE = 'en_US.UTF-8'
       CONNECTION LIMIT = -1;

DROP TABLE IF EXISTS chromosome CASCADE;
CREATE TABLE chromosome (
    chromosome_id serial PRIMARY KEY,
    chr_number int,
    chr_name varchar,
    chr_length int
);

DROP TABLE IF EXISTS gene CASCADE;
CREATE TABLE gene (
    gene_id serial PRIMARY KEY,
    gene_identifier varchar(64) UNIQUE,
    gene_name varchar,
    chromosome_id int,
    start_position int,
    end_position int,
    strand boolean,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS transcript CASCADE;
CREATE TABLE transcript (
    transcript_id serial PRIMARY KEY,
    transcript_identifier varchar(64) UNIQUE,
    gene_id int,
    start_position int,
    end_position int,
    FOREIGN KEY (gene_id) references gene
);

DROP TABLE IF EXISTS cds CASCADE;
CREATE TABLE cds (
    cds_id serial PRIMARY KEY,
    transcript_id int,
    start_position int,
    end_position int,
    FOREIGN KEY (transcript_id) references transcript
);

DROP TABLE IF EXISTS exon;
CREATE TABLE exon (
    exon_id serial PRIMARY KEY,
    transcript_id int,
    start_position int,
    end_position int,
    FOREIGN KEY (transcript_id) references transcript
);

DROP TABLE IF EXISTS intron;
CREATE TABLE intron (
    intron_id serial PRIMARY KEY,
    transcript_id int,
    start_position int,
    end_position int,
    FOREIGN KEY (transcript_id) references transcript
);

DROP TABLE IF EXISTS annotation_scope CASCADE;
CREATE TABLE annotation_scope (
    annotation_scope_id serial PRIMARY KEY,
    scope varchar(64)
);

DROP TABLE IF EXISTS annotation_value_type CASCADE;
CREATE TABLE annotation_value_type (
    annotation_value_type_id serial PRIMARY KEY,
    value_type varchar(64)
);

DROP TABLE IF EXISTS annotation_name CASCADE;
CREATE TABLE annotation_name (
    annotation_name_id serial PRIMARY KEY,
    name varchar,
    abbrev varchar(64),
    annotation_scope_id int,
    annotation_value_type_id int,
    FOREIGN KEY (annotation_scope_id) references annotation_scope,
    FOREIGN KEY (annotation_value_type_id) references annotation_value_type
);

DROP TABLE IF EXISTS boolean_gene_annotation;
CREATE TABLE boolean_gene_annotation (
    annotation_name_id int,
    gene_id int,
    value boolean,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (gene_id) references gene,
    PRIMARY KEY (annotation_name_id, gene_id)
);

DROP TABLE IF EXISTS int_gene_annotation;
CREATE TABLE int_gene_annotation (
    annotation_name_id int,
    gene_id int,
    value int,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (gene_id) references gene,
    PRIMARY KEY (annotation_name_id, gene_id)
);

DROP TABLE IF EXISTS double_gene_annotation;
CREATE TABLE double_gene_annotation (
    annotation_name_id int,
    gene_id int,
    value double precision,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (gene_id) references gene,
    PRIMARY KEY (annotation_name_id, gene_id)
);

DROP TABLE IF EXISTS text_gene_annotation;
CREATE TABLE text_gene_annotation (
    annotation_name_id int,
    gene_id int,
    value varchar(64),
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (gene_id) references gene,
    PRIMARY KEY (annotation_name_id, gene_id)
);

DROP TABLE IF EXISTS boolean_transcript_annotation;
CREATE TABLE boolean_transcript_annotation (
    annotation_name_id int,
    transcript_id int,
    value boolean,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (transcript_id) references transcript,
    PRIMARY KEY (annotation_name_id, transcript_id)
);

DROP TABLE IF EXISTS int_transcript_annotation;
CREATE TABLE int_transcript_annotation (
    annotation_name_id int,
    transcript_id int,
    value int,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (transcript_id) references transcript,
    PRIMARY KEY (annotation_name_id, transcript_id)
);

DROP TABLE IF EXISTS double_transcript_annotation;
CREATE TABLE double_transcript_annotation (
    annotation_name_id int,
    transcript_id int,
    value double precision,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (transcript_id) references transcript,
    PRIMARY KEY (annotation_name_id, transcript_id)
);

DROP TABLE IF EXISTS text_transcript_annotation;
CREATE TABLE text_transcript_annotation (
    annotation_name_id int,
    transcript_id int,
    value varchar(64),
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (transcript_id) references transcript,
    PRIMARY KEY (annotation_name_id, transcript_id)
);

DROP TABLE IF EXISTS boolean_cds_annotation;
CREATE TABLE boolean_cds_annotation (
    annotation_name_id int,
    cds_id int,
    value boolean,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (cds_id) references cds,
    PRIMARY KEY (annotation_name_id, cds_id)
);

DROP TABLE IF EXISTS int_cds_annotation;
CREATE TABLE int_cds_annotation (
    annotation_name_id int,
    cds_id int,
    value int,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (cds_id) references cds,
    PRIMARY KEY (annotation_name_id, cds_id)
);

DROP TABLE IF EXISTS double_cds_annotation;
CREATE TABLE double_cds_annotation (
    annotation_name_id int,
    cds_id int,
    value double precision,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (cds_id) references cds,
    PRIMARY KEY (annotation_name_id, cds_id)
);

DROP TABLE IF EXISTS text_cds_annotation;
CREATE TABLE text_cds_annotation (
    annotation_name_id int,
    cds_id int,
    value varchar(64),
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (cds_id) references cds,
    PRIMARY KEY (annotation_name_id, cds_id)
);

DROP TABLE IF EXISTS boolean_chr_range_annotation;
CREATE TABLE boolean_chr_range_annotation (
    annotation_name_id int,
    chromosome_id int,
    start_position int,
    end_position int,
    value boolean,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS int_chr_range_annotation;
CREATE TABLE int_chr_range_annotation (
    annotation_name_id int,
    chromosome_id int,
    start_position int,
    end_position int,
    value int,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS double_chr_range_annotation;
CREATE TABLE double_chr_range_annotation (
    annotation_name_id int,
    chromosome_id int,
    start_position int,
    end_position int,
    value double precision,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS text_chr_range_annotation;
CREATE TABLE text_chr_range_annotation (
    annotation_name_id int,
    chromosome_id int,
    start_position int,
    end_position int,
    value varchar(64),
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS boolean_site_annotation;
CREATE TABLE boolean_site_annotation (
    annotation_name_id int,
    chromosome_id int,
    position int,
    value boolean,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS int_site_annotation;
CREATE TABLE int_site_annotation (
    annotation_name_id int,
    chromosome_id int,
    position int,
    value int,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS double_site_annotation;
CREATE TABLE double_site_annotation (
    annotation_name_id int,
    chromosome_id int,
    position int,
    value double precision,
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS text_site_annotation;
CREATE TABLE text_site_annotation (
    annotation_name_id int,
    chromosome_id int,
    position int,
    value varchar(64),
    FOREIGN KEY (annotation_name_id) references annotation_name,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS allele_freq_scope CASCADE;
CREATE TABLE allele_freq_scope (
    allele_freq_scope_id serial PRIMARY KEY,
    scope varchar(64)
);

DROP TABLE IF EXISTS allele_freq;
CREATE table allele_freq (
    allele_freq_scope_id int,
    chromosome_id int,
    position int,
    freq real,
    derived_allele boolean,
    FOREIGN KEY (allele_freq_scope_id) references allele_freq_scope,
    FOREIGN KEY (chromosome_id) references chromosome
);

DROP TABLE IF EXISTS gwas_analysis CASCADE;
CREATE TABLE gwas_analysis (
    gwas_analysis_id serial PRIMARY KEY,
    trait varchar(64),
    germplasm varchar(64),
    markers varchar(64),
    description varchar
);

DROP TABLE IF EXISTS gwas_effect;
CREATE TABLE gwas_effect (
    gwas_analysis_id int,
    position int,
    effect real,
    in_model boolean,
    FOREIGN KEY (gwas_analysis_id) references gwas_analysis
);

DROP TABLE IF EXISTS qtl_analysis CASCADE;
CREATE TABLE qtl_analysis (
    qtl_analysis_id serial PRIMARY KEY,
    trait varchar(64),
    population varchar(64),
    markers varchar(64),
    description varchar
);

DROP TABLE IF EXISTS qtl_peak;
CREATE TABLE qtl_peak (
    qtl_peak_id serial PRIMARY KEY,
    qtl_analysis_id int,
    chromosome_id int,
    start_position int,
    peak_position int,
    end_position int,
    p_value real,
    FOREIGN KEY (qtl_analysis_id) references qtl_analysis,
    FOREIGN KEY (chromosome_id) references chromosome
);
