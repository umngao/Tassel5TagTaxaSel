-- Table: tag
CREATE TABLE tag (
    tagid    INTEGER PRIMARY KEY,
    tagName VARCHAR,
    sequence BLOB NOT NULL,
    seqlen INTEGER NOT NULL,
    qualityScore INTEGER,
    numTagInstances,
    UNIQUE (sequence, seqlen)
);

-- Table: tagPhysicalPosition
-- Junction (link) table between tag, and physicalposition
CREATE TABLE tagMapping (
    tagid       INTEGER NOT NULL,
    position_id  INTEGER NOT NULL,
    method_id    INTEGER NOT NULL,
    bp_error	INTEGER,
    cm_error	FLOAT(2),
    PRIMARY KEY (tagid, position_id)
);

-- Table: referenceGenome
-- Identifies reference genomes for physicalMapPosition table
CREATE TABLE referenceGenome (
	refid	INTEGER PRIMARY KEY,
	refname	TEXT
);

-- Table: physicalMapPosition
CREATE TABLE physicalMapPosition (
    posid INTEGER   PRIMARY KEY,
    reference_genome_id	INTEGER,
    chromosome TEXT      NOT NULL,
    physical_position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
);
CREATE UNIQUE INDEX physpos_idx ON physicalMapPosition(chromosome,physical_position,strand);
CREATE INDEX phychrpos_idx ON physicalMapPosition(chromosome);

-- Table: mappingApproach
CREATE TABLE mappingApproach (
    mapappid INTEGER   PRIMARY KEY,
    approach TEXT NOT NULL UNIQUE,
    software   TEXT NOT NULL,
    protocols   TEXT NOT NULL
);

-- Table: allelepair used for holding Alleles
-- Junction (link) table between tag and alleles
CREATE TABLE allelepair (
    allelepairid        INTEGER PRIMARY KEY,
    tagid_1     INTEGER NOT NULL,
    tagid_2	INTEGER,
    qualityscore INTEGER (1),
    dist_bp	INTEGER,
    dist_cm FLOAT(2)
);
CREATE INDEX newalleleidta_idx on allelepair(tagid_1);

-- Table: tagsnp used for holding Alleles
-- Junction (link) table between tag and snpposition
CREATE TABLE allele (
  alleleid        INTEGER   PRIMARY KEY,
  snpid     INTEGER NOT NULL,
  allelecall         INTEGER(1) NOT NULL,
  qualityscore INTEGER (1)
);
CREATE UNIQUE INDEX snpidallcall_idx on allele (snpid, allelecall);
CREATE INDEX newsnpidallcall_idx on allele (snpid);

-- Table: SNP Position
-- Is the refAllele "NOT NULL" required?
CREATE TABLE snpposition (
    snpid INTEGER   PRIMARY KEY,
    chromosome TEXT      NOT NULL,
    position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL,
    qualityscore FLOAT(1),
    refAllele         INTEGER(1) NOT NULL
);
CREATE UNIQUE INDEX snppos_idx ON snpposition(chromosome,position,strand);

-- Table: tagtaxadistribution
CREATE TABLE tagtaxadistribution (
    tagtxdstid  INTEGER   PRIMARY KEY,
    tagid      INTEGER NOT NULL,
    depthsRLE  BLOB,
    totalDepth INTEGER
);
CREATE INDEX tagid_idx ON tagtaxadistribution(tagid);


-- Table: taxa
CREATE TABLE taxa (
    taxonid INTEGER PRIMARY KEY,
    name    TEXT NOT NULL
);

-- Table: SNP Quality
--
CREATE TABLE snpQuality (
  snpid INTEGER   NOT NULL,
  taxasubset TEXT      NOT NULL,
  avgDepth REAL NOT NULL,
  minorDepthProp REAL NOT NULL,
  minor2DepthProp REAL NOT NULL,
  gapDepthProp REAL NOT NULL,
  propCovered REAL NOT NULL,
  propCovered2 REAL NOT NULL,
  taxaCntWithMinorAlleleGE2 REAL NOT NULL,
  minorAlleleFreqGE2      REAL NOT NULL,
  inbredF_DGE2    REAL,
  externalPositive REAL,
  externalNegative REAL
);
CREATE UNIQUE INDEX snpqual_idx ON snpQuality(snpid, taxasubset);
