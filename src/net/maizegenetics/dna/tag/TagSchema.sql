-- Table: tag
CREATE TABLE tag (
    tagid    INTEGER PRIMARY KEY,
    sequence BLOB NOT NULL,
    seqlen INTEGER NOT NULL,
    UNIQUE (sequence, seqlen)
);

-- Table: tagPosition
-- Junction (link) table between tag, position, and mapping approach
CREATE TABLE tagCutPosition (
    tagid       INTEGER NOT NULL,
    positionid  INTEGER NOT NULL,
    mapappid    INTEGER NOT NULL,
    bestmapping BOOLEAN,
    cigar       TEXT,
    supportval  INTEGER(1),
    PRIMARY KEY (tagid, positionid, mapappid)
);


-- Table: cutposition
CREATE TABLE cutposition (
    positionid INTEGER   PRIMARY KEY,
    chromosome TEXT      NOT NULL,
    position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
);
CREATE UNIQUE INDEX cutpos_idx ON cutposition(chromosome,position,strand);

-- Table: mappingApproach
CREATE TABLE mappingApproach (
    mapappid INTEGER   PRIMARY KEY,
    approach TEXT NOT NULL UNIQUE,
    software   TEXT NOT NULL,
    protocols   TEXT NOT NULL
);

-- Table: tagsnp used for holding Alleles
-- Junction (link) table between tag and snpposition
CREATE TABLE tagallele (
    tagid        INTEGER NOT NULL,
    alleleid     INTEGER NOT NULL,
    qualityscore INTEGER (1),
    PRIMARY KEY (tagid, alleleid)
);

-- Table: tagsnp used for holding Alleles
-- Junction (link) table between tag and snpposition
CREATE TABLE allele (
  alleleid        INTEGER   PRIMARY KEY,
  snpid     INTEGER NOT NULL,
  allelecall         INTEGER(1) NOT NULL,
  qualityscore INTEGER (1)
);


-- Table: SNP Position
CREATE TABLE snpposition (
    snpid INTEGER   PRIMARY KEY,
    chromosome TEXT      NOT NULL,
    position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
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

