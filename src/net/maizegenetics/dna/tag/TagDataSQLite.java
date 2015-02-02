package net.maizegenetics.dna.tag;

import com.google.common.collect.*;
import com.google.common.io.CharStreams;

import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import org.json.simple.JSONObject;
import org.sqlite.SQLiteConfig;

import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Defines xxxx
 * TAS-480
 * @author Ed Buckler
 */
public class TagDataSQLite implements TagDataWriter, AutoCloseable {
    private Connection connection = null;

    /*These maps contain  objects that are most queried by users.  This is a not the simplest way to do this, which
    would probably be done cleaner just with queries against the databases.  However, there are large performance
    gains in the case of SQLite (or at least with my ability to optimize).

    The logic behind this most of the datasets are relatively small (except tagTagIDMap), and this prevents creation
    of these objects over and over again.
     */
    private BiMap<Tag,Integer> tagTagIDMap;
    private Map<String,Integer> mappingApproachToIDMap;
    private SortedMap<Position,Integer> cutPosToIDMap;
    private BiMap<Position,Integer> snpPosToIDMap;
    private BiMap<Allele,Integer> alleleToIDMap;

    private TaxaList myTaxaList;

    PreparedStatement tagTaxaDistPS;
    PreparedStatement tagAlleleWhereTagPS;
    PreparedStatement tagidWhereSNPidPS;
    PreparedStatement tagidWhereAlleleidPS;
    PreparedStatement posTagInsertPS;
    PreparedStatement taxaDistWhereCutPositionIDPS;
    PreparedStatement snpPositionsForChromosomePS;
    PreparedStatement alleleTaxaDistForSnpidPS;
    PreparedStatement snpQualityInsertPS;




    public TagDataSQLite(String filename) {
        try{
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        // create a database connection

        try {
            boolean doesDBExist= Files.exists(Paths.get(filename));

            SQLiteConfig config=new SQLiteConfig();
            //config.setSynchronous(SQLiteConfig.SynchronousMode.OFF);
            /*Optimization ideas
            sqlite3_exec(mDb, "PRAGMA synchronous=OFF", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA count_changes=OFF", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA journal_mode=MEMORY", NULL, NULL, &errorMessage);
            sqlite3_exec(mDb, "PRAGMA temp_store=MEMORY", NULL, NULL, &errorMessage);
             */
            //config.setTempStore(SQLiteConfig.TempStore.MEMORY);
            connection = DriverManager.getConnection("jdbc:sqlite:"+filename,config.toProperties());
            connection.setAutoCommit(true);  //This has massive performance effects
            Statement statement = connection.createStatement();
            statement.setQueryTimeout(30);  // set timeout to 30 sec.
//            System.out.println(schema);
            if(doesDBExist==false) {
                String schema = CharStreams.toString(new InputStreamReader(TagDataSQLite.class.getResourceAsStream("TagSchema.sql")));
                statement.executeUpdate(schema);
            }
            initPreparedStatements();
            loadTagHash();
            loadMappingApproachHash();
            loadTaxaList();
        }
        catch(Exception e)
        {
            // if the error message is "out of memory",
            // it probably means no database file is found
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }

    @Override
    public void close() throws Exception {
        System.out.println("Closing SQLDB");
        connection.close();
    }

    private void initPreparedStatements() {
        try{
            posTagInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into tagCutPosition (tagid, positionid, mapappid, bestmapping, forward, cigar, supportval)" +
                            " values(?,?,?,?,?,?,?)");
            tagTaxaDistPS=connection.prepareStatement("select depthsRLE from tagtaxadistribution where tagid=?");
            tagAlleleWhereTagPS=connection.prepareStatement("select * from tagallele where tagid=?");
            tagidWhereSNPidPS=connection.prepareStatement(
                    "select tagid from allele, tagallele where allele.snpid=? and allele.alleleid=tagallele.alleleid");
            tagidWhereAlleleidPS=connection.prepareStatement(
                    "select tagid from tagallele where alleleid=?");
            taxaDistWhereCutPositionIDPS=connection.prepareStatement(
                    "select tagtaxadistribution.* from tagCutPosition, tagtaxadistribution where tagCutPosition.positionid=? and " +
                            "tagCutPosition.tagid=tagtaxadistribution.tagid and tagCutPosition.bestmapping=1");
            snpPositionsForChromosomePS=connection.prepareStatement(
            		"select position from snpposition where chromosome=?");
            alleleTaxaDistForSnpidPS =connection.prepareStatement("select a.*, td.* from allele a, tagallele ta, tagtaxadistribution td\n" +
                    "where a.alleleid=ta.alleleid and ta.tagid=td.tagid and a.snpid=?");
            snpQualityInsertPS=connection.prepareStatement(
                    "INSERT into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
                            "propCovered, propCovered2, taxaCntWithMinorAlleleGE2, minorAlleleFreq, inbredF_DGE2)" +
                            " values(?,?,?,?,?,?,?,?,?,?,?)");
//            snpQualityInsertPS=connection.prepareStatement(
//                    "INSERT into snpQuality (snpid, taxasubset)" +
//                            " values(?,?)");
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTagHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from tag");
            int size=rs.getInt(1);
            System.out.println("size of all tags in tag table="+size);
            if(tagTagIDMap==null || size/(tagTagIDMap.size()+1)>3) tagTagIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from tag");
            while(rs.next()) {
                tagTagIDMap.putIfAbsent(TagBuilder.instance(rs.getBytes("sequence"),rs.getShort("seqlen")).build(),rs.getInt("tagid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadCutPositionHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from cutPosition");
            int size=rs.getInt(1);
            System.out.println("size of all positions in cutPosition table="+size);
            if(cutPosToIDMap==null) {cutPosToIDMap=new TreeMap<>();}
            else if(size==cutPosToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from cutPosition");
            while(rs.next()) {
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("position"))
                        .strand(rs.getByte("strand"))
                        .build();
                cutPosToIDMap.putIfAbsent(p, rs.getInt("positionid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadSNPPositionHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from snpposition");
            int size=rs.getInt(1);
            System.out.println("size of all positions in snpPosition table="+size);
            if(snpPosToIDMap==null) {snpPosToIDMap=HashBiMap.create(size);}
            else if(size==snpPosToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from snpposition");
            while(rs.next()) {
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("position"))
                        .strand(rs.getByte("strand"))
                        .build();
                snpPosToIDMap.putIfAbsent(p, rs.getInt("snpid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadAlleleHash() {
        try{
            loadSNPPositionHash();
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from allele");
            int size=rs.getInt(1);
            System.out.println("size of all alleles in allele table="+size);
            if(alleleToIDMap==null) {alleleToIDMap=HashBiMap.create(size);}
            if(size==alleleToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from allele");
            while(rs.next()) {
                int snpid=rs.getInt("snpid");
                Position p=snpPosToIDMap.inverse().get(snpid);
                Allele a=new SimpleAllele(rs.getByte("allelecall"),p);
                alleleToIDMap.putIfAbsent(a, rs.getInt("alleleid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadMappingApproachHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from mappingApproach");
            int size=rs.getInt(1);
            System.out.println("size of all tags in mappingApproach table="+size);
            if(size==0) {
                connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                        "values('unknown','unknown','unknown')");
                size=1;
            }
            mappingApproachToIDMap=new HashMap<>(size);
            rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                mappingApproachToIDMap.put(rs.getString("approach"),rs.getInt("mapappid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTaxaList() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from taxa");
            int size=rs.getInt(1);
            System.out.println("size of all taxa in taxa table="+size);
            TaxaListBuilder tlb=new TaxaListBuilder();
            rs=connection.createStatement().executeQuery("select * from taxa");
            while(rs.next()) {
                tlb.add(new Taxon(rs.getString("name")));
            }
            myTaxaList=tlb.build();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public TaxaList getTaxaList() {
        if(myTaxaList==null) loadTaxaList();
        return myTaxaList;
    }

    @Override
    public boolean putAllTag(Set<Tag> tags) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tag (sequence, seqlen) values(?,?)");
            for (Tag tag : tags) {
                if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                tagInsertPS.setShort(2, tag.seqLength());
                tagInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                    System.out.println("tagInsertPS.executeBatch() "+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        if(totalCount>0) loadTagHash();
        return true;
    }

    @Override
    public void putTaxaList(TaxaList taxaList) {
        try {
            connection.createStatement().execute("delete from taxa");
            connection.setAutoCommit(false);
            PreparedStatement taxaInsertPS=connection.prepareStatement("insert into taxa (taxonid, name) values(?,?)");
            for (int i = 0; i < taxaList.size(); i++) {
                taxaInsertPS.setInt(1, i);
                taxaInsertPS.setString(2, taxaList.get(i).getName());
                taxaInsertPS.addBatch();

            }
            taxaInsertPS.executeBatch();
            connection.setAutoCommit(true);
            loadTaxaList();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap) {
        int batchCount=0;
        try {
            int numTaxa=myTaxaList.numberOfTaxa();
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tagtaxadistribution (tagid, depthsRLE, totalDepth) values(?,?,?)");
            for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaDistributionMap.entrySet()) {
                int tagID=tagTagIDMap.get(entry.getKey());
                tagInsertPS.setInt(1,tagID);
                if(entry.getValue().maxTaxa()!=numTaxa) throw new IllegalStateException("Number of taxa does not agree with taxa distribution");
                tagInsertPS.setBytes(2, entry.getValue().encodeTaxaDepth());
                tagInsertPS.setInt(3, entry.getValue().totalDepth());
                tagInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                    System.out.println("putTaxaDistribution next"+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    @Override
    public void putTagAlignments(Multimap<Tag, Position> tagAnnotatedPositionMap) {
        int batchCount=0;
        try {
            putAllTag(tagAnnotatedPositionMap.keySet());
            putCutPositionsIfAbsent(tagAnnotatedPositionMap.values());
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Position> entry : tagAnnotatedPositionMap.entries()) {
                Position p=entry.getValue();
                int ind=1;
                posTagInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                posTagInsertPS.setInt(ind++, cutPosToIDMap.get(p));
                posTagInsertPS.setInt(ind++, getMappingApproachID(p));
                posTagInsertPS.setBoolean(ind++, true);  //todo this needs to be derived from the position or set later.
                boolean forward=true;
                try{
                    if(p.getTextAnnotation("forward")[0].toLowerCase().equals("false")) forward=false;
                } catch (Exception e) {
                    System.err.println(p.toString());
                    System.err.println("Error with forward annotation");
                    //no valid cigarValue
                }
                posTagInsertPS.setBoolean(ind++,forward);
                String cigarValue="";
                try{
                    cigarValue=p.getTextAnnotation("cigar")[0];
                } catch (Exception e) {
                    System.err.println(p.toString());
                    System.err.println("Error with cigar");
                    //no valid cigarValue
                }
                posTagInsertPS.setString(ind++, cigarValue);
                short supportVal=0;
                try{
                    String[] svS=p.getTextAnnotation("supportvalue");
                    if(svS.length>0) {
                        supportVal=Short.parseShort(svS[0]);
                    }
                } catch (Exception e) {
                    System.err.println("Error with supportVal");
                    //no valid supportVal
                }
                posTagInsertPS.setByte(ind++, (byte) supportVal);
                //System.out.println(posTagInsertPS.toString());
                posTagInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putTagAlignments next"+batchCount);
                    posTagInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posTagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    //@Override
    public void putSNPQualityProfile(Map<Position, Map<String,Double>> tagAnnotatedPositionMap, String taxaSubset) {
        int batchCount=0;
        try {
            putSNPPositionsIfAbsent(tagAnnotatedPositionMap.keySet());
            connection.setAutoCommit(false);
            for (Map.Entry<Position, Map<String,Double>> entry : tagAnnotatedPositionMap.entrySet()) {
                Position p=entry.getKey();
                Map<String,Double> vals=entry.getValue();
                int ind=1;

//                snpQualityInsertPS=connection.prepareStatement(
//                        "INSERT OR IGNORE into snpQuality (snpid, taxasubset ,avgDepth, minorDepthProp, minor2DepthProp, gapDepthProp, " +
//                                "propCovered, propCovered2, taxaCntWithMinorAlleleGE2,minorAlleleFreq, inbredF_DGE2)" +
//                                " values(?,?,?,?,?,?,?,?,?,?,?)");
                System.out.println("vals = " + vals.size());
                System.out.println(vals.get("inbredF_DGE2").toString());
                snpQualityInsertPS.setInt(ind++, snpPosToIDMap.get(p));
                snpQualityInsertPS.setString(ind++, taxaSubset);
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("avgDepth",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minor2DepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("gapDepthProp",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("propCovered2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("taxaCntWithMinorAlleleGE2",0.0));
                System.out.println("MAF:"+vals.getOrDefault("minorAlleleFreqGE2",-1.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("minorAlleleFreqGE2",0.0));
                snpQualityInsertPS.setDouble(ind++,vals.getOrDefault("inbredF_DGE2",null));
                //snpQualityInsertPS.get();
                //System.out.println(posTagInsertPS.toString());
                snpQualityInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putSNPQualityProfile next"+batchCount);
                    snpQualityInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            snpQualityInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    private int getMappingApproachID(Position p) throws SQLException{
        String mapApp=p.getTextAnnotation("mappingapproach")[0];
        if(mapApp==null) return mappingApproachToIDMap.get("unknown");
        Integer val=mappingApproachToIDMap.get(mapApp);
        if(val==null) {
            connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                    "values('"+mapApp+"','unknown','unknown')");
            loadMappingApproachHash();
            return mappingApproachToIDMap.get(mapApp);
        } else return val;
    }

    @Override
    public void setTagAlignmentBest(Tag tag, Position position, boolean isBest) {

    }

    @Override
    public boolean putTagAlleles(Multimap<Tag, Allele> tagAlleleMap) {
        int batchCount=0;
        try {
            PreparedStatement alleleTagInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into tagallele (alleleid, tagid) values(?,?)");
            putAllTag(tagAlleleMap.keySet());
            loadSNPPositionHash();
            putSNPPositionsIfAbsent(tagAlleleMap.values().stream()
                    .map(a -> a.position())
                    .distinct()
                    .collect(Collectors.toSet()));
            putAlleleIfAbsent(tagAlleleMap.values().stream()
                    .distinct()
                    .collect(Collectors.toSet()));
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Allele> tagAlleleEntry : tagAlleleMap.entries()) {
                int ind=1;
                alleleTagInsertPS.setInt(ind++, alleleToIDMap.get(tagAlleleEntry.getValue()));
                alleleTagInsertPS.setInt(ind++, tagTagIDMap.get(tagAlleleEntry.getKey()));
                alleleTagInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("alleleTagInsertPS next"+batchCount);
                    alleleTagInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            alleleTagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    @Override
    public boolean putTagAlignmentApproach(String tagAlignmentName, String protocol) {
        return false;
    }

    @Override
    public TaxaDistribution getTaxaDistribution(Tag tag) {
        int tagid=tagTagIDMap.get(tag);
        try {
            tagTaxaDistPS.setInt(1,tagid);
            ResultSet rs=tagTaxaDistPS.executeQuery();
            return TaxaDistBuilder.create(rs.getBytes(1));
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public Set<Allele> getAlleles(Tag tag) {if(alleleToIDMap==null) loadAlleleHash();
        ImmutableSet.Builder<Allele> alleleBuilder=new ImmutableSet.Builder<>();
        try{
            tagAlleleWhereTagPS.setInt(1,tagTagIDMap.get(tag));
            ResultSet rs=tagAlleleWhereTagPS.executeQuery();
            while(rs.next()) {
                alleleBuilder.add(alleleToIDMap.inverse().get(rs.getInt("alleleid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return alleleBuilder.build();
    }

    @Override
    public Multimap<Tag, Allele> getAlleleMap() {
        //if slow consider caching the hash and equal codes
        if(alleleToIDMap==null) loadAlleleHash();
        ImmutableMultimap.Builder<Tag, Allele> tagAlleleBuilder=new ImmutableMultimap.Builder<>();
        try{
            ResultSet rs=connection.createStatement().executeQuery("select * from tagallele");;
            while(rs.next()) {
                tagAlleleBuilder.put(tagTagIDMap.inverse().get(rs.getInt("tagid")), alleleToIDMap.inverse().get(rs.getInt("alleleid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagAlleleBuilder.build();
    }

    @Override
    public Set<Tag> getTagsForSNPPosition(Position position) {
        ImmutableSet.Builder<Tag> tagBuilder=new ImmutableSet.Builder<>();
        try{
            tagidWhereSNPidPS.setInt(1,snpPosToIDMap.get(position));
            ResultSet rs=tagAlleleWhereTagPS.executeQuery();
            while(rs.next()) {
                tagBuilder.add(tagTagIDMap.inverse().get(rs.getInt("tagid")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public Set<Tag> getTagsForAllele(Position position, byte allele) {
        return getTagsForAllele(new SimpleAllele(allele,position));
    }

    @Override
    public Set<Tag> getTagsForAllele(Allele allele) {
        return null;
    }

    public Multimap<Allele,TaxaDistribution> getAllelesTaxaDistForSNP(Position position) {
        ImmutableMultimap.Builder<Allele,TaxaDistribution> atdBuilder=ImmutableMultimap.builder();
        try{
            alleleTaxaDistForSnpidPS.setInt(1, snpPosToIDMap.get(position));
            ResultSet rs= alleleTaxaDistForSnpidPS.executeQuery();
            while(rs.next()) {
                Allele allele=new SimpleAllele((byte)rs.getInt("allelecall"),position);
                atdBuilder.put(allele,TaxaDistBuilder.create(rs.getBytes("depthsRLE")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return atdBuilder.build();
    }

    @Override
    public Set<Tag> getTags() {
        return tagTagIDMap.keySet();
    }

    @Override
    public PositionList getSNPPositions() {
        if(snpPosToIDMap==null) loadSNPPositionHash();
        return new PositionListBuilder().addAll(snpPosToIDMap.keySet()).build();
    }

    @Override
    public PositionList getSNPPositions(int minSupportValue) {
        return null;
    }

    @Override
    public Set<Tag> getTagsForTaxon(Taxon taxon) {
        ImmutableSet.Builder<Tag> tagBuilder=new ImmutableSet.Builder<>();
        int taxonIndex=myTaxaList.indexOf(taxon);
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from tagtaxadistribution");
            while(rs.next()) {
                if(TaxaDistBuilder.create(rs.getBytes("depthsRLE")).depths()[taxonIndex]>0) {
                    tagBuilder.add(tagTagIDMap.inverse().get(rs.getInt("tagid")));
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public Map<Tag, Integer> getTagDepth(Taxon taxon, Position position) {
        return null;
    }

    @Override
    public Map<Tag, Integer> getTagsWithDepth(int minimumDepth) {
        ImmutableMap.Builder<Tag, Integer> tagBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery(
                    "select tagid, totalDepth from tagtaxadistribution where totalDepth >= "+minimumDepth);
            while(rs.next()) {
                tagBuilder.put(tagTagIDMap.inverse().get(rs.getInt(1)),rs.getInt(2));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public PositionList getTagCutPositions(boolean onlyBest) {
        if(cutPosToIDMap.size()==0) loadCutPositionHash();
        PositionListBuilder plb=new PositionListBuilder();
        cutPosToIDMap.keySet().stream()
                //.filter(p -> p.isAnnotatedWithValue("isbest","true"))
                .forEach(p -> plb.add(p));  //todo only best not implemented here
        plb.sortPositions();
        return plb.build();
    }

    @Override
    public PositionList getTagCutPositions(Chromosome chromosome, int firstPosition, int lastPosition, boolean onlyBest) {
        PositionListBuilder plb=new PositionListBuilder();
        plb.addAll(getPositionSubMap(chromosome,firstPosition,lastPosition).keySet());  //todo only best not implemented here
        return plb.build();
    }

    private Map<Position,Integer> getPositionSubMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        if(cutPosToIDMap==null) loadCutPositionHash();
        Position startPos=new GeneralPosition.Builder(chromosome,firstPosition).build();
        if(lastPosition<0) lastPosition=Integer.MAX_VALUE;
        Position lastPos=new GeneralPosition.Builder(chromosome,lastPosition).build();
        return cutPosToIDMap.subMap(startPos,lastPos);
    }

    @Override
    public Map<String, String> getTagAlignmentApproaches() {
        ImmutableMap.Builder<String,String> appBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                appBuilder.put(rs.getString("approach"), rs.getString("software") + ":" + rs.getString("approach"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return appBuilder.build();
    }

//    @Override
//    public Map<Position, Map<Tag, TaxaDistribution>> getCutPositionTagTaxaMapX(Chromosome chromosome, int firstPosition, int lastPosition) {
//        //consider doing this all by SQL if performance suffers
//        PositionList pl=getTagCutPositions(chromosome,firstPosition,lastPosition,true);
//        ImmutableMap.Builder<Position, Map<Tag, TaxaDistribution>> positionMapBuilder=new ImmutableMap.Builder<>();
//        pl.stream().forEach(p -> positionMapBuilder.put(p,getTagsTaxaMap(p)));
//        //this is slow as each position is a separate transaction
//        return positionMapBuilder.build();
//    }

    //TODO need to add the forward direction to the resulting map somehow.  Perhaps  Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>>
    //alternatively there could be a tag alignment object.

    @Override
    public Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        String sqlQuery="select p.positionid, forward, chromosome, position, strand, t.tagid, depthsRLE  " +
                "from tag t, cutposition p, tagCutPosition tc, tagtaxadistribution ttd " +
                "where p.positionid=tc.positionid and tc.tagid=t.tagid and t.tagid=ttd.tagid " +
                "and chromosome="+chromosome.toString()+//" and position>"+firstPosition+" " + //todo position would need to be index to make fast
                " order by position";
        Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> positionTagTaxaMap=new HashMap<>();
        Map<Integer,Position> tempPositionMap=new HashMap<>();  //reverse the map
        getPositionSubMap(chromosome,firstPosition,lastPosition).entrySet().stream()
                .forEach(entry -> tempPositionMap.put(entry.getValue(),entry.getKey()));
        try{
            ResultSet rs=connection.createStatement().executeQuery(sqlQuery);
            while(rs.next()) {
                Position position=tempPositionMap.get(rs.getInt("positionid"));
                Tag tag=tagTagIDMap.inverse().get(rs.getInt("tagid"));
                TaxaDistribution taxaDistribution=TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                Boolean forwardAlignDirection=rs.getBoolean("forward");
                Map<Tag, Tuple<Boolean,TaxaDistribution>> tagTaxaMap=positionTagTaxaMap.computeIfAbsent(position, k -> new HashMap<>());
                tagTaxaMap.put(tag,new Tuple<>(forwardAlignDirection,taxaDistribution));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        System.out.println("positionTagTaxaMap = " + positionTagTaxaMap.size());
        return positionTagTaxaMap;
    }


    @Override
    public Map<Tag, TaxaDistribution> getTagsTaxaMap(Position cutPosition) {
        ImmutableMap.Builder<Tag, TaxaDistribution> tagTaxaDistributionBuilder=new ImmutableMap.Builder<>();
        try{
            taxaDistWhereCutPositionIDPS.setInt(1,cutPosToIDMap.get(cutPosition));
            ResultSet rs=taxaDistWhereCutPositionIDPS.executeQuery();
            while(rs.next()) {
                tagTaxaDistributionBuilder.put(tagTagIDMap.inverse().get(rs.getInt("tagid")), TaxaDistBuilder.create(rs.getBytes("depthsRLE")));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagTaxaDistributionBuilder.build();
    }

    private void putCutPositionsIfAbsent(Collection<Position> positions) {
        try {
        int batchCount=0;
        if(cutPosToIDMap==null) loadCutPositionHash();
        connection.setAutoCommit(false);
        PreparedStatement posInsertPS=connection.prepareStatement(
                "INSERT OR IGNORE into cutposition (chromosome, position, strand) values(?,?,?)");
        for (Position p : positions) {
            if(cutPosToIDMap.containsKey(p)) continue;
            posInsertPS.setString(1, p.getChromosome().toString());
            posInsertPS.setInt(2, p.getPosition());
            posInsertPS.setByte(3, p.getStrand());
            posInsertPS.addBatch();
            batchCount++;
            if(batchCount>10000) {
                System.out.println("putCutPositionsIfAbsent next"+batchCount);
                posInsertPS.executeBatch();
                batchCount=0;
            }
        }
        posInsertPS.executeBatch();
        if(batchCount>0) loadCutPositionHash();
        connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    private void putSNPPositionsIfAbsent(Collection<Position> positions) {
        try {
            int batchCount=0;
            if(snpPosToIDMap==null) loadSNPPositionHash();
            connection.setAutoCommit(false);
            PreparedStatement snpPosInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into snpposition (chromosome, position, strand) values(?,?,?)");
            for (Position p : positions) {
                if(snpPosToIDMap.containsKey(p)) continue;
                snpPosInsertPS.setString(1, p.getChromosome().toString());
                snpPosInsertPS.setInt(2, p.getPosition());
                snpPosInsertPS.setByte(3, p.getStrand());
                snpPosInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putSNPPositionsIfAbsent next"+batchCount);
                    snpPosInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            snpPosInsertPS.executeBatch();
            if(batchCount>0) loadSNPPositionHash();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    private void putAlleleIfAbsent(Collection<Allele> alleles) throws IllegalStateException {
        try {
            int batchCount=0;
            PreparedStatement alleleInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into allele (snpid, allelecall, qualityscore) values(?,?,?)");
            connection.setAutoCommit(false);
            for (Allele allele : alleles) {
                Integer snpID=snpPosToIDMap.get(allele.position());
                if(snpID==null) throw new IllegalStateException("SNP position missing for allele");
                int index=1;
                alleleInsertPS.setInt(index++, snpID);
                alleleInsertPS.setByte(index++, allele.allele());
                alleleInsertPS.setByte(index++, (byte) 0);  //todo set quality scores once annotation pattern is set
                alleleInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putAlleleIfAbsent next"+batchCount);
                    alleleInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            alleleInsertPS.executeBatch();
            if(batchCount>0) loadAlleleHash();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
    
    // Return SNP positions for specified list of chromosomes
    @Override
    public  PositionList getSNPPositionsForChromosomes(Integer startChr, Integer endChr) {
    	PositionListBuilder plb = new PositionListBuilder();
    	// Verify good chromsome values
    	if (startChr < 1 ||
    	    endChr < 1 ||
    	    startChr > endChr ) {
    		System.err.printf("getSNPPOsitionsForChromosomes:  bad Chromosome values: startChr %d, endChr %d\n",
    				startChr, endChr);
    		return null;
    	}
        try{
        	for (int chrom = startChr; chrom <= endChr; chrom++ ){
        		snpPositionsForChromosomePS.setString(1, Integer.toString(chrom));
                ResultSet rs=snpPositionsForChromosomePS.executeQuery();
                while(rs.next()) {
                	Chromosome chr= new Chromosome(Integer.toString(chrom));
                	Position position = new GeneralPosition.Builder(chr,rs.getInt("position")).build();
                    plb.add(position);
                }
        	} 
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return plb.build();
    }
}
