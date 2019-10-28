package edu.caltech.lncrna.barcode.core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.barcode.core.Read;
import edu.caltech.lncrna.barcode.core.TagCategory;
import edu.caltech.lncrna.bio.io.FastqParser;
import edu.caltech.lncrna.bio.sequence.FastqSequence;
import edu.caltech.lncrna.bio.sequence.PhredEncoding;

public class BarcodeIdentification {

    private static final String VERSION = "1.2.0";
    private static final String NOT_FOUND = "[NOT_FOUND]";
    
    private final Path inputPath1;
    private final Path inputPath2;
    private final Path outputPath1;
    private final Path outputPath2;
    
    private final List<TagCategory> layout1;
    private final List<TagCategory> layout2;
    
    private final Map<TagCategory, Map<String, Tag>> tags;
    private final Map<TagCategory, int[]> tagLengths;
    
    private final int spacerLength;
    private final int maxAdvance;
    
    private static final Logger LOGGER =
            Logger.getLogger("BarcodeIdentification");
    
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        CommandLine cmd = parseArgs(args);

        new BarcodeIdentification(cmd).findBarcodes();
        LOGGER.info("Program complete.");
        LOGGER.info((System.currentTimeMillis() - startTime) + " milliseconds elapsed.");
    }
    
    public BarcodeIdentification(CommandLine cmd) {
        inputPath1 = Paths.get(cmd.getOptionValue("input1"));
        inputPath2 = Paths.get(cmd.getOptionValue("input2"));
        outputPath1 = Paths.get(cmd.getOptionValue("output1"));
        outputPath2 = Paths.get(cmd.getOptionValue("output2"));
       
        LOGGER.info("Parsing configuration file.");
        Configuration config = new Configuration(Paths.get(cmd.getOptionValue("config")));
       
        layout1 = Collections.unmodifiableList(config.layout1);
        layout2 = Collections.unmodifiableList(config.layout2);
       
        tags = Collections.unmodifiableMap(config.tags);
        tagLengths = Collections.unmodifiableMap(config.lengths);
        
        spacerLength = config.spacerLength;
        maxAdvance = config.laxity;
    }
    
    private static CommandLine parseArgs(String[] args) {
        
        Option versionOption = Option.builder("v")
                .longOpt("version")
                .desc("show version information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option helpOption = Option.builder("h")
                .longOpt("help")
                .desc("show usage information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option configOption = Option.builder()
                .longOpt("config")
                .desc("configuration file")
                .hasArg()
                .required()
                .build();
        
        Option input1Option = Option.builder()
                .longOpt("input1")
                .desc("input FASTQ file 1")
                .hasArg()
                .required()
                .build();
        
        Option input2Option = Option.builder()
                .longOpt("input2")
                .desc("input FASTQ file 2")
                .hasArg()
                .required()
                .build();
        
        Option output1Option = Option.builder()
                .longOpt("output1")
                .desc("output FASTQ file 1")
                .hasArg()
                .required()
                .build();

        Option output2Option = Option.builder()
                .longOpt("output2")
                .desc("output FASTQ file 2")
                .hasArg()
                .required()
                .build();
        
        Options helpOptions = new Options().addOption(helpOption);
        Options versionOptions = new Options().addOption(versionOption);
        
        Options mainOptions = new Options()
                .addOption(configOption)
                .addOption(input1Option)
                .addOption(input2Option)
                .addOption(output1Option)
                .addOption(output2Option);
        
        Options allOptions = new Options();
        helpOptions.getOptions().forEach(allOptions::addOption);
        versionOptions.getOptions().forEach(allOptions::addOption);
        mainOptions.getOptions().forEach(allOptions::addOption);
        HelpFormatter formatter = new HelpFormatter();
                
        CommandLineParser parser = new DefaultParser();
        CommandLine rtrn = null;
        
        try {
            CommandLine cmds = parser.parse(helpOptions, args, true);
            if (cmds.getOptions().length == 1) {
                formatter.printHelp("test", allOptions);
                System.exit(0);
            }
            cmds = parser.parse(versionOptions, args, true);
            if (cmds.getOptions().length == 1) {
                System.out.println(VERSION);
                System.exit(0);
            }
            rtrn = parser.parse(mainOptions, args);
        } catch (ParseException e) {
            formatter.printHelp("test", allOptions);
            System.exit(1);
        }
        
        if (rtrn == null) {
            LOGGER.severe("An unknown error occurred while parsing command line arguments");
            System.exit(1);
        }
        
        return rtrn;
    }
    
    private BarcodeIdentification findBarcodes() {

        try (FastqParser fq1 = new FastqParser(inputPath1, PhredEncoding.SANGER);
             FastqParser fq2 = new FastqParser(inputPath2, PhredEncoding.SANGER);
             BufferedWriter writer1 = new BufferedWriter(new OutputStreamWriter(
                     new GZIPOutputStream(new FileOutputStream(
                     outputPath1.toFile()))));
             BufferedWriter writer2 = new BufferedWriter(new OutputStreamWriter(
                     new GZIPOutputStream(new FileOutputStream(
                     outputPath2.toFile()))));
             ) {
            
             //BufferedWriter writer1 = Files.newBufferedWriter(outputPath1,
             //        StandardCharsets.US_ASCII, StandardOpenOption.CREATE);
             //BufferedWriter writer2 = Files.newBufferedWriter(outputPath2,
             //        StandardCharsets.US_ASCII, StandardOpenOption.CREATE)) {

            int fqCount = 0;
            while (fq1.hasNext() && fq2.hasNext()) {

                if ((++fqCount) % 500000 == 0) {
                    LOGGER.info("Processing read " + fqCount + ".");
                }

                FastqSequence read1 = fq1.next();
                FastqSequence read2 = fq2.next();
                
                StringBuilder sb = new StringBuilder();
                checkRead1(read1, sb);
                checkRead2(read2, sb);

                FastqSequence barcodedRead1 = appendBarcodesToName(read1, sb);
                FastqSequence barcodedRead2 = appendBarcodesToName(read2, sb);
                
                writer1.write(barcodedRead1.toFormattedString(PhredEncoding.SANGER));
                writer2.write(barcodedRead2.toFormattedString(PhredEncoding.SANGER));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        
        return this;
    }
    
    private void checkRead1(FastqSequence seq, StringBuilder sb) {
        checkRead(seq, Read.READ1, sb);
    }
    
    private void checkRead2(FastqSequence seq, StringBuilder sb) {
        checkRead(seq, Read.READ2, sb);
    }
    
    private void checkRead(FastqSequence seq, Read read, StringBuilder sb) {
        String bases = seq.getBases();
        List<TagCategory> layout = read.equals(Read.READ1) ? layout1 : layout2;
        boolean startOfRead = true;
        for (TagCategory category : layout) {
            if (category.equals(TagCategory.SPACER)) {
                bases = skipSpacer(bases);
            } else {
                bases = startOfRead
                      ? checkReadForTagAndReturnRemainder(bases, category, 0, sb)
                      : checkReadForTagAndReturnRemainder(bases, category, maxAdvance, sb);
            }
            startOfRead = false;
        }
    }

    private String checkReadForTagAndReturnRemainder(String bases, TagCategory category,
            int maxAdvance, StringBuilder sb) {
        
        Tag match = null;
        int[] lengths = tagLengths.get(category);
        assert lengths != null;
        int advance = -1;
        
        while (match == null && advance < maxAdvance) {
            advance++;
            for (int i = 0; i < lengths.length; i++) {
                int len = lengths[i];
                
                if (advance + len <= bases.length()) {
                    String possibleTag = bases.substring(advance, advance + len);
                    match = tags.get(category).get(possibleTag);
                    if (match != null) {
                        break;
                    }
                }
            }
        }
        
        if (match == null) {
            sb.append(NOT_FOUND);
            return advance + lengths[lengths.length - 1] > bases.length() ?
                "" : bases.substring(advance + lengths[lengths.length - 1]);
        } else {
            sb.append("[" + match.name() + "]");
            return bases.substring(advance + match.seq().length());
        }
    }
    
    private String skipSpacer(String bases) {
        return bases.length() < spacerLength ? "" : bases.substring(spacerLength);
    }
    
    private FastqSequence appendBarcodesToName(FastqSequence fq, StringBuilder sb) {
        return fq.changeName(fq.getName().split(" ")[0] + "::" + sb.toString());   
    }
    
    /**
     * A builder-like class for initializing a BarcodeIdentification instance
     * from a configuration file.
     */
    
    public static class Configuration {
        
        private final List<TagCategory> layout1 = new ArrayList<>();
        private final List<TagCategory> layout2 = new ArrayList<>();
        
        private final Map<TagCategory, Map<String, Tag>> tags = new HashMap<>();
        private final Map<TagCategory, int[]> lengths = new HashMap<>();
        
        private int spacerLength = UNINITIALIZED;
        private int laxity = UNINITIALIZED;
        
        private static final int NUM_SPACER_LENGTH_FIELDS = 3;
        private static final int NUM_LAXITY_FIELDS = 3;
        private static final int NUM_TAG_FIELDS = 4;
        private static final int UNINITIALIZED = -1;
        
        private static final int DEFAULT_SPACER_LENGTH = 6;
        private static final int DEFAULT_LAXITY = 6;
        
        private static final Logger LOGGER = Logger.getLogger("Configuration");
        
        public Configuration(Path configPath) {
            parseConfigFile(configPath);
        }
        
        private void parseConfigFile(Path configFile) {
            try (BufferedReader br = Files.newBufferedReader(configFile, StandardCharsets.US_ASCII)) {
                String line = null;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#") || line.isEmpty()) {
                        continue;
                    }

                    if (line.startsWith("READ1 = ")) {
                        parseTagLayout(layout1, line);
                    } else if (line.startsWith("READ2 = ")) {
                        parseTagLayout(layout2, line);
                    } else if (line.startsWith("SPACER = ")) {
                        parseSpacer(line);
                    } else if (line.startsWith("LAXITY = ")) {
                        parseLaxity(line);
                    } else {
                        parseTag(line);
                    }
                }
            } catch (IOException ex) {
                LOGGER.log(Level.SEVERE, "An error occurred while reading the "
                        + "configuration file.", ex);
            }
            
            populateLengths();
            
            if (laxity == UNINITIALIZED) {
                laxity = DEFAULT_LAXITY;
            }
            
            if (spacerLength == UNINITIALIZED) {
                spacerLength = DEFAULT_SPACER_LENGTH;
            }
            
            if (layout1.isEmpty() && layout2.isEmpty()) {
                LOGGER.severe("No tag layouts found in configuration file. Exiting.");
                System.exit(0);
            }
            
            if (tags.isEmpty()) {
                LOGGER.severe("No tags found in configuration file. Exiting.");
                System.exit(0);
            }
        }
        
        //////////////////////////////////////
        // Code for parsing individual tags //
        //////////////////////////////////////
        
        private void parseTag(String line) {
            String [] fields = line.split("\t");
            
            if (fields.length != NUM_TAG_FIELDS) {
                logConfigErrorAndExit(line);
            }
            
            TagCategory category = null;
            try {
                category = TagCategory.valueOf(fields[0].toUpperCase());
            } catch (IllegalArgumentException e) {
                LOGGER.log(Level.SEVERE, "Encountered invalid tag category " +
                        "in configuration file. See line: " + line, e);
                System.exit(0);
            }
            
            String name = fields[1];
            String seq = fields[2];
            
            int mismatches = 0;
            try {
                mismatches = Integer.parseInt(fields[3]);
            } catch (NumberFormatException e) {
                logConfigErrorAndExit(e, line);
            }
            
            if (mismatches < 0) {
                logConfigErrorAndExit(line);
            }
            
            Tag tag = new Tag(category, name, seq);
            
            addTag(tag, mismatches);
        }
        
        private void addTag(Tag tag, int mismatches) {
            
            assert mismatches >= 0;
            
            // Add empty map for tag category if it doesn't exist yet.
            TagCategory category = tag.category();
            if (!tags.containsKey(category)) {
                tags.put(category, new HashMap<String, Tag>());
            }
            
            if (mismatches == 0) {
                addExactMatch(tag);
            } else {
                addWithMismatches(tag, mismatches);
            }
        }
        
        private void addExactMatch(Tag tag) {
            TagCategory category = tag.category();
            Map<String, Tag> tagsByCat = tags.get(category);
            warnIfAmbiguity(tagsByCat, tag);
            tagsByCat.put(tag.seq(), tag); 
        }
        
        private void addWithMismatches(Tag tag, int mismatches) {
            Set<Tag> hammingTags = tag.generateTagsWithinHammingOf(mismatches);
            TagCategory category = tag.category();
            Map<String, Tag> tagsByCat = tags.get(category);
            for (Tag hammingTag : hammingTags) {
                warnIfAmbiguity(tagsByCat, hammingTag);
                tagsByCat.put(hammingTag.seq(), hammingTag);
            }                    
            
        }
        
        private void warnIfAmbiguity(Map<String, Tag> map, Tag barcode) {
            if (map.containsKey(barcode.seq())) {
                String otherName = map.get(barcode.seq()).name();
                if (!otherName.equals(barcode.name())) {
                    LOGGER.warning("Warning, ambiguity: " + barcode.name() +
                            " is too close in sequence to " + otherName);
                }
            }
        }

        //////////////////////////////////
        // Code for parsing tag layouts //
        //////////////////////////////////
        
        private void parseTagLayout(List<TagCategory> layout, String line) {
            if (!layout.isEmpty()) {
                warnValueAppearsTwice("LAYOUT");
                layout.clear();
            }
            
            String[] tags = line.split("= ")[1].split("\\|");

            try {
                int length = tags.length;
                for (int i = 0; i < length; i++) {
                    layout.add(TagCategory.valueOf(tags[i].toUpperCase()));
                }
            } catch (IllegalArgumentException ex) {
                LOGGER.log(Level.SEVERE, "Invalid tag category", ex);
                System.exit(0);
            }
        }
        
        
        private void parseSpacer(String line) {
            
            String[] fields = line.split(" ");
            if (fields.length != NUM_SPACER_LENGTH_FIELDS) {
                logConfigErrorAndExit(line);
            }
            
            int length = 0;
            
            try {
                length = Integer.parseInt(fields[2]);
            } catch (NumberFormatException ex) {
                logConfigErrorAndExit(ex, line);
            }
            
            if (length < 0) {
                logConfigErrorAndExit(line);
            }
            
            if (spacerLength != UNINITIALIZED) {
                warnValueAppearsTwice("SPACER");
            }
            
            spacerLength = length;
        }
        
        private void parseLaxity(String line) {
            
            String[] fields = line.split(" ");
            if (fields.length != NUM_LAXITY_FIELDS) {
                logConfigErrorAndExit(line);
            }
            
            int length = 0;
            
            try {
                length = Integer.parseInt(fields[2]);
            } catch (NumberFormatException ex) {
                logConfigErrorAndExit(ex, line);
            }
            
            if (length < 0) {
                logConfigErrorAndExit(line);
            }
            
            if (laxity != UNINITIALIZED) {
                warnValueAppearsTwice("LAXITY");
            }
            
            laxity = length;
        }
        
        //////////////////////////////////
        // Code for storing tag lengths //
        //////////////////////////////////
        
        private void populateLengths() {
            for (Entry<TagCategory, Map<String, Tag>> tagsByCategory : tags.entrySet()) {
                Set<Integer> lengthsSet = new TreeSet<>();
                TagCategory category = tagsByCategory.getKey();
                Map<String, Tag> tags = tagsByCategory.getValue();
                for (String tag : tags.keySet()) {
                    lengthsSet.add(tag.length());
                }

                // Convert the set of lengths to an array of (primitive) lengths.
                int[] lengthsArray = new int[lengthsSet.size()];
                int idx = 0;
                for (int len : lengthsSet) {
                    lengthsArray[idx++] = len;
                }
                
                lengths.put(category, lengthsArray);
            }
        }
        
        //////////////////////////
        // Misc logging methods //
        //////////////////////////
        
        private void logConfigErrorAndExit(Exception ex, String line) {
            LOGGER.log(Level.SEVERE, "Encountered error while parsing the " +
                    "configuration file. See line: " + line, ex);
            System.exit(0);
        }
        
        private void logConfigErrorAndExit(String line) {
            LOGGER.log(Level.SEVERE, "Encountered error while parsing the " +
                    "configuration file. See line: " + line);
            System.exit(0);
        }
        
        private void warnValueAppearsTwice(String var) {
            LOGGER.warning(var + " found twice in configuration file. " +
                    "Overwriting old value");
        }

    }
}