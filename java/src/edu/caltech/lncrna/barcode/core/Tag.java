package edu.caltech.lncrna.barcode.core;

import java.util.HashSet;
import java.util.Set;

import edu.caltech.lncrna.barcode.core.TagCategory;

public class Tag {
    private final TagCategory cat;
    private final String name;
    private final String seq;
    
    public Tag(TagCategory cat, String name, String seq) {
        this.cat = cat;
        this.name = name;
        this.seq = seq;
    }
    
    public TagCategory category() {
        return cat;
    }
    
    public String name() {
        return name;
    }
    
    public String seq() {
        return seq;
    }
    
    /**
     * Generates the set of barcodes between 0 and <code>k</code>
     * Hamming-distance (inclusive) from this barcode.
     * @param k - the maximum Hamming-distance
     */
    public Set<Tag> generateTagsWithinHammingOf(int k) {
        Set<Tag> rtrn = new HashSet<Tag>();
        rtrn.add(this);
        
        for (int i = 0; i < k; i++) {
            rtrn = generateTagsWithinHammingOfOne(rtrn);
        }
        
        return rtrn;
    }
    
    /**
     * Generates the set of barcodes 0 or 1 Hamming-distance from
     * any barcode within the given set.
     * @param barcodes - the set of barcodes
     */
    private Set<Tag> generateTagsWithinHammingOfOne(Set<Tag> tags) {
        Set<Tag> rtrn = new HashSet<Tag>();
        
        for (Tag tag : tags) {
            rtrn.add(tag);
            String sequence = tag.seq();
            for (int i = 0; i < sequence.length(); i++) {
                String left = sequence.substring(0, i);
                String right = sequence.substring(i + 1);
                
                rtrn.add(new Tag(tag.category(),
                                     tag.name,
                                     left + "A" + right));
                rtrn.add(new Tag(tag.category(),
                                     tag.name,
                                     left + "C" + right));
                rtrn.add(new Tag(tag.category(),
                                     tag.name,
                                     left + "G" + right));
                rtrn.add(new Tag(tag.category(),
                                     tag.name,
                                     left + "T" + right));
                rtrn.add(new Tag(tag.category(),
                                     tag.name,
                                     left + "N" + right));
            }
        }
        return rtrn;
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        
        if (!(o instanceof Tag)) {
            return false;
        }
        
        Tag other = (Tag) o;
        
        return cat.equals(other.cat) &&
               name.equals(other.name) &&
               seq.equals(other.seq);
    }
    
    @Override
    public int hashCode() {
        int hashCode = 17;
        hashCode = 37 * hashCode + name.hashCode();
        hashCode = 37 * hashCode + cat.hashCode();
        hashCode = 37 * hashCode + seq.hashCode();
        return hashCode;
    }
}