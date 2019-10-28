package edu.caltech.lncrna.barcode.core;

/**
 * An enumeration of what type of molecule a read could have originated from.
 * <p>
 * For example, <code>Origin.RNA</code> indicates that a read came from an RNA
 * molecule.
 */
public enum Origin {
    RNA(), DNA(), AMBIGUOUS();
}