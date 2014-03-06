package edu.vt.vbi.ci.util;

import java.util.BitSet;

public class ExtendedBitSet extends BitSet implements Comparable {

	public ExtendedBitSet getAnd(BitSet bs) {
		ExtendedBitSet r = (ExtendedBitSet) this.clone();
		r.and(bs);
		return r;
	}
	
	public ExtendedBitSet getOr(BitSet bs) {
		ExtendedBitSet r = (ExtendedBitSet) this.clone();
		r.or(bs);
		return r;
	}
	
	public ExtendedBitSet getXor(BitSet bs) {
		ExtendedBitSet r = (ExtendedBitSet) this.clone();
		r.xor(bs);
		return r;
	}
	
	public ExtendedBitSet getComplement(int length) {
		ExtendedBitSet r = (ExtendedBitSet)this.clone();
		r.flip(0, length);
		return r;
	}
	
	public ExtendedBitSet getAndNot(BitSet bs) {
		ExtendedBitSet r = (ExtendedBitSet)this.clone();
		r.andNot(bs);
		return r;
	}
	
	public int[] getSetValues() {
		int[] r = new int[cardinality()];
		int index = 0;
		for(int set = nextSetBit(0); set >= 0; set = nextSetBit(set+1)) {
			r[index] = set;
			index++;
		}
		
		return r;
	}
	
	public String getString(int size) {
		String r = null;
			char present = '*';
			char absent = '.';
			
			StringBuffer sb = new StringBuffer();
			for(int i = 0; i < size; i++) {
				if(this.get(i)) {
					sb.append(present);
				} else {
					sb .append(absent);
				}
			}
			r = sb.toString();
		return r;
	}

	public int compareTo(Object o) {
		int r = 0;
		ExtendedBitSet obs = (ExtendedBitSet)o;
		r = this.cardinality() - obs.cardinality();
		return r;
	}

	
}
