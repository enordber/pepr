package edu.vt.vbi.ci.util;

public class StringPair {

	private String a;
	private String b;
	private int hc;

	public StringPair(String a, String b) {
		if(a.compareTo(b) > 0) {
			this.a = a;
			this.b = b;
		} else {
			this.a = b;
			this.b = a;
		}

		this.hc = this.a.hashCode() * 10 + this.b.hashCode();
	}

	public int hashCode() {
		return hc;
	}

	public boolean equals(Object o) {
		boolean r = false;
		try {
			StringPair osp = (StringPair)o;
			r = osp.a.equals(this.a) && osp.b.equals(this.b);
		} catch(ClassCastException cce) {}
		
		return r;
	}
	
	public String toString() {
		String r = "{" + a + "," + b + "}";
		return r;
	}

	public String getA() {
		return a;
	}
	
	public String getB() {
		return b;
	}
}
