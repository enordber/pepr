package edu.vt.vbi.ci.util;

public class JSONUtilities {
//
//	private static final char openBracket = '[';
//	private static final char closeBracket = ']';
//	private static final char openCurly = '{';
//	private static final char closeCurly = '}';
//	private static final char colon = ':';
//	private static final char comma = ',';
//	private static final char quote = '"';
//
//	public static String getJSONString(double[] d) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(d != null && d.length > 0) {
//			sb.append(d[0]);
//			for(int i = 1; i < d.length; i++) {
//				sb.append(comma);
//				sb.append(d[i]);
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//	
//	public static String getJSONString(double[][] d) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(d != null && d.length > 0) {
//			sb.append(getJSONString(d[0]));
//			for(int i = 1; i < d.length; i++) {
//				sb.append(comma);
//				sb.append(getJSONString(d[i]));
//			}
//		}
//		sb.append(closeBracket);
//		return r;
//	}
//
//	public static String getJSONString(int[] a) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(a != null && a.length > 0) {
//			sb.append(a[0]);
//			for(int i = 1; i < a.length; i++) {
//				sb.append(comma);
//				sb.append(a[i]);
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//	
//	public static String getJSONString(int[][] a) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(a != null && a.length > 0) {
//			sb.append(getJSONString(a[0]));
//			for(int i = 1; i < a.length; i++) {
//				sb.append(comma);
//				sb.append(getJSONString(a[i]));
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//	
//	public static String getJSONString(String[] a) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(a != null && a.length > 0) {
//			sb.append(quote);
//			sb.append(a[0]);
//			sb.append(quote);
//			for(int i = 1; i < a.length; i++) {
//				sb.append(comma);
//				sb.append(quote);
//				sb.append(a[i]);
//				sb.append(quote);
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//
//	public static String getJSONString(float[][] f) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(f != null && f.length > 0) {
//			sb.append(getJSONString(f[0]));
//			for(int i = 1; i < f.length; i++) {
//				sb.append(comma);
//				sb.append(getJSONString(f[i]));
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//
//	public static String getJSONString(float[] f) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openBracket);
//		if(f != null && f.length > 0) {
//			sb.append(quote);
//			sb.append(f[0]);
//			sb.append(quote);
//			for(int i = 1; i < f.length; i++) {
//				sb.append(comma);
//				sb.append(quote);
//				sb.append(f[i]);
//				sb.append(quote);
//			}
//		}
//		sb.append(closeBracket);
//		r = sb.toString();
//		return r;
//	}
//
//	/**
//	 * names and values arrays must be the same length.
//	 * 
//	 * @param names
//	 * @param values
//	 * @return
//	 */
//	public static String assembleJSONObject(String[] names, String[] values) {
//		String r = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(openCurly);
//		if(names != null && names.length > 0) {
//			sb.append(quote);
//			sb.append(names[0]);
//			sb.append(quote);
//			sb.append(colon);
//			sb.append(values[0]);
//			for(int i = 1; i < names.length; i++) {
//				sb.append(comma);
//				sb.append(quote);
//				sb.append(names[i]);
//				sb.append(quote);
//				sb.append(colon);
//				sb.append(values[i]);				
//			}
//		}
//		sb.append(closeCurly);
//		r = sb.toString();
//		return r;
//	}
//
}
