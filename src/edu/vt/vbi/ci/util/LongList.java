
package edu.vt.vbi.ci.util;
import java.io.Serializable;
import java.nio.LongBuffer;




/**
 *  A class for storing and manipulating an
 *  ordered list of longs.
 */
public class LongList implements Serializable {
	static final long serialVersionUID = 10019855232l;

	private long[] longArray;
	private int   size;
	private int   initCapacity;

	private transient LongBuffer lb;

	/**
	 *  Creates a new LongList. The default initial
	 *  capacity is 10.
	 */
	public LongList() {

		size         = 0;
		initCapacity = 10;
		longArray     = new long[initCapacity];
	}

	/**
	 *  Creates a new LongList with the
	 *  specified initial capacity.
	 *
	 *  @param initialCapacity int - initial size of the internal array
	 */
	public LongList(int initialCapacity) {

		size         = 0;
		initCapacity = initialCapacity;
		longArray     = new long[initCapacity];
	}

	/**
	 *  Creates a new IntList with the specified
	 *  initial values.
	 *
	 *  @param initialList int[] - initial values for list
	 */
	public LongList(long[] initialList) {

		size         = 0;
		initCapacity = initialList.length;
		longArray     = new long[initCapacity];

		add(initialList);
	}

	/**
	 *  Adds a long to the end of the list.
	 *
	 *  @param long
	 */
	public final void add(long newLong) {

		if (size + 1 >= longArray.length) {
			int increase = 0;

			if (longArray.length < 2) {
				increase = 1;
			} else {
				increase = longArray.length / 2;
			}

			int   newCapacity = longArray.length + increase;
			long[] newLongArray = new long[newCapacity];

			System.arraycopy(longArray, 0, newLongArray, 0, longArray.length);

			longArray = newLongArray;
		}

		longArray[size] = newLong;

		size++;
	}

	/**
	 *  Adds the specified list of values
	 *  to the end of the list.
	 *
	 *  @param long[] new values
	 */
	public final void add(long[] newLongs) {

		for (int i = 0; i < newLongs.length; i++) {
			add(newLongs[i]);
		}
	}

	/**
	 *  Sets the value of the specified index
	 *  to value.
	 *
	 *  @param index int
	 *  @param value long
	 */
	public final void set(int index, long value) {

		if ((index > size) || (index < 0)) {
			return;
		}

		longArray[index] = value;
	}

	/**
	 *  Sets a range of values starting at index
	 *  to the specified list of values.
	 *
	 *  @param index int
	 *  @param values long[]
	 */
	public final void set(int index, long[] values) {

		if ((index + values.length) > size || (index < 0)) {
			return;
		}

		System.arraycopy(values, 0, longArray, index, values.length);
		for (int i = 0; i < values.length; i++) {
			longArray[index + i] = values[i];
		}
	}

	/**
	 *  Inserts the specified value at the specified
	 *  index.
	 *  @param index for new value
	 *  @param value
	 */
	public final void insert(int index, int value) {

		if ((index > size) || (index < 0)) {
			return;
		}
		if (size + 1 >= longArray.length) {
			int increase = 0;

			if ((longArray.length / 2) > (size + 1)) {
				increase = longArray.length / 2;
			} else {
				increase = size + 1;
			}

			int   newCapacity = longArray.length + increase;
			long[] newLongArray = new long[newCapacity];

			System.arraycopy(longArray, 0, newLongArray, 0, longArray.length);

			longArray = newLongArray;
		}

		//shift all elements down one
		for (int i = size; i > index; i--) {
			longArray[i] = longArray[i - 1];
		}

		longArray[index] = value;
		size++;
	}

	/**
	 *  Inserts the specified values beginning at
	 *  the specified index.
	 *  @param index for new value
	 *  @param values long[]
	 */
	public final void insert(int index, long[] values) {

		if ((index > size) || (index < 0)) {
			return;
		}

		int newCapacity = size + values.length;

		if (size + values.length >= longArray.length) {
			int increase = 0;

			if ((longArray.length / 2) > (size + values.length * 2)) {
				increase = longArray.length / 2;
			} else {
				increase = size + values.length * 2;
			}

			newCapacity = longArray.length + increase;
		}

		long[] newLongArray = new long[newCapacity];

		System.arraycopy(longArray, 0, newLongArray, 0, index);
		System.arraycopy(values, 0, newLongArray, index, values.length);
		System.arraycopy(longArray, index, newLongArray, index + values.length,
				size - index);

		longArray = newLongArray;
		size     += values.length;
	}

	/**
	 *  Gets the long at the specified index.
	 *
	 *  @param index
	 *  @return long value at index
	 *  @throws ArrayIndexOutOfBoundsException
	 */
	public final long get(int index) {

		long r;

		if ((index >= size) || (index < 0)) {
			throw new ArrayIndexOutOfBoundsException(index);
		} else {
			r = longArray[index];
		}

		return r;
	}

	/**
	 *  Returns the range of values from index
	 *  through (index + length - 1).
	 *
	 *  @param index int - starting index
	 *  @param length int - length of range
	 *  @return long[]
	 */
	public final long[] get(int index, int length) {

		if ((index < 0) || (index + length) > size) {
			return null;
		}

		long[] r = new long[length];

		System.arraycopy(longArray, index, r, 0, length);
		return r;
	}

	/**
	 *  Removes the value at the specified index and compresses
	 *  the internal array.
	 *
	 *  @param index to be removed
	 *  @return long value that was removed
	 */
	public final long remove(int index) {

		//System.out.println("Entering IntList.remove(int index)");
		long r;

		if ((index > size) || (index < 0)) {
			r = -1;
		} else {
			r = longArray[index];

			for (int i = index; i < size - 1; i++) {
				longArray[i] = longArray[i + 1];
			}

			size--;
		}

		return r;
	}

	/**
	 *  Removes length number of elements,
	 *  begining at index.
	 *
	 *  @param int index
	 *  @param int length
	 *  @return long[] the list of removed ints
	 */
	public final long[] remove(int index, int length) {

		if ((index + length) > longArray.length) {
			length = longArray.length - index;
		}

		//make sure index is valid
		if ((index >= longArray.length) || (index < 0)) {
			index  = 0;
			length = 0;
		}

		long[] r = new long[length];

		System.arraycopy(longArray, index, r, 0, length);

		for (int i = index; i < size - length; i++) {
			longArray[i] = longArray[i + length];
		}

		size -= length;

		if (size < 0) {
			size = 0;
		}

		return r;
	}

	/**
	 *  Returns the size of the list.
	 *
	 *  @return int number of elements in list
	 */
	public final int size() {
		return size;
	}

	/**
	 *  Returns the list.
	 *
	 *  @return long[]
	 */
	public final long[] getLongs() {

		//Change made by nishant on 08-01-04
		// to avoid exception occuring after removing 
		//entire sequence under particular circumstances.
		if (size < 0) {
			size = 0;
		}

		long[] returnArray = new long[size];

		System.arraycopy(longArray, 0, returnArray, 0, size);

		return returnArray;
	}

	/**
	 *  Returns the first index of the specified value.
	 *
	 *  @param value long
	 *  @return int index of value, or -1 if value is not in list
	 */
	public final int firstIndexOf(long value) {

		for (int i = 0; i < size; i++) {
			if (longArray[i] == value) {
				return i;
			}
		}

		return -1;
	}

	/**
	 *  Returns the last index of the specified value.
	 *
	 *  @param value int
	 *  @return int index of value, or -1 if value is not in list
	 */
	public final int lastIndexOf(int value) {

		for (int i = size - 1; i >= 0; i--) {
			if (longArray[i] == value) {
				return i;
			}
		}

		return -1;
	}

	/**
	 *  Returns every index of the specified value.
	 *
	 *  @param value int
	 *  @return int[] indices of value - empty if value is not found
	 */
	public final int[] everyIndexOf(int value) {
		IntList r = new IntList();

		for (int i = 0; i < size; i++) {
			if (longArray[i] == value) {
				r.add(i);
			}
		}

		return r.getInts();
	}

	/**
	 *  Resizes the internal array to exactly fit number of elements.
	 */
	public final void trimToSize() {
		long[] newArray = new long[size];
		System.arraycopy(longArray, 0, newArray, 0, size);
		longArray = newArray;
	}

	/**
	 *  Clears the list - removes all elements.
	 *  After clearing, the IntList will have size = 0
	 *  and capacity = initial capacity.
	 */
	public final void clear() {
		longArray = new long[initCapacity];
		size     = 0;
	}   
}