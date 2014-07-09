package edu.vt.vbi.ci.util;
import java.io.Serializable;
import java.util.Arrays;

/**
 *  A class for storing and manipulating an
 *  ordered list of ints.
 */
public class IntList implements Serializable {
	static final long serialVersionUID = 10019855261l;

	private int[] intArray;
	private int   size;
	private int   initCapacity;

	//when sorted, search operations can be sped up using Arrays.binarySearch()
	private boolean sorted = false;

	/**
	 *  Creates a new IntList. The default initial
	 *  capacity is 10.
	 */
	public IntList() {

		size         = 0;
		initCapacity = 10;
		intArray     = new int[initCapacity];
	}

	/**
	 *  Creates a new IntList with the
	 *  specified initial capacity.
	 *
	 *  @param initialCapacity int - initial size of the internal array
	 */
	public IntList(int initialCapacity) {

		size         = 0;
		initCapacity = initialCapacity;
		intArray     = new int[initCapacity];
	}

	/**
	 *  Creates a new IntList with the specified
	 *  initial values.
	 *
	 *  @param initialList int[] - initial values for list
	 */
	public IntList(int[] initialList) {

		size         = 0;
		initCapacity = initialList.length;
		intArray     = new int[initCapacity];

		add(initialList);
	}

	/**
	 *  Adds an int to the end of the list.
	 *
	 *  @param int
	 */
	public final void add(int newInt) {

		if (size + 1 >= intArray.length) {
			int increase = 0;

			if (intArray.length < 2) {
				increase = 1;
			} else {
				increase = intArray.length / 2;
			}

			int   newCapacity = intArray.length + increase;
			int[] newIntArray = new int[newCapacity];

			System.arraycopy(intArray, 0, newIntArray, 0, intArray.length);

			intArray = newIntArray;
		}

		intArray[size] = newInt;

		size++;
		sorted = false;
	}

	/**
	 *  Adds the specified list of values
	 *  to the end of the list.
	 *
	 *  @param int[] new values
	 */
	public final void add(int[] newInts) {

		ensureCapacity(size+newInts.length);
		for (int i = 0; i < newInts.length; i++) {
			add(newInts[i]);
		}
		sorted = false;
	}

	/**
	 * Makes sure the underlying int[] is at least large enough
	 * to hold the specified number of ints. 
	 * @param capacity
	 */
	public void ensureCapacity(int capacity) {
		if(intArray.length < capacity) {
			int[] newArray = new int[capacity];
			System.arraycopy(intArray, 0, newArray, 0, intArray.length);
			intArray = newArray;
		}
	}

	/**
	 *  Sets the value of the specified index
	 *  to value.
	 *
	 *  @param index int
	 *  @param value int
	 */
	public final void set(int index, int value) {

		if ((index > size) || (index < 0)) {
			return;
		}

		intArray[index] = value;
		sorted = false;
	}

	/**
	 *  Sets a range of values starting at index
	 *  to the specified list of values.
	 *
	 *  @param index int
	 *  @param values int[]
	 */
	public final void set(int index, int[] values) {

		if ((index + values.length) > size || (index < 0)) {
			return;
		}

		System.arraycopy(values, 0, intArray, index, values.length);
		for (int i = 0; i < values.length; i++) {
			intArray[index + i] = values[i];
		}
		sorted = false;
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

		if (size + 1 >= intArray.length) {
			int increase = 0;

			if ((intArray.length / 2) > (size + 1)) {
				increase = intArray.length / 2;
			} else {
				increase = size + 1;
			}

			int   newCapacity = intArray.length + increase;
			int[] newIntArray = new int[newCapacity];

			System.arraycopy(intArray, 0, newIntArray, 0, intArray.length);

			intArray = newIntArray;
		}

		//shift all elements down one
		for (int i = size; i > index; i--) {
			intArray[i] = intArray[i - 1];
		}

		intArray[index] = value;

		size++;

		sorted = false;
	}

	/**
	 *  Inserts the specified values beginning at
	 *  the specified index.
	 *  @param index for new value
	 *  @param values int[]
	 */
	public final void insert(int index, int[] values) {

		if ((index > size) || (index < 0)) {
			return;
		}

		int newCapacity = size + values.length;

		if (size + values.length >= intArray.length) {
			int increase = 0;

			if ((intArray.length / 2) > (size + values.length * 2)) {
				increase = intArray.length / 2;
			} else {
				increase = size + values.length * 2;
			}

			newCapacity = intArray.length + increase;
		}

		int[] newIntArray = new int[newCapacity];

		System.arraycopy(intArray, 0, newIntArray, 0, index);
		System.arraycopy(values, 0, newIntArray, index, values.length);
		System.arraycopy(intArray, index, newIntArray, index + values.length,
				size - index);

		intArray = newIntArray;
		size     += values.length;
		sorted = false;
	}

	/**
	 *  Gets the int at the specified index.
	 *
	 *  @param index
	 *  @return int value at index
	 *  @throws ArrayIndexOutOfBoundsException
	 */
	public final int get(int index) {

		int r;

		if ((index >= size) || (index < 0)) {
			throw new ArrayIndexOutOfBoundsException(index);
		} else {
			r = intArray[index];
		}

		return r;
	}

	/**
	 *  Returns the range of values from index
	 *  through (index + length - 1).
	 *
	 *  @param index int - starting index
	 *  @param length int - length of range
	 *  @return int[]
	 */
	public final int[] get(int index, int length) {

		if ((index < 0) || (index + length) > size) {
			return null;
		}

		int[] r = new int[length];

		System.arraycopy(intArray, index, r, 0, length);

		return r;
	}

	/**
	 *  Removes the value at the specified index and compresses
	 *  the internal array.
	 *
	 *  @param index to be removed
	 *  @return int value that was removed
	 */
	public final int remove(int index) {

		//System.out.println("Entering IntList.remove(int index)");
		int r;

		if ((index > size) || (index < 0)) {
			r = -1;
		} else {
			r = intArray[index];

			for (int i = index; i < size - 1; i++) {
				intArray[i] = intArray[i + 1];
			}

			size--;
			//when sorted, binary search is used for firstIndexOf(). 
			//This requires the actual int[] to be kept in sorted order,
			//otherwise the binary search may not work. To make sure this
			//	happens, any values after the end of the list should be set to 
			//Integer.MAX_VALUE
			intArray[size] = Integer.MAX_VALUE;
		}

		return r;
	}

	/**
	 *  Removes length number of elements,
	 *  begining at index.
	 *
	 *  @param int index
	 *  @param int length
	 *  @return int[] the list of removed ints
	 */
	public final int[] remove(int index, int length) {

		//System.out.println("IntList.remove(int index, int length)");
		if ((index + length) > intArray.length) {
			length = intArray.length - index;
		}

		//make sure index is valid
		if ((index >= intArray.length) || (index < 0)) {
			index  = 0;
			length = 0;
		}

		int[] r = new int[length];

		System.arraycopy(intArray, index, r, 0, length);

		for (int i = index; i < size - length; i++) {
			intArray[i] = intArray[i + length];
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
	 *  @return int[]
	 */
	public final int[] getInts() {

		//Change made by nishant on 08-01-04
		// to avoid exception occuring after removing 
		//entire sequence under particular circumstances.
		if (size < 0) {
			size = 0;
		}

		int[] returnArray = new int[size];

		System.arraycopy(intArray, 0, returnArray, 0, size);

		return returnArray;
	}

	/**
	 *  Returns the first index of the specified value.
	 *
	 *  @param value int
	 *  @return int index of value, or -1 if value is not in list
	 */
	public final int firstIndexOf(int value) {
		int r = -1;

		if(sorted) {
			r = Arrays.binarySearch(intArray, value);
		} else {
			for (int i = 0; i < size && r < 0 ; i++) {
				if (intArray[i] == value) {
					r = i;
				}
			}
		}
		return r;
	}

	/**
	 *  Returns the last index of the specified value.
	 *
	 *  @param value int
	 *  @return int index of value, or -1 if value is not in list
	 */
	public final int lastIndexOf(int value) {

		for (int i = size - 1; i >= 0; i--) {
			if (intArray[i] == value) {
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
			if (intArray[i] == value) {
				r.add(i);
			}
		}

		return r.getInts();
	}

	/**
	 *  Resizes the internal array to exactly fit number of elements.
	 */
	public final void trimToSize() {

		if(intArray.length > size) {
			int[] newArray = new int[size];

			System.arraycopy(intArray, 0, newArray, 0, size);

			intArray = newArray;
		}
	}

	/**
	 * Sorts the ints being stored. This allows faster searches. This also
	 * trims the underlying int[], to avoid problems with binary search.
	 * Any add(), set(), or insert() calls will invalidate the sorting, and
	 * sort() must be called again afterwards (even if you KNOW that the
	 * insert(), set(), or add() did not alter the sorting). Currently, 
	 * only firstIndexOf() takes advantage of sorting.
	 */
	public void sort() {
		if(!sorted) { 
			trimToSize();
			Arrays.sort(intArray);
			sorted = true;
		}
	}

	/**
	 *  Clears the list - removes all elements.
	 *  After clearing, the IntList will have size = 0
	 *  and capacity = initial capacity.
	 */
	public final void clear() {
		intArray = new int[initCapacity];
		size     = 0;
	}
}
