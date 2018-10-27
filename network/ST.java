package network;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeMap;

public class ST<Key extends Comparable<Key>, Value> implements Iterable<Key> {

    private TreeMap<Key, Value> st;

    /**
     * Initializes an empty symbol table.
     */
    public ST() {
        st = new TreeMap<Key, Value>();
        Collections.synchronizedSortedMap(st);
    }


    /**
     * Returns the value associated with the given key in this symbol table.
     *
     * @param  key the key
     * @return the value associated with the given key if the key is in this symbol table;
     *         <tt>null</tt> if the key is not in this symbol table
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public Set<Key> getKeys() {
    	return st.keySet();
    }
    
    
    public Value get(Key key) {
        if (key == null) throw new NullPointerException("called get() with null key");
        return st.get(key);
    }

    /**
     * Inserts the specified key-value pair into the symbol table, overwriting the old 
     * value with the new value if the symbol table already contains the specified key.
     * Deletes the specified key (and its associated value) from this symbol table
     * if the specified value is <tt>null</tt>.
     *
     * @param  key the key
     * @param  val the value
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public void put(Key key, Value val) {
        if (key == null) throw new NullPointerException("called put() with null key");
        if (val == null) st.remove(key);
        else             st.put(key, val);
    }

    /**
     * Removes the specified key and its associated value from this symbol table     
     * (if the key is in this symbol table).
     *
     * @param  key the key
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public void delete(Key key) {
        if (key == null) throw new NullPointerException("called delete() with null key");
        st.remove(key);
    }

    /**
     * Returns true if this symbol table contain the given key.
     *
     * @param  key the key
     * @return <tt>true</tt> if this symbol table contains <tt>key</tt> and
     *         <tt>false</tt> otherwise
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public boolean contains(Key key) {
        if (key == null) throw new NullPointerException("called contains() with null key");
        return st.containsKey(key);
    }

    /**
     * Returns the number of key-value pairs in this symbol table.
     *
     * @return the number of key-value pairs in this symbol table
     */
    public int size() {
        return st.size();
    }

    /**
     * Returns true if this symbol table is empty.
     *
     * @return <tt>true</tt> if this symbol table is empty and <tt>false</tt> otherwise
     */
    public boolean isEmpty() {
        return size() == 0;
    }

    /**
     * Returns all keys in this symbol table.
     * <p>
     * To iterate over all of the keys in the symbol table named <tt>st</tt>,
     * use the foreach notation: <tt>for (Key key : st.keys())</tt>.
     *
     * @return all keys in this symbol table
     */
    public Iterable<Key> keys() {
        return st.keySet();
    }

    /**
     * Returns all of the keys in this symbol table.
     * To iterate over all of the keys in a symbol table named <tt>st</tt>, use the
     * foreach notation: <tt>for (Key key : st)</tt>.
     * <p>
     * This method is provided for backward compatibility with the version from
     * <em>Introduction to Programming in Java: An Interdisciplinary Approach.</em>
     *
     * @return     an iterator to all of the keys in this symbol table
     * @deprecated Replaced by {@link #keys()}.
     */
    public Iterator<Key> iterator() {
        return st.keySet().iterator();
    }

    /**
     * Returns the smallest key in this symbol table.
     *
     * @return the smallest key in this symbol table
     * @throws NoSuchElementException if this symbol table is empty
     */
    public Key min() {
        if (isEmpty()) throw new NoSuchElementException("called min() with empty symbol table");
        return st.firstKey();
    }

    /**
     * Returns the largest key in this symbol table.
     *
     * @return the largest key in this symbol table
     * @throws NoSuchElementException if this symbol table is empty
     */
    public Key max() {
        if (isEmpty()) throw new NoSuchElementException("called max() with empty symbol table");
        return st.lastKey();
    }

    /**
     * Returns the smallest key in this symbol table greater than or equal to <tt>key</tt>.
     *
     * @param  key the key
     * @return the smallest key in this symbol table greater than or equal to <tt>key</tt>
     * @throws NoSuchElementException if there is no such key
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public Key ceiling(Key key) {
        if (key == null) throw new NullPointerException("called ceiling() with null key");
        Key k = st.ceilingKey(key);
        if (k == null) throw new NoSuchElementException("all keys are less than " + key);
        return k;
    }

    /**
     * Returns the largest key in this symbol table less than or equal to <tt>key</tt>.
     *
     * @param  key the key
     * @return the largest key in this symbol table less than or equal to <tt>key</tt>
     * @throws NoSuchElementException if there is no such key
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public Key floor(Key key) {
        if (key == null) throw new NullPointerException("called floor() with null key");
        Key k = st.floorKey(key);
        if (k == null) throw new NoSuchElementException("all keys are greater than " + key);
        return k;
    }

}

