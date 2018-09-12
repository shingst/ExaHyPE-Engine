package eu.exahype.fromtosable;

/**
 * This class represents a dictionary entry.
 *
 *
 * A quick naive way to store structured JSON-like information
 * with absolute paths.
 *
 * -> There is LinkedHashMap as ordereddict in java to use.
 **/
public interface Structured {

  /** Make a new one from this type */
  public Structured makeInstance();
  
  /** Store anything right hand side typed */
  public void put(Object value);
  
  /**
   * Declare this pointing as structure, get access to child.
   * This is a shorthand for:
   *
   *    Structured rhs = thing.makeInstance();
   *    think.put(rhs);
   */
  public Structured descend(String subkey);

  /** Shall return only the path in a nice manner such as "/a/b/c" */
  public String toString();
}
