package eu.exahype.fromtosable;

/**
 * A quick naive way to store structured JSON-like information
 * with absolute paths.
 
 * -> There is LinkedHashMap as ordereddict in java to use.
 **/
public interface Structured {
  public void put(String key, String value);
  public void put(String key, int value);
  public void put(String key, double value);
  public void put(String key, boolean value);
}
