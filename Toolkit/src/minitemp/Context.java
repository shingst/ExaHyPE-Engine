package minitemp;

import javax.script.ScriptEngineManager;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import java.util.Collection;

/**
 * Context class
 *
 * Store the pair of key->value and use it to evaluate expressions.
 * Type allowed for the value: String, boolean, int, double
 * A key can have only one value, redefining it (even with a value of a different type) override the
 * previous value
 *
 * Use put(key, value) to add a key and its value to the context.
 * evaluateX evaulate an expression and return a value of type X.
 */
public class Context {
  
  /** Use a Javascript engine to store the context and evaluate expressions */
  private ScriptEngine engine;
  /** 
   * If corrupted this context cannot perform any evaluation
   * Used to avoid throwing exceptions to handle during put operations
   */
  private boolean corrupted = false;
  /** epsilon = 10^-11, use find int returned as double */
  private static final double EPS = 0.00000000001;
  
  public Context() {    
    engine = (new ScriptEngineManager()).getEngineByName("JavaScript");
  }
  
  /** 
   * Put a key->value pair in the context 
   * Errors are catched but trigger warnings and corrupt the Context (=unusable to render)
   */
  public void put(String key, Object value) {
    try {
      validateKey(key);
      engine.put(key, value); 
    } catch (Exception e) {
      System.err.println("ERROR in Context: put("+key+", "+value+") failed. Context is corrupted.");
      corrupted = true;
    }
  }
  
  /** Put only if no value is set for this key */
  public void putIfNotSet(String key, Object value) {
    if(getValueOrNull(key) == null) {
      put(key, value);
    }
  }
  
  /** Verify that a key is one word with normal characters ([a-zA-Z_0-9]+) */
  private void validateKey(String key) throws IllegalArgumentException {
    if(key == null || !key.matches("^\\w+$")) {
      System.err.println("ERROR in Context: invalid key: "+key);
      throw new IllegalArgumentException("Invalid key: "+key);
    }
  }
  
  private void checkCorruption() throws IllegalArgumentException {
    if(corrupted) {
      throw new IllegalArgumentException("Cannot evaluate corrupted context");
    }
  }
  
  /** Get a String from the context, null => "" */
  public String evaluateString(String expression) throws IllegalArgumentException { 
    checkCorruption();
    try {
      Object valueRaw = engine.eval(expression);
      if(valueRaw == null) {
        return "";
      }
      if(valueRaw instanceof Double) { //need to detect int returned as double to not print .0
         Double value = (Double)valueRaw;
        if(Math.abs(Math.floor(value)-value) < EPS && Math.abs(value) < Integer.MAX_VALUE) {
          return Integer.toString(value.intValue()); //return double as an int
        } else {
          return value.toString(); 
        }
      }
      return valueRaw.toString();
    } catch (ScriptException e) {
      throw new IllegalArgumentException("Cannot evaluate token's expression: "+expression);
    }
  }
  
  /** Get a boolean from the context */
  public boolean evaluateBoolean(String expression) throws IllegalArgumentException { 
    checkCorruption();
    try {
      Object valueRaw = engine.eval(expression);
      if(valueRaw instanceof Boolean) {
        return ((Boolean)valueRaw).booleanValue();
      }
    } catch (ScriptException e) {} //throw an exception later reach that point
    throw new IllegalArgumentException("Cannot evaluate token's expression as boolean: "+expression);
  }
  
  /** Get a collection from the context */
  @SuppressWarnings("unchecked") //cast to Collection<Object> is unchecked but normal
  public Collection<Object> getCollection(String key) throws IllegalArgumentException {
    checkCorruption();
    try {
      Object valueRaw = engine.get(key);
      if(valueRaw instanceof Collection) {
        return (Collection<Object>)valueRaw;
      }
    } catch (IllegalArgumentException e) {
      throw new IllegalArgumentException("Cannot find a value to the key "+key);
    }
    throw new IllegalArgumentException("The key "+key+" doesn't match to a collection");
  }
  
  public Object getValueOrNull(String key) {
    checkCorruption();
    try {
      return engine.get(key);
    } catch (IllegalArgumentException e) {
      return null; //key don't exist
    }
  }
}
