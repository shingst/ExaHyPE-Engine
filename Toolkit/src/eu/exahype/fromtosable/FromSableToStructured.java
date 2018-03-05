package eu.exahype.fromtosable;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.*;

import eu.exahype.node.*;

import eu.exahype.variables.Variables;

import java.lang.reflect.*;

import java.util.regex.Pattern;
import java.util.regex.Matcher;


/**
 * A SableCC visitor to transform a specfile into a structured
 * container.
 *
 * Which subsequently can be exported to a format such as JSON.
 *
 * This can be helpful to process specfiles without rough
 * regexpes or similiar.
 *
 **/
public class FromSableToStructured  {
  /// This maps an [hierarchical path key] -> [value]
  // private java.util.Map<String, String> _data;
  Structured _data;
  
  /// Wheter to include null references to optional members 
  boolean _includeMissingOptionals;
  
  public FromSableToStructured() {
    //_data = new java.util.HashMap<String,String>();
    _data = new StructuredDumper();
  }

  /// Chainable setter (argument factory like)
  public FromSableToStructured setIncludeMissingOptionals(boolean includeMissingOptionals) {
    _includeMissingOptionals = includeMissingOptionals;
    return this;
  }
  
  public void dump(Node root) {
    dump(root, "");
  }
  
  /**
   * A recursive serialzer for a SableCC Node tree (AST).
   *
   *  obj: A Node subclass from eu.exahype.node.*
   *  path: A path such as /foo/bar/baz
   **/
  public void dump(Object obj, String path) {
    if(obj == null) {
      // the null/None type
      _data.put(path, "<None>");
      return;
    }
    String qualifiedname = obj.getClass().getName(); // w namespace
    String clsname = obj.getClass().getSimpleName(); // w.o. namespace
    boolean sableNode = qualifiedname.startsWith("eu.exahype.node");
    if(sableNode && clsname.startsWith("A") || clsname.equals("Start")) {
      Field[] allFields = obj.getClass().getDeclaredFields();
      for (Field field : allFields) {
	  field.setAccessible(true);
          Pattern p = Pattern.compile("^_(.*)_$");
          Matcher m = p.matcher(field.getName());
          if (m.find()) {
          //if(field.getName().matches("_.*_")) {
             // this is one of the sableCC guys
             String fieldName = m.group(1);
             try {
	       Object target = field.get(obj);
               if(target != null || _includeMissingOptionals) {
                  dump(target, path + "/" + fieldName);
	       } else {
	          System.out.println(path + ": Skipping field " + field.getName());
	       }
	     } catch(IllegalAccessException e) {
	       System.out.println(e.toString());
	     }
          } else {
             System.out.println(path + ": Field "+field.toGenericString()+" not interesting");
          }
       }
    } else if(sableNode && clsname.startsWith("P")) {
      // Pointer: Look it up
      try {
        Class aclass = Class.forName("A" + clsname.substring(1, clsname.length()));
        dump(aclass.cast(obj), path);
      } catch(ClassNotFoundException e) {
        System.out.println(e.toString());
      }
    } else if(sableNode && clsname.startsWith("T")) {
      // Token: Finalize it
      _data.put(path, ((Token) obj).getText());
    } else if (obj instanceof Collection<?>){
      // it is something like LinkedList<PSolver>
      int i = 0;
      for(Object child : (Collection)obj) {
         dump(child, path+"["+Integer.toString(i)+"]");
         i++;
      }
    } else if(clsname.equals("EOF")) {
      // the SableCC end node.
    } else {
      throw new RuntimeException("Unknown class: " + clsname);
    }
  }
  
  public void mandatory(String target, Token token) {
    _data.put(target, token.getText());
  }
  
  public boolean optional(String target, Token token) {
    if(token != null)
      _data.put(target, token.getText());
    return (token != null);
  }
}
