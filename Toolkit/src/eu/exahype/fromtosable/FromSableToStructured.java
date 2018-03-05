package eu.exahype.fromtosable;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.*;

import eu.exahype.node.*;

import eu.exahype.variables.Variables;

import java.lang.reflect.*;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

import java.io.OutputStream;

import org.json.*;


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
//    _data = new StructuredDumper();
  }

  /// Chainable setter (argument factory like)
  public FromSableToStructured setIncludeMissingOptionals(boolean includeMissingOptionals) {
    _includeMissingOptionals = includeMissingOptionals;
    return this;
  }
  
  public String toJSON(Node document) {
    try{
      int indentFactor = 4;
      return (dumpDocument((Start)document)).toString(indentFactor);
    } catch(JSONException e) {
      System.err.println("Could not create JSON: " + e.toString());
      System.exit(-1);
      return "";
    }
  }
  
  /**
   * Jump already to the main sub node in the head node.
   * This sugar makes use of exahype.grammar knowledge. If you want to
   * be generic, just use dump() with any Node.
   **/
  public JSONObject dumpDocument(Start document) {
    return (JSONObject) dump(document.getPProject());
  }
  
  /**
   * Try to make sense of a Token by exploiting knowledge of the exahype.grammar.
   *
   **/
  public Object parseToken(Token t) {
    String s = t.getText();
    String clsname = t.getClass().getSimpleName(); // w.o. namespace
    
    // ExaHyPE booleans
    if(clsname.equals("TTokenOnOff")) {
      Boolean truth = s.toLowerCase().equals("on");
      return truth;
    }
    
    try {
      Integer i = Integer.parseInt(s);
      return i;
    } catch(NumberFormatException nfe) {}
    
    try {
      Double d = Double.parseDouble(s);
      return d;
    } catch(NumberFormatException nfe) {}
    
    return s;
  }
  
  /**
   * A recursive serialzer for a SableCC Node tree (AST).
   *
   *  obj: A Node subclass from eu.exahype.node.*
   *  path: A path such as /foo/bar/baz
   **/
  public Object dump(Object obj) {
    if(obj == null) {
      // the null/None type
      //path.put("<None>");
      return JSONObject.NULL;
    }
    String qualifiedname = obj.getClass().getName(); // w namespace
    String clsname = obj.getClass().getSimpleName(); // w.o. namespace
    boolean sableNode = qualifiedname.startsWith("eu.exahype.node");

    if(sableNode && clsname.startsWith("A") || clsname.equals("Start")) {
      // A node can be mapped to a structure (dictionary)
      Field[] allFields = obj.getClass().getDeclaredFields();
      JSONObject out_obj = new JSONObject();
      for (Field field : allFields) {
	  field.setAccessible(true);
          Pattern p = Pattern.compile("^_(.*)_$");
          Matcher m = p.matcher(field.getName());
          if (m.find()) {
             // this is one of the sableCC guys
             String fieldName = m.group(1);
             try {
	       Object target = field.get(obj);
               if(target != null || _includeMissingOptionals) {
                  out_obj.put(fieldName, dump(target));
	       } else {
	          //System.out.println(path.pathtoString() + ": Skipping field " + field.getName());
	       }
	     } catch(IllegalAccessException e) {
	       System.out.println(e.toString());
	     }
          } else {
             //System.out.println(path.pathtoString() + ": Field "+field.toGenericString()+" not interesting");
          }
       }
       return out_obj;
    } else if(sableNode && clsname.startsWith("P")) {
      // Pointer: Look up and pass throught
      try {
        Class aclass = Class.forName("A" + clsname.substring(1, clsname.length()));
        return dump(aclass.cast(obj));
      } catch(ClassNotFoundException e) {
        throw new RuntimeException("SableCC inconsistent", e);
      }
    } else if(sableNode && clsname.startsWith("T")) {
      // Token: Finally a primitive value
      return parseToken((Token) obj);
    } else if (obj instanceof Collection<?>){
      // a list (such as LinkedList<PSolver>) can be mapped to an array
      //Structured lst = path.put_list(path);
      JSONArray out_lst = new JSONArray();
      for(Object child : (Collection)obj) {
         out_lst.put(dump(child));
      }
      return out_lst;
    } else if(clsname.equals("EOF")) {
      // the SableCC end node.
      return JSONObject.NULL;
    } else {
      throw new RuntimeException("Unknown class: " + clsname);
    }
  }
}
