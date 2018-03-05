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

import javax.json.*;


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
  
  JsonBuilderFactory factory;
  
  public FromSableToStructured() {
    //_data = new java.util.HashMap<String,String>();
//    _data = new StructuredDumper();

     Map<String, Object> config = new HashMap<String, Object>();
     // config.put(JsonGenerator.PRETTY_PRINTING, true);
     factory = Json.createBuilderFactory(config);
  }

  /// Chainable setter (argument factory like)
  public FromSableToStructured setIncludeMissingOptionals(boolean includeMissingOptionals) {
    _includeMissingOptionals = includeMissingOptionals;
    return this;
  }
  
  public void printJson(Node document, OutputStream out) {
    JsonWriter writer = Json.createWriter(out);
    JsonObject data = (JsonObject) dump(document);
    writer.writeObject(data);
    writer.close();	
  }
  
  /**
   * A recursive serialzer for a SableCC Node tree (AST).
   *
   *  obj: A Node subclass from eu.exahype.node.*
   *  path: A path such as /foo/bar/baz
   **/
  public JsonValue dump(Object obj) {
    if(obj == null) {
      // the null/None type
      //path.put("<None>");
      return JsonValue.NULL;
    }
    String qualifiedname = obj.getClass().getName(); // w namespace
    String clsname = obj.getClass().getSimpleName(); // w.o. namespace
    boolean sableNode = qualifiedname.startsWith("eu.exahype.node");

    if(sableNode && clsname.startsWith("A") || clsname.equals("Start")) {
      // A node can be mapped to a structure (dictionary)
      Field[] allFields = obj.getClass().getDeclaredFields();
      JsonObjectBuilder out_obj = factory.createObjectBuilder();
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
                  out_obj.add(fieldName, dump(target));
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
       return out_obj.build();
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
      String primitive = ((Token) obj).getText();
      // the most ugly way to cast java.lang.String to javax.json.JsonString:
      String key="tmp";
      return factory.createObjectBuilder().add(key, primitive).build().getJsonString(key);
    } else if (obj instanceof Collection<?>){
      // a list (such as LinkedList<PSolver>) can be mapped to an array
      //Structured lst = path.put_list(path);
      JsonArrayBuilder out_lst = factory.createArrayBuilder();
      for(Object child : (Collection)obj) {
         out_lst.add(dump(child));
      }
      return out_lst.build();
    } else if(clsname.equals("EOF")) {
      // the SableCC end node.
      return JsonValue.NULL;
    } else {
      throw new RuntimeException("Unknown class: " + clsname);
    }
  }
}
