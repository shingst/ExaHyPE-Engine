package minitemp.syntaxtree;

import java.lang.StringBuilder;
import java.util.Collection;
import minitemp.Context;
import minitemp.TemplateEngine;

/**
 * Syntax tree node for Loop logic.
 * Store what's contained between its delimiter in a subtree and repeatidly render it with
 * the loop variable changed using value from the Collection in argument
 */
public class ForNode extends SyntaxTree {
  
  private String variable;       //variable of the for
  private String collectionName; //name of the sat of values taken by variable
  private RootNode forBlock;     //subtree for the for block

  public ForNode(Token token, SyntaxTree parent) {
    String[] tokenContent = token.getContentClean().split(" "+TemplateEngine.LOGIC_FOR_SET_TAG+" ");
    variable = tokenContent[0].trim();
    collectionName = tokenContent[1].trim();
    this.parent = parent;
    forBlock = new RootNode();
  }
  
  /** Add a node to the current subtree */
  @Override
  public void addNode(SyntaxTree node) {
    forBlock.addNode(node);
  }
  
  /** Render the subtree corresponding to the correct branch (eventually empty) */
  @Override
  public String render(Context context) throws IllegalArgumentException {
    Object previousValue = context.getValueOrNull(variable); //store the previous value of variable
    Collection<Object> values = context.getCollection(collectionName);
    String result = "";
    for(Object value : values) {
      context.put(variable, value);
      result += forBlock.render(context);
    }
    context.put(variable, previousValue); //restore the previous value of variable
    
    return result;
  }
  
  /** Render the subtree corresponding to the correct branch (eventually empty) */
  @Override
  public void renderWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {
    Object previousValue = context.getValueOrNull(variable); //store the previous value of variable
    Collection<Object> values = context.getCollection(collectionName);
    for(Object value : values) {
      context.put(variable, value);
      forBlock.renderWithStringBuilder(context, sb);
    }
    context.put(variable, previousValue); //restore the previous value of variable
  }
  
}
