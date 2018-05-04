package minitemp.syntaxtree;

import java.lang.StringBuilder;
import minitemp.Context;

/**
 * Syntax tree node for Branch logic.
 * Store what's contained between its branches in subtree and render the correct one
 */
public class IfNode extends SyntaxTree {
  
  private Token condition;       //condition for the if
  private RootNode ifBlock;      //subtree for the if block
  private RootNode elseBlock;    //subtree for the else block
  private RootNode currentBlock; //pointer to the subtree being filled
  
  public IfNode(Token token, SyntaxTree parent) {
    condition = token;
    this.parent = parent;
    ifBlock = new RootNode();
    elseBlock = new RootNode();
    currentBlock = ifBlock;
  }
  
  /** Add a node to the current subtree */
  @Override
  public void addNode(SyntaxTree node) {
    currentBlock.addNode(node);
  }
  
  /** Trigger to start filling the else subtree, triggered when encountering the relevant else */
  public void startElseBlock() {
    currentBlock = elseBlock;
  }
  
  /** Render the subtree corresponding to the correct branch (eventually empty) */
  @Override
  public String render(Context context) throws IllegalArgumentException {
    boolean value = context.evaluateBoolean(condition.getContentClean());
    if(value) {
      return ifBlock.render(context);
    } else {
      return elseBlock.render(context);
    }
  }
  
  /** Render the subtree corresponding to the correct branch (eventually empty) */
  @Override
  public void renderWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {
    boolean value = context.evaluateBoolean(condition.getContentClean());
    if(value) {
      ifBlock.renderWithStringBuilder(context, sb);
    } else {
      elseBlock.renderWithStringBuilder(context, sb);
    }
  }
  
}
