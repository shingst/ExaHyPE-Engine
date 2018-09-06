package minitemp.syntaxtree;

import java.lang.StringBuilder;
import java.util.ArrayList;
import java.util.List;
import minitemp.Context;

/**
 * Syntax tree node with multiple children. Do not contain a token.
 * Also the root of the syntax tree or any subtree
 */
public class RootNode extends SyntaxTree {
  
  protected List<SyntaxTree> children;
  
  public RootNode() {
    this.children = new ArrayList<SyntaxTree>();
  }
  
  @Override
  public void addNode(SyntaxTree node) {
    this.children.add(node);
  }
  
  /** Performs a DFS */
  @Override
  public String render(Context context) throws IllegalArgumentException {
    return renderChildren(context);
  }
  
  /** Performs a DFS */
  @Override
  public void renderWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {
    renderChildrenWithStringBuilder(context, sb);
  }
  
  /**
   * Render the children of this node with String concatenation (inefficient for large files).
   * It is expected that the node classes will use a DFS when performing render
   *
   * @param the context of the rendering used to evaluate nodes and leafs
   * @return the rendered syntax tree
   */
  protected String renderChildren(Context context) throws IllegalArgumentException {
    String result = "";
    for(SyntaxTree child : children) {
      result += child.render(context);
    }
    return result;
  }
  
  /**
   * Render the children of this node in a Stringbuilder.
   * It is expected that the node classes will use a DFS when performing render
   *
   * @param the context of the rendering used to evaluate nodes and leafs
   * @param the StringBuilder to accumulate
   */
  protected void renderChildrenWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {
    for(SyntaxTree child : children) {
      child.renderWithStringBuilder(context, sb);
    }
  }
}
