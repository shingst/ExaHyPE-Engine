package minitemp.syntaxtree;

import java.lang.StringBuilder;
import minitemp.Context;

/**
 * Syntax tree leaf for variable '{{ var }}'
 * Currently only work with simple string replacement
 */
public class VariableLeaf extends SyntaxTree {
  
  /** Store the token for later rendering */
  private Token token;
  
  public VariableLeaf(Token token) {
    this.token = token;
  }

  @Override
  public String render(Context context) throws IllegalArgumentException {    
    return context.evaluateString(token.getContentClean());
  }

  @Override
  public void renderWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {
    sb.append(context.evaluateString(token.getContentClean()));
  }
  
}
