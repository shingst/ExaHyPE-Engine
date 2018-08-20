package minitemp.syntaxtree;

import java.lang.StringBuilder;
import minitemp.Context;

/**
 * Abstract class to build a syntax tree from the token list.
 * The template can then be rendered by the syntax tree performing a DFS
 */
public abstract class SyntaxTree {
  
  /** parent node if available */
  public SyntaxTree parent = null;
  
  /** Add a node to this node */
  public void addNode(SyntaxTree node) {}

  /** render this node naively using String and String concatenation (inefficient for large files) */
  public String render(Context context) throws IllegalArgumentException { return "";}
  
  /** render this node in a StringBuilder */
  public void renderWithStringBuilder(Context context, StringBuilder sb) throws IllegalArgumentException {}
  
  /**
   * Build a SyntaxTree from a template (String) and a regex to tokenize the template
   * Use the static Token.tokenize to tokenize the template with the regex
   *
   * @param the template to transform into a SyntaxTree
   * @param the regex to tokenize the template
   * @return a SyntaxTree ready to be rendered
   */
  static public final SyntaxTree buildTree(String template, String regex) 
      throws IllegalArgumentException
  {
    Token[] tokens = Token.tokenize(template, regex);
    
    SyntaxTree currentNode = new RootNode(); //start with a root node
    for(Token t : tokens) { //build the tree token by token
      switch(t.type) {
        case COMMENT: //comment is ignored
          break;
        case VAR: //variable is a Leaf
          currentNode.addNode(new VariableLeaf(t));
          break;
        case TEXT: //text is a Leaf
          currentNode.addNode(new TextLeaf(t));
          break;
        case IF_OPEN: //if logic start a subtree
          SyntaxTree newIfNode = new IfNode(t, currentNode);
          currentNode.addNode(newIfNode);
          currentNode = newIfNode;
          break;
        case IF_ELSE: //else stay in the if subtree but notify the node
          if (currentNode instanceof IfNode) { // the current node has to be an IfNode
            IfNode p = (IfNode)(currentNode);
            p.startElseBlock();
          } else {
            throw new IllegalArgumentException("Read an else block without being in an if block");
          }
          break;
        case IF_CLOSE: //endif, go back to the parent of the subtree
          if(!(currentNode instanceof IfNode)) {
            throw new IllegalArgumentException("Read an endif block without being in an if block");
          }
          currentNode = currentNode.parent;
          break;
        case FOR_OPEN: //for logic start a subtree
          SyntaxTree newForNode = new ForNode(t, currentNode);
          currentNode.addNode(newForNode);
          currentNode = newForNode;
          break;
        case FOR_CLOSE: //endfor, go back to the parent of the subtree
          if(!(currentNode instanceof ForNode)) {
            throw new IllegalArgumentException("Read an endfor block without being in an for block");
          }
          currentNode = currentNode.parent;
          break;
        case INVALID:
          throw new IllegalArgumentException("Invalid token read, something went wrong at initialization");
        default:
          throw new IllegalArgumentException("Unrecognized token read, something went wrong at initialization");
      }
    }
    
    return currentNode;
  }
  
}
