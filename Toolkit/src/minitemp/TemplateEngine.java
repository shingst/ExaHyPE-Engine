package minitemp;

import minitemp.syntaxtree.Token;
import minitemp.syntaxtree.SyntaxTree;

/**
 * Main class
 *
 * Once the context is build and initialized, use render(template, context) to run the template
 * engine on the template with the given context
 *
 * You can change the grammar by editing the constant from this class
 */
public class TemplateEngine {
  
  //Configuration parameters
  //------------------------
  /** Var tokens contain variable and come alone */
  public static final String VAR_TOKEN_START   = "{{";
  public static final String VAR_TOKEN_END     = "}}";
  /** Block tokens contain logic and come in group, requiring a syntax tree to evaluate */
  public static final String BLOCK_TOKEN_START = "{%";
  public static final String BLOCK_TOKEN_END   = "%}";
  /** Comment tokens, everything inside is ignored */
  public static final String COMMENT_TOKEN_START = "{#";
  public static final String COMMENT_TOKEN_END   = "#}";
  
  /** Special block and comment token start delimiter that strip whitespaces around it */
  public static final String STRIP_BLOCK_TOKEN_START = "{%-";
  public static final String STRIP_COMMENT_TOKEN_START = "{#-";
  
  /** Tags of the logic block for Branches */
  public static final String LOGIC_IF_TAG    = "if";
  public static final String LOGIC_ELSE_TAG  = "else";
  public static final String LOGIC_ENDIF_TAG = "endif";
  
  /** Tags of the logic block for Loops */
  public static final String LOGIC_FOR_TAG     = "for";
  public static final String LOGIC_FOR_SET_TAG = "in"; // {% for value in collection %}
  public static final String LOGIC_ENDFOR_TAG  = "endfor";
  
  /** Regex used to tokenize the template */
  private String regex;
  
  public TemplateEngine() {
    buildRegex(); //set this.regex
  }
  
  public String getRegex() {
    return regex;
  }
  
  /**
   * Build a regex pattern to split a string (the template) into tokens (grammar and text ones)
   *
   * A grammar token is delimiter start, content, delimiter end
   * All the text between two grammar token (or before the first/ after the last) is one text token
   *
   * The base regex find the grammar tokens delimiters, the regex then use it with the lookahead and lookbehind
   * trick to keep the grammar token matches in the split result
   *
   * The tokens should then be rebuild from the identified delimiter start, content, delimiter end
   */
  private void buildRegex() {
    // tokens
    final String varTokenStart     = escRegex(VAR_TOKEN_START);
    final String varTokenEnd       = escRegex(VAR_TOKEN_END);
    final String blockTokenStart   = escRegex(BLOCK_TOKEN_START);
    final String blockTokenEnd     = escRegex(BLOCK_TOKEN_END);
    final String commentTokenStart = escRegex(COMMENT_TOKEN_START);
    final String commentTokenEnd   = escRegex(COMMENT_TOKEN_END);
    
    //match one of the grammar tokens
    final String tokenPattern = varTokenStart+"|"+varTokenEnd+"|"+blockTokenStart+"|"+blockTokenEnd+"|"+commentTokenStart+"|"+commentTokenEnd;
    
    //lookahead+lookbehind trick
    this.regex = "(?="+tokenPattern+")|(?<="+tokenPattern+")";
  }
  
  /** escape regex special char '{' and '}' by adding '\' */
  private String escRegex(String regex) {
    // ReplaceAll uses regex itself hence the '\\' and '\\\\'
    return regex.replaceAll("\\{", "\\\\{").replaceAll("\\}", "\\\\}");
  }
  
  /**
   * Evaluate a template with a given context.
   *
   * @param the template to render as a String
   * @param the context to use to render the template (of Context class)
   * @return a string with the rendered template
   */
  public String render(String template, Context context) throws IllegalArgumentException {
    SyntaxTree st = SyntaxTree.buildTree(template, this.regex);
    StringBuilder sb = new StringBuilder(template.length()); //initial guess
    st.renderWithStringBuilder(context, sb);
    return sb.toString();
  }
  
}
