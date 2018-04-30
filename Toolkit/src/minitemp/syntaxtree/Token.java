package minitemp.syntaxtree;

import java.util.ArrayList;
import java.lang.StringBuilder;
import minitemp.TemplateEngine;

public class Token {
  
  /** Possible types of a token */
  public static enum TokenType {
    INVALID,   /** default type, should be changed at initialization */
    TEXT,      /** type of text token */
    VAR,       /** type of variable grammar token */
    IF_OPEN,   /** type of logic grammar token: if */
    IF_ELSE,   /** type of logic grammar token: else */
    IF_CLOSE,  /** type of logic grammar token: endif */
    FOR_OPEN,  /** type of logic grammar token: for */
    FOR_CLOSE  /** type of logic grammar token: endfor */
  }
  
  public TokenType type = TokenType.INVALID; /** type of token, by default Invalid */
  private String rawContent; /** the token full string */
  
  
  public Token(String tokenAsString) {
    rawContent = tokenAsString;
    defineType();
  }
  
  
  private void defineType() {
    if(rawContent.startsWith(TemplateEngine.VAR_TOKEN_START)) {
      type = TokenType.VAR;
    } else if(rawContent.startsWith(TemplateEngine.BLOCK_TOKEN_START)) {
      String tag = rawContent;
      if(tag.startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START)) {
        tag = tag.substring(TemplateEngine.STRIP_BLOCK_TOKEN_START.length()).trim();
      } else {
        tag = tag.substring(TemplateEngine.BLOCK_TOKEN_START.length()).trim();
      }
      if(tag.startsWith(TemplateEngine.LOGIC_ENDIF_TAG)) {
        type = TokenType.IF_CLOSE;
      } else if(tag.startsWith(TemplateEngine.LOGIC_ELSE_TAG)) {
        type = TokenType.IF_ELSE;
      } else if(tag.startsWith(TemplateEngine.LOGIC_IF_TAG)) {
        type = TokenType.IF_OPEN;
      } else if(tag.startsWith(TemplateEngine.LOGIC_FOR_TAG)) {
        type = TokenType.FOR_OPEN;
      } else if(tag.startsWith(TemplateEngine.LOGIC_ENDFOR_TAG)) {
        type = TokenType.FOR_CLOSE;
      }
    } else {
      type = TokenType.TEXT;
    }
  }
  
  @Override
  public String toString() {
    return type+" | \""+rawContent+"\"";
  }
  
  /**
   * Return the clean content of a token
   * Text token: everything in the token
   * Grammar token: the inside of the delimiters - tag + trim()
   */
  public String getContentClean() {
    switch(type) {
      case TEXT:
        return rawContent;
      case VAR:
        String cleanVar = rawContent.substring(TemplateEngine.BLOCK_TOKEN_START.length());
        cleanVar = cleanVar.substring(0, cleanVar.length()-TemplateEngine.BLOCK_TOKEN_END.length());
        return cleanVar.trim();
      case IF_OPEN:
      case FOR_OPEN:
        String cleanOpen = rawContent;
        if(cleanOpen.startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START))
        {
          cleanOpen = cleanOpen.substring(TemplateEngine.STRIP_BLOCK_TOKEN_START.length());
        } else {
          cleanOpen = cleanOpen.substring(TemplateEngine.BLOCK_TOKEN_START.length());
        }
        cleanOpen = cleanOpen.substring(0, cleanOpen.length()-TemplateEngine.BLOCK_TOKEN_END.length());
        //remove the tag
        if(type == TokenType.IF_OPEN) {
          cleanOpen = cleanOpen.trim().substring(TemplateEngine.LOGIC_IF_TAG.length()).trim();
        } else if(type == TokenType.FOR_OPEN) {
          cleanOpen = cleanOpen.trim().substring(TemplateEngine.LOGIC_FOR_TAG.length()).trim();
        }
        return cleanOpen;
      default:
        return "";
    }
  }
  
  /**
   * Breaks a string into tokens using a given regex.
   * Apply strip token logic before building the token so return value != input splitted
   *
   * @param the template to evaluate as String
   * @param the regex used to split the template into tokens
   * @return the string splitted into token
   */
  static public Token[] tokenize(String template, String regex){
    String[] tokens_splitted = template.split(regex); //split the template around delimiters
    
    //assemble the tokens
    ArrayList<String> tokens_assembled = new ArrayList<String>(tokens_splitted.length/2+10); //initial guess
    StringBuilder chunk = null;
    boolean inChunk = false;
    for(int i=0; i<tokens_splitted.length; i++) {
      switch(tokens_splitted[i]){
        case TemplateEngine.VAR_TOKEN_START:
            inChunk = true;
            chunk = new StringBuilder(TemplateEngine.VAR_TOKEN_START);
            break;
        case TemplateEngine.BLOCK_TOKEN_START:
            inChunk = true;
            chunk = new StringBuilder(TemplateEngine.BLOCK_TOKEN_START);
            break;
        case TemplateEngine.VAR_TOKEN_END:
            inChunk = false;
            chunk.append(TemplateEngine.VAR_TOKEN_END);
            tokens_assembled.add(chunk.toString());
            break;
        case TemplateEngine.BLOCK_TOKEN_END:
            inChunk = false;
            chunk.append(TemplateEngine.BLOCK_TOKEN_END);
            tokens_assembled.add(chunk.toString());
            break;
        default:
          if(inChunk) {
            chunk.append(tokens_splitted[i]);
          } else {
            tokens_assembled.add(tokens_splitted[i]);
          }
      }
    }
    
    String[] tokens_raw = new String[tokens_assembled.size()];
    tokens_raw = tokens_assembled.toArray(tokens_raw);
    
    //apply strip block
    for(int i=0; i<tokens_raw.length; i++) {
      if(tokens_raw[i].startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START)) {
        if(i>0) {
          //remove trailing vertical whitespace from previous token
          tokens_raw[i-1] = tokens_raw[i-1].replaceAll("\\h*$", ""); 
        }
        if(i<tokens_raw.length-1) {
          //remove leading whitespace + newline from next token
          tokens_raw[i+1] = tokens_raw[i+1].replaceAll("^\\h*\\R?", "");
        }
      }
    }
    
    //build Token object
    Token[] tokens = new Token[tokens_raw.length];
    for(int i=0; i<tokens_raw.length; i++) {
      tokens[i] = new Token(tokens_raw[i]);
    }
    
    return tokens;
  }
  
}
