# Minimalist Template Engine

A minimalist template engine in java with no dependency outside the standard java libraries (java >= 1.8).

/!\ Some inputs are not sanitized, use this template engine responsibly.

Github: https://github.com/gallardjm/minitemp

## Quick start

Hello World!

    String template = "Hello {{place}}{% if excited %}!{% if mad %}!!!1!!{% endif %}{% else %}.{% endif %}";

    TemplateEngine engine = new TemplateEngine();
    Context context = new Context();
    context.put("place", "world");
    context.put("excited", true);
    context.put("mad", false);

    String result = engine.render(template, context);

## Supported template syntax

The context contains the key->values pairs used by the expressions in the grammar tokens, defined below. A context's key has to be a valid variable name (only characters from [a-zA-Z_0-9], at least one). Values can be of any java standard object type (e.g. String, Boolean, Double, Integer) or array of a standard object type. Note that lists are allowed as value and that primitive types will be casted to the corresponding object type automatically (e.g. int to Integer).

The context uses a Javascript engine (from the package javax.script) to evaluate expression in the grammar tokens using the allocated key->value pairs. All expression inside the grammar tokens have to be valid javascript expressions that can use the keys as a variable with value the one associated to the key. 

/!\ Expressions and values are not sanitized.

### Text replacement token

Default delimiters: ```{{``` and ```}}```.

Any expression resulting in a String, an int or a double is valid inside these delimiters. The token is replaced by the output of the expression. Javascript formating tools may be used.

### Logic token

Default delimiters: ```{%``` and ```%}```.

These delimiters allow some logic in the template processing. Nested logic is supported. Logic requires a start and end token.

The following logic is implemented:

* Branch ```{% if foo %} {% else %} {% endif %}```
  - the content of the if (or else) branch between the tokens will be rendered depending on the condition foo
  - foo is any boolean expression (javascript syntax)
  - else is optional

  
* Loop ```{% for foo in bar %}  {% endfor %}```
  - the content between the tokens will be rendered for each different value of foo in bar
  - bar is any Collection of object
  - if foo was already defined, its value is overridden within the loop and restored afterward
  

### Comment token

Default delimiters: ```{#``` and ```#}```.

Everything inside a comment token is ignored during rendering.

Comment delimiters work like java/C++ ```/* */``` delimiters: the first occurence of an end delimiter end the comment token (no nested comments). 

Comments cannot be used inside other tokens, as this would break the parser, e.g. ```{{ foo {# this render foo #} }}``` is incorrect, use ```{{ foo }}{# this render foo #}``` instead.

However a comment token can contain other tokens, that will be therefore ignored (comment out behavior), e.g. ```{# {{ foo }} #}``` is correct and the evaluation of foo will be ignored.

### Misc

#### Strip block

With ```{%-``` as start delimiter of a logic token, strip all whitespaces before the token + all whitespaces and one newline if present after the token. This allow more readable templates without unecessary whitespaces and newlines in the result.

Example: ``` foo \n {%- x %}   \n bar``` is evaluated as ``` foo \n{% x %} bar```.

Likewise a comment block can use this behavior by using ```{#-``` as start delimiter

