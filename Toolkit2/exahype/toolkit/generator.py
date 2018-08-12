##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
import os
import jinja2

import helper

##
# @section DESCRIPTION
#
# Abstract base class for the generators (models)
#
# initialise paths and provide a render method
class Generator():
  ##
  # Contructor method
  #
  # In child class constructor, call: `super().__init__(spec,spec_file,verbose)
  #
  # @param spec_file path to the spec_file
  # @param spec      a validated dictionary parsed from a JSON file
  def __init__(self,spec,spec_file,verbose):
    self._exahype_root         = os.path.dirname(os.path.abspath(__file__+"/../../..")) # adjust when file is moved
    self._spec_file            = str(spec_file)
    self._spec                 = spec
    self._output_directory     = self._exahype_root+"/"+spec["paths"]["output_directory"]
    self._dimensions           = spec["computational_domain"]["dimension"]
    self._verbose              = verbose
    self._jinja2_env           = jinja2.Environment(loader=jinja2.FileSystemLoader(self._exahype_root+"/Toolkit2/exahype/toolkit/templates"),trim_blocks=True)
    
    if not os.path.exists(self._output_directory):
      print("Created output directory '"+self._output_directory+"'")
      os.mkdir(self._output_directory)
    else:
      print("Output directory '"+self._output_directory+"' already exists.")
    
  ##
  # Render a template
  #
  # @param template_file_path path to a template file relative to 'Toolkit2/exahype/toolkit/templates'
  # @param context            context for rendering the spec file (content accessible via 'data' in templates)
  # @param output_file_path   path to a template file relative to the project output directory
  # @param message_success    begin of success message, e.g. "Generated file ".
  # @param overwrite          overwrite existing files
  def render_template(self,template_file_path,context,output_file_path,message_success,overwrite=False):
    if not overwrite and os.path.exists(self._output_directory+"/"+output_file_path):
      print("File '"+output_file_path+"' already exists. Is not overwritten.")
    else:
      try:
        template = self._jinja2_env.get_template(template_file_path)
        rendered_output = template.render(context)
        with open(self._output_directory+"/"+output_file_path,"w") as file_handle:
          file_handle.write(rendered_output)
        print(message_success+" '"+output_file_path+"'")
      except Exception as e:
        raise helper.TemplateNotFound("Could not generate '"+output_file_path+"' from template file '"+template_file_path+"'")
