import os

import generator

import helper

class PlotterGenerator(generator.Generator):
  def __init__(self,spec,spec_file,verbose):
    generator.Generator.__init__(self,spec,spec_file,verbose)
    
    self._output_directory += "/"+spec["paths"].get("plotter_subdirectory","").strip() # TODO(Dominic): Overwrite with plotter directory for plotters
    if not os.path.exists(self._output_directory):
      print("Create plotter directory '"+self._output_directory+"'")
      os.mkdirs(self.__output_directory)
    else:
      print("Plotter directory '"+self._output_directory+"' already exists")
  
  def generatePlotter(self, context):
      user_defined = context["plotterType"]=="user::defined"
      
      template_map = {
        True : {  context["plotterName"]+".h" : "plotters/UserDefinedDeviceHeader.template",
                  context["plotterName"]+".cpp" : "plotters/UserDefinedDeviceImplementation.template" },
        False: {  context["plotterName"]+".h" : "plotters/UserOnTheFlyPostProcessingHeader.template",
                  context["plotterName"]+".cpp" : "plotters/UserOnTheFlyPostProcessingImplementation.template" }
      }
      
      for file_path in template_map.get(user_defined,{}):
        generator.Generator.render_template(self,template_map[user_defined][file_path],context,file_path,"Generated plotter file",False)
