import sys
import configparser
import re

##
# A reader for the original ExaHyPE spec file format
# 
# Sample usage:
#```
#r = SpecFile1Reader()
#r.read(path_to_spec_file) 
#
class SpecFile1Reader():
  def __init__(self):
    pass
  ##
  # Converts the original ExaHyPE spec file into an ini file.
  #
  # In the process, it performs the following changes
  # * group 'exahype-project' -> 'project'
  # * removes all 'const' tokens
  # * replaces dashes (-) in group names and keys of options by underscore (_)
  #
  # @note the SableCC parser for the original spec file did ignore whitespaces.
  #       We are thus free to place any whitespaces where it makes our work easier.
  # 
  # @return tuple containing the spec file converted to INI format (str),
  #          the number of solvers (int) and the number of plotters per solver (list of int)
  def spec_file_1_to_ini(self,spec_file_1):
    solver=0
    plotter=[]
    plotter.append(0)
    spec_file_1_ini = ""
    reads_multiline_comment = 0
    reads_singline_comment  = False
    for line in spec_file_1:
      reads_singline_comment = line.strip().startswith("//") or line.strip().startswith("#")
      reads_multiline_comment += 1 if line.strip().startswith("/*") else 0
      if reads_multiline_comment==0 and not reads_singline_comment:
        m_project = re.match(r"\s*exahype-project\s+(\w+)",line)
        m_group   = re.match(r"\s*(computational-domain|shared-memory|distributed-memory|global-optimisation)",line)
        m_solver  = re.match(r"\s*solver\s+([^\s]+)\s+(\w+)",line)
        m_plotter = re.match(r"\s*plot\s+([^\s]+)\s+(\w+)",line)
        if m_project:
          spec_file_1_ini += "[project]\n"
          spec_file_1_ini += "name="+m_project.group(1)+"\n"
        elif m_group:
          spec_file_1_ini += "[%s]" % m_group.group(1).replace("-","_").replace("global","")+"\n"
        elif m_solver:
          spec_file_1_ini += "[solver%d]" % (solver)  +"\n"
          spec_file_1_ini += "solver_type="+m_solver.group(1)+"\n"
          spec_file_1_ini += "name="+m_solver.group(2)+"\n"
        elif re.match(r"^\s*end\s*solver",line):
          solver+=1
          plotter.append(0)
        elif m_plotter:
          spec_file_1_ini += "[solver%dplotter%d]\n" % (solver,plotter[solver])
          spec_file_1_ini += "plotter_type="+m_plotter.group(1)+"\n"
          spec_file_1_ini += "name="+m_plotter.group(2)+"\n"  
          plotter[solver]+=1
        elif line.strip().startswith("end"):
          pass
        else:
          m_option = re.match(r"\s*((\w|-)+)\s*(const)?\s*=(.+)",line)
          if m_option:
            spec_file_1_ini += m_option.group(1).replace("-","_")+"="+m_option.group(4).strip()+"\n"
          else: # multiline options need to be indented
            spec_file_1_ini += "  "+line.strip()+"\n"
          
      reads_multiline_comment -= 1 if line.strip().endswith("*/") else 0 
    return (spec_file_1_ini, solver, plotter[0:-1])

  ## 
  # Convert ini file as created by method `spec_file_1_to_ini` to nested list of dict - dist of list structure.
  '
  def convert_to_dict(self,config):
    context = {}
    # root node
    root = "project"
    context[root] = {}
    for option in config.options(root):
      context[root][option.replace("-","_")] = config.get(root, option)
    # other sections
    for section in config.sections():
      if not section.startswith("solver") and not section is "project":
        context[root][section] = {}
        for option in config.options(section):
            context[root][section][option.replace("-","_")] = config.get(section, option)
    context[root]["solvers"]=[]
    for i in range(0,n_solvers):
      solver="solver%d" % i
      context[root]["solvers"].append({})
      for option in config.options(solver):
        context[root]["solvers"][i][option.replace("-","_")] = config.get(solver,option)
      context[root]["solvers"][i]["plotters"]=[]
      for j in range(0,n_plotters[i]):
        plotter="solver%dplotter%d" % (i,j)
        context[root]["solvers"][i]["plotters"].append({})
        for option in config.options(plotter):
          context[root]["solvers"][i]["plotters"][j][option] = config.get(plotter,option)
    return context
  
  ##
  # Post processes result of `convert_to_dict`, i.e. 
  # changes the structure of the dict.
  #
  # This allows it to be passed to the jinja2 template.
  # The template assumes a structure similar to the json file.
  #
  # @param context the specfile as dict
  # 
  def post_process_options(self,context):
    root = "project"
    # paths section
    context[root]["paths"]={}
    for option in context[root]:
      if option in ["log_file","peano_kernel_path","peano_toolbox_path","exahype_path","output_directory"]:
        context[root]["paths"] = context[root].pop(option)
    if "distributed_memory" in context[root]:
      context[root]["distributed_memory"]["load_balancing_type"]=context[root]["distributed_memory"].pop("identifier")
      context[root]["distributed_memory"]["buffer_size"]        =context[root]["distributed_memory"].pop("buffersize")
      # configure
      configure = context[root]["distributed_memory"].pop("configure")
      m_ranks_per_node          = re.search(r"(^|,|\s)ranks-per-node:([0-9]+)",configure)
      m_primary_ranks_per_node  = re.search(r"(^|,|\s)primary-ranks-per-node:([0-9]+)",configure)
      m_node_pool_strategy      = re.search(r"(hotspot|FCFS|sfc_diffusion)",configure)
      m_load_balancing_strategy = re.search(r"(fair|greedy_naive|greedy_regular)",configure)
      if m_ranks_per_node:
        context[root]["distributed_memory"]["ranks_per_node"]=m_ranks_per_node.group(2)
      if m_node_pool_strategy:
        context[root]["distributed_memory"]["node_pool_strategy"]=m_node_pool_strategy.group(1)
      if m_load_balancing_strategy:
        context[root]["distributed_memory"]["load_balancing_strategy"]=m_load_balancing_strategy.group(1)
    if "shared_memory" in context[root]:
      context[root]["shared_memory"]["autotuning_strategy"]=context[root]["shared_memory"].pop("identifier")
       # configure
       configure = context[root]["shared_memory"].pop("configure")
       m_background_job_consumer_threads          = re.search(r"(^|,|\s)ranks-per-node:([0-9]+)",configure)
       m_primary_ranks_per_node  = re.search(r"(^|,|\s)primary-ranks-per-node:([0-9]+)",configure)
      "background_job_consumers":{{shared_memory["background_job_consumer_threads"]}},
      "manual_pinning"                 :{{shared_memory["manual_pinning"]}},
    return context
  
  ##
  # @return dictionary representing spec file 1 content
  def read(spec_file_path):
    try:
      with open(spec_file_path, "r") as input_file:
        spec_file_1 = input_file.readlines()
    except Exception:
      raise
    
    (spec_file_1_ini, n_solvers, n_plotters) = spec_file_1_to_ini(spec_file_1)
    
    config = configparser.ConfigParser(delimiters=('='))
    config.read_string(spec_file_1_ini)
    
    return convert_to_dict(config)
