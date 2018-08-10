import sys
import configparser
import re
import json
import collections

##
# A reader for the original ExaHyPE specification file format
# 
# Sample usage:
#```
# r = SpecFile1Reader()
# my_json_obj = r.read()
# my_json_file_content = json.dumps(my_json_ob)
#```
class SpecFile1Reader():
  ##
  # Converts the original ExaHyPE spec file into an ini file.
  #
  # In the process, it performs the following changes
  # * group 'exahype-project'      -> 'project'
  # * group 'global-optimisations' -> 'optimisations'
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
          spec_file_1_ini += "project_name="+m_project.group(1)+"\n"
        elif m_group:
          spec_file_1_ini += "[%s]" % m_group.group(1).replace("-","_").replace("global_","")+"\n"
        elif m_solver:
          spec_file_1_ini += "[solver%d]" % (solver)  +"\n"
          spec_file_1_ini += "solver_type="+m_solver.group(1)+"\n" # will be replaced later on
          spec_file_1_ini += "name="+m_solver.group(2)+"\n"
        elif re.match(r"^\s*end\s*solver",line):
          solver+=1
          plotter.append(0)
        elif m_plotter:
          spec_file_1_ini += "[solver%dplotter%d]\n" % (solver,plotter[solver])
          spec_file_1_ini += "type="+m_plotter.group(1)+"\n"
          spec_file_1_ini += "name="+m_plotter.group(2)+"\n"  
          plotter[solver]+=1
        elif line.strip().startswith("end "):
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
  # Tries to convert certain options to the expected type.
  # Returns string if it fails
  def convert(self,option,value):
    integers=[\
      "dimension",\
      "time_steps",\
      "buffer_size",\
      "timeout",\
      "cores",\
      "order",\
      "patch_size",\
      "maximum_mesh_depth",\
      "dmp_observables",\
      "steps_till_cured",\
      "helper_layers"\
    ]
    numbers=[\
      "end_time",\
      "maximum_mesh_size",\
      "time",\
      "repeat",\
      "fuse_algorithmic_steps_factor",\
      "time_step_batch_factor",\
      "double_compression",\
      "dmp_relaxation_parameter",\
      "dmp_difference_scaling"\
    ]
    try:
      if option in integers:
        return int(value)
      elif option in numbers:
        return float(value)
      else:
        return value 
    except:
      return value # let the JSON schema validation deal with the error handling

  ## 
  # Convert ini file as created by method `spec_file_1_to_ini` to nested list of dict - dist of list structure.
  def convert_to_dict(self,config,n_solvers,n_plotters):
    # root node
    root = "project"
    context = collections.OrderedDict()
    for option in config.options(root):
      context[option.replace("-","_")] = config.get(root, option)
    # other sections
    for section in config.sections():
      if not section.startswith("solver") and not section==root:
        context[section] = collections.OrderedDict()
        for option in config.options(section):
            context[section][option.replace("-","_")] = self.convert(option, config.get(section, option))
    # solvers
    context["solvers"]=[]
    for i in range(0,n_solvers):
      solver="solver%d" % i
      context["solvers"].append(collections.OrderedDict())
      for option in config.options(solver):
        context["solvers"][i][option.replace("-","_")] = self.convert(option, config.get(solver,option))
      # plotters
      context["solvers"][i]["plotters"]=[]
      for j in range(0,n_plotters[i]):
        plotter="solver%dplotter%d" % (i,j)
        context["solvers"][i]["plotters"].append(collections.OrderedDict())
        for option in config.options(plotter):
          context["solvers"][i]["plotters"][j][option] = self.convert(option, config.get(plotter,option))
    return context
  
  ##
  # TODO
  def map_computational_domain(self,domain):
    for option in ["offset","width"]:
      token = domain.pop(option)
      domain[option]=[]
      for value in token.split(","):
        domain[option].append(float(value))
  ##
  # TODO
  def map_distributed_memory(self,distributed_memory):
    distributed_memory["load_balancing_type"]=distributed_memory.pop("identifier").replace("_load_balancing","")
    distributed_memory["buffer_size"]        =distributed_memory.pop("buffer_size")
    # configure
    configure = distributed_memory.pop("configure")
    m_ranks_per_node          = re.search(r"(^|,|\s)ranks_per_node:([0-9]+)",configure)
    m_primary_ranks_per_node  = re.search(r"(^|,|\s)primary_ranks_per_node:([0-9]+)",configure)
    m_node_pool_strategy      = re.search(r"(hotspot|FCFS|sfc_diffusion)",configure)
    m_load_balancing_strategy = re.search(r"(fair|greedy_naive|greedy_regular)",configure)
    if m_ranks_per_node:
      distributed_memory["ranks_per_node"]        =int(m_ranks_per_node.group(2))
    if m_primary_ranks_per_node:
      distributed_memory["primary_ranks_per_node"]=int(m_ranks_per_node.group(2))
    if m_node_pool_strategy:
      distributed_memory["node_pool_strategy"]=m_node_pool_strategy.group(1)
    if m_load_balancing_strategy:
      distributed_memory["load_balancing_strategy"]=m_load_balancing_strategy.group(1)
  
  ##
  # TODO
  def map_shared_memory(self,shared_memory):
    shared_memory["autotuning_strategy"]=shared_memory.pop("identifier")
    # configure
    configure = shared_memory.pop("configure")
    m_background_job_consumers = re.search(r"background_task:([0-9]+)",configure)
    if m_background_job_consumers:
     shared_memory["background_job_consumers"] = m.background_job_consumers.group(1)
    if re.search(r"manual_pinning",configure)!=None:
     shared_memory["manual_pinning"] = True
    
  ##
  # Converts the kernel terms of a solver entry in the original spec file
  #
  # @return tuple consisting of parsed context and the number of point sources
  def map_kernel_terms(self,kernel_terms):
    n_point_sources = 0
    context = []
    for token in kernel_terms.split(","):
      token_s = token.strip() 
      for term in ["flux","source","ncp","pointsources","materialparameters"]:
        if token_s.startswith(term):
          context.append(term)
          if term=="pointsources":
            try:
              n_point_sources = int(token_s.split(":")[-1])
              print("WARNING: Found 'pointsources' term. Ensure that you specify field 'point_sources' in the generated JSON file.",file=sys.stderr)
            except:
              print("ERROR: Number of point sources could not be parsed in original ExaHyPE specification file (is: '%s'. expected: 'pointsources':<int>)!" % token_s,file=sys.stderr)
              sys.exit()
      return (context,n_point_sources)
  
  ##
  # Converts the "optimisation" string of a ADER-DG and Limiting-ADERDG solver found 
  # in the original spec file
  def map_aderdg_kernel_opts(self,kernel_opts):
    context = collections.OrderedDict()
    stp     = "space_time_predictor"
    opt     = "optimised_terms"
    opt_dbg = "optimised_kernel_debugging"
    context[stp]=collections.OrderedDict()
    context[opt]=[]
    context[opt_dbg]=[]
    for token in kernel_opts.split(","):
      token_s = token.strip() 
      for term in ["generic","optimised","user"]:
        if token_s==term:
          context["implementation"]=term
      if token_s=="patchwiseadjust":
        context["adjust_solution"]="patchwise"
      if token_s=="usestack":
        context["allocate_temporary_arrays"]="stack"
      for term in ["fusedsource","fluxvect","fusedsourcevect"]:
        if token_s==term:
          context[opt].append(term)
    for term in ["converter","flops"]:
        if token_s==term:
          context[opt_debug].append(term)
    for term in ["cerkguess","notimeavg","maxpicarditer"]:
      if token_s.startswith(term):
        if term=="maxpicarditer":
          try:
            context[stp][term]=int(token_s.split(":")[-1])
          except:
            print("ERROR: Number of point sources could not be parsed in original ExaHyPE specification file (is: '%s'. expected: 'pointsources':<int>)!" % token_s,file=sys.stderr)
            sys.exit()
        else:
          context[stp][term]=term
    return context
    
  ##
  # Converts the "optimisation" string of a Finite-Volumes and the "limiter-optimisation" string of Limiting-ADERDG solver 
  # found in the original spec file
  def map_fv_kernel_opts(self,kernel_opts):
    context = collections.OrderedDict()
    for token in kernel_opts.split(","):
      token_s = token.strip()
      for term in ["generic","optimised","user"]:
        if token_s==term:
          context["implementation"]=term
      if token_s=="patchwiseadjust":
        context["adjust_solution"]="patchwise"
      if token_s=="usestack":
        context["allocate_temporary_arrays"]="stack"
    return context
  
  ##
  # TODO
  def map_variables(self,variables):
    if re.match("\s*([0-9]+)\s*",variables):
      return int(variables.strip())
    else:
      result=[]
      it = re.finditer("(\w+)\s*:\s*([0-9]+)",variables)
      for m in it:
        result.append(collections.OrderedDict([ ( "name",m.group(1) ), ( "multiplicity",int(m.group(2)) ) ])) 
      if result:
        return result
      else:
        return [variables]
  
  ##
  # TODO
  def map_constants(self,constants):
    result=[]
    it = re.finditer("(\s|,|^)([^\s,]+)\s*:\s*([^\s,]+)",constants)
    for m in it:
      result.append({ m.group(2) : m.group(3) }) 
    if result:
      return result
    else:
      return [constants]
    
  ##
  # Post processes result of `convert_to_dict`, i.e. 
  # changes the structure of the dict.
  #
  # This allows it to be passed to the jinja2 template.
  # The template assumes a structure similar to the json file.
  #
  # @param context the specfile as dict
  # 
  def map_options(self,context):
    # paths, optimisation, distributed_memory, shared_memory
    context["paths"]=collections.OrderedDict()
    context.move_to_end("paths",last=False)        # put on top
    context.move_to_end("project_name",last=False) # put on top
    for option in list(context.keys()):
      if option in ["log_file","peano_kernel_path","peano_toolbox_path","exahype_path","output_directory"]:
        context["paths"][option] = context.pop(option)
    self.map_computational_domain(context["computational_domain"])
    if "optimisation" in context:
      for option in context["optimisation"]:
        if context["optimisation"][option] in ["on","off"]:
          context["optimisation"][option]=False if context["optimisation"][option]=="off" else True
    if "distributed_memory" in context:
      self.map_distributed_memory(context["distributed_memory"])  
    if "shared_memory" in context:
      self.map_shared_memory(context["shared_memory"])
    # solvers
    for i,solver in enumerate(context["solvers"]):
      # type
      old_type = solver["type"]
      solver["type"]=solver.pop("solver_type")
      solver.move_to_end("type",last=False)  # put on top
      
      # kernels
      n_point_sources = 0
      aderdg_kernel_type, aderdg_kernel_terms, aderdg_kernel_opts  = "", "", ""
      fv_kernel_type, fv_kernel_terms, fv_kernel_opts  = "", "", ""
      # aderdg 
      if solver["type"] in ["Limiting-ADER-DG","ADER-DG"]:
          solver["aderdg_kernel"]=collections.OrderedDict()
          if "language" in solver:
            solver["aderdg_kernel"]["language"]=solver.pop("language")
          if "type" in solver:
            aderdg_kernel_type  = old_type
          if "terms" in solver:
            aderdg_kernel_terms = solver.pop("terms")
          if "optimisation" in solver:
            aderdg_kernel_opts  = solver.pop("optimisation")
          # type
          for token in aderdg_kernel_type.split(","):
            token_s = token.strip() 
            if token_s in ["linear","nonlinear"]:
              solver["aderdg_kernel"]["nonlinear"]=token_s=="nonlinear"
            if token_s in ["Legendre","Lobatto"]:
              solver["aderdg_kernel"]["basis"]=token_s
          # terms
          result, n_point_sources = self.map_kernel_terms(aderdg_kernel_terms)
          solver["aderdg_kernel"]["terms"]=result
          # opts
          solver["aderdg_kernel"].update(self.map_aderdg_kernel_opts(aderdg_kernel_opts))
      # limiter
      if solver["type"]=="Limiting-ADER-DG":
        solver["limiter"]=collections.OrderedDict()
        for option in [ "dmp_observables", "dmp_relaxation_parameter", "dmp_difference_scaling", "helper_layers", "steps_till_cured" ]:
          if option in solver:
            solver["limiter"][option] = solver.pop(option)
        if "implementation" in solver["aderdg_kernel"]:
          solver["limiter"]["implementation"]=solver["aderdg_kernel"]["implementation"]
      # fv
      if solver["type"]=="Finite-Volumes":
          solver["fv_kernel"]=collections.OrderedDict()
          if "language" in solver:
            solver["fv_kernel"]["language"]=solver.pop("language")
          if "type" in solver:
            fv_kernel_type  = old_type
          if "terms" in solver:
           fv_kernel_terms = solver.pop("terms")
           # fv terms
          result, n_point_sources = self.map_kernel_terms(fv_kernel_terms)
          solver["fv_kernel"]["terms"]=result
          if "optimisation" in solver:
            fv_kernel_opts  = solver.pop("optimisation")
      if solver["type"]=="Limiting-ADER-DG":
          solver["fv_kernel"]=collections.OrderedDict()
          if "limiter_language" in solver:
            solver["fv_kernel"]["language"]=solver.pop("limiter_language")
          if "limiter_type" in solver:
            fv_kernel_type  = solver.pop("limiter_type") 
          if "terms" in solver["aderdg_kernel"]:
            solver["fv_kernel"]["terms"]=solver["aderdg_kernel"]["terms"] # copy ADER-DG terms
          if "limiter_optimisation" in solver:
            fv_kernel_opts  = solver.pop("limiter_optimisation")
      # fv type
      for token in fv_kernel_type.split(","):
        token_s = token.strip() 
        if token_s in ["godunov","musclhancock"]:
          solver["fv_kernel"]["scheme"]=token_s
      # fv opts
      if "fv_kernel" in solver: 
        solver["fv_kernel"].update(self.map_fv_kernel_opts(fv_kernel_opts))
      
      # variables, parameters, and more
      solver["variables"]=self.map_variables(solver.pop("variables"))
      if "parameters" in solver:
        solver["material_parameters"]=self.map_variables(solver.pop("parameters"))
      if "global_observables" in solver:
        solver["global_observables"]=self.map_variables(solver.pop("global_observables"))
      if "constants" in solver:
        solver["parameters"]=self.map_constants(solver.pop("constants"))
      
      # plotters
      solver["plotters"] = solver.pop("plotters") # put at end
      for j,plotter in enumerate(solver["plotters"]):
        plotter["variables"]=self.map_variables(plotter.pop("variables"))
        if "select" in plotter:
          plotter["select"]=self.map_constants(plotter.pop("select"))
    
    return context
  
  ##
  # @return dictionary storing spec file 1 content expressed
  # in terms of a spec file 2 JSON object
  #
  # Sample usage:
  #```
  # r = SpecFile1Reader()
  # my_json_obj = r.read()
  # my_json_file_content = json.dumps(my_json_ob)
  # ```
  def read(self,spec_file_path):
    try:
      with open(spec_file_path, "r") as input_file:
        spec_file_1 = input_file.readlines()
    except Exception:
      raise
    (spec_file_1_ini, n_solvers, n_plotters) = self.spec_file_1_to_ini(spec_file_1)
    
    config = configparser.ConfigParser(delimiters=('='))
    config.read_string(spec_file_1_ini)
    
    return self.map_options(self.convert_to_dict(config,n_solvers,n_plotters))
