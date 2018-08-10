import sys
import configparser
import re
import json

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
          spec_file_1_ini += "solver_type="+m_solver.group(1)+"\n"
          spec_file_1_ini += "name="+m_solver.group(2)+"\n"
        elif re.match(r"^\s*end\s*solver",line):
          solver+=1
          plotter.append(0)
        elif m_plotter:
          spec_file_1_ini += "[solver%dplotter%d]\n" % (solver,plotter[solver])
          spec_file_1_ini += "type="+m_plotter.group(1)+"\n"
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
  def convert_to_dict(self,config,n_solvers,n_plotters):
    context = {}
    # root node
    root = "project"
    context = {}
    for option in config.options(root):
      context[option.replace("-","_")] = config.get(root, option)
    # other sections
    for section in config.sections():
      if not section.startswith("solver") and not section=="project":
        context[section] = {}
        for option in config.options(section):
            context[section][option.replace("-","_")] = config.get(section, option)
    context["solvers"]=[]
    for i in range(0,n_solvers):
      solver="solver%d" % i
      context["solvers"].append({})
      for option in config.options(solver):
        context["solvers"][i][option.replace("-","_")] = config.get(solver,option)
      context["solvers"][i]["plotters"]=[]
      for j in range(0,n_plotters[i]):
        plotter="solver%dplotter%d" % (i,j)
        context["solvers"][i]["plotters"].append({})
        for option in config.options(plotter):
          context["solvers"][i]["plotters"][j][option] = config.get(plotter,option)
    return context
  
  ##
  # Convets the kernel terms of a solver entry in the original spec file
  #
  # @return tuple consisting of parsed context and the number of point sources
  def map_kernel_terms(kernel_terms):
    n_point_sources
    context = {}
    for token in kernel_terms.split(","):
    token_s = token.strip() 
    for term in ["flux","source","ncp","pointsources","materialparameters"]:
      if token_s.startswith(terms):
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
  def map_aderdg_kernel_opts(kernel_opts):
    context = {}
    stp     = "space_time_predictor"
    opt     = "optimised_terms"
    opt_dbg = "optimised_kernel_debugging"
    context[stp]={}
    context[opt]=[]
    context[opt_dbg]=[]
    for token in kernel_terms.split(","):
      token_s = token.strip() 
      for term in ["generic","optimised","user"]:
        if token_s==term:
          context["implementation"]=term
      if token_s=="patchwiseadjust":
        context["adjust_solution"]="patchwise"
      if token_s=="usestack"
        context["allocate_temporary_arrays"]="stack"
      for term in ["fusedsource","fluxvect","fusedsourcevect"]
        if token_s==term:
          context[opt].append(term)
    for term in ["converter","flops"]:
        if token_s==term:
          context[opt_debug].append(term)
    for term in ["cerkguess","notimeavg","maxpicarditer"]
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
  def map_fv_kernel_opts(kernel_opts):
    context = {}
    for token in kernel_terms.split(","):
      token_s = token.strip() 
      for term in ["generic","optimised","user"]:
        if token_s==term:
          context["implementation"]=term
      if token_s=="patchwiseadjust":
        context["adjust_solution"]="patchwise"
      if token_s=="usestack"
        context["allocate_temporary_arrays"]="stack"
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
    # paths 
    context["paths"]={}
    for option in context:
      if option in ["log_file","peano_kernel_path","peano_toolbox_path","exahype_path","output_directory"]:
        context["paths"] = context.pop(option)
    # optimisation
    if "optimisation" in context:
      for option in context["optimisation"]:
        if context["optimisation"][option] in ["on","off"]:
          context["optimisation"][option]="false" if context["optimisation"][option]=="off" else "true"
    # distributed_memory
    if "distributed_memory" in context:
      context["distributed_memory"]["load_balancing_type"]=context["distributed_memory"].pop("identifier")
      context["distributed_memory"]["buffer_size"]        =context["distributed_memory"].pop("buffersize")
      # configure
      configure = context["distributed_memory"].pop("configure")
      m_ranks_per_node          = re.search(r"(^|,|\s)ranks-per-node:([0-9]+)",configure)
      m_primary_ranks_per_node  = re.search(r"(^|,|\s)primary-ranks-per-node:([0-9]+)",configure)
      m_node_pool_strategy      = re.search(r"(hotspot|FCFS|sfc_diffusion)",configure)
      m_load_balancing_strategy = re.search(r"(fair|greedy_naive|greedy_regular)",configure)
      if m_ranks_per_node:
        context["distributed_memory"]["ranks_per_node"]=m_ranks_per_node.group(2)
      if m_node_pool_strategy:
        context["distributed_memory"]["node_pool_strategy"]=m_node_pool_strategy.group(1)
      if m_load_balancing_strategy:
        context["distributed_memory"]["load_balancing_strategy"]=m_load_balancing_strategy.group(1)
    # shared_memory
    if "shared_memory" in context:
       context["shared_memory"]["autotuning_strategy"]=context["shared_memory"].pop("identifier")
       # configure
       configure = context["shared_memory"].pop("configure")
       m_background_job_consumers = re.search(r"background_task:([0-9]+)",configure)
       if m_background_job_consumers:
         context["shared_memory"]["background_job_consumers"] = m.background_job_consumers.group(1)
       if re.search(r"manual_pinning",configure)!=None:
         context["shared_memory"]["manual_pinning"] = True
     for i,solver in enumerate(context["solvers"]):
        if solver["solver_type"]=="Limiting-ADER-DG":
          context["solvers"][i]["limiter"]={}
          for option in [ "dmp_observables", "dmp_relaxation_parameter", "dmp_difference_scaling", "helper_layers", "steps_till_cured" ]
            if option in solver:
              context["solvers"][i]["limiter"][options] = context["solvers"][i].pop(options)
          # kernels
          n_point_sources = 0
          aderdg_kernel_type, aderdg_kernel_terms, aderdg_kernel_opts  = ""
          fv_kernel_type, fv_kernel_terms, fv_kernel_opts  = ""
          # aderdg 
          if solver["solver_type"] in ["Limiting-ADER-DG","ADER-DG"]:
              context["solvers"][i]["aderdg_kernel"]={}
              if "language" in context["solvers"][i]:
                context["solvers"][i]["aderdg_kernel"]["language"]=context["solvers"][i].pop("language")
              if "type" in context["solvers"][i]:
                aderdg_kernel_type  = context["solvers"][i].pop("type") 
              if "terms" in context["solvers"][i]:
                aderdg_kernel_terms = context["solvers"][i].pop("terms")
              if "optimisation" in context["solvers"][i]:
                aderdg_kernel_opts  = context["solvers"][i].pop("optimisation")
              # type
              for token in aderdg_kernel_type.split(","):
                token_s = token.strip() 
                if token_s in ["linear","nonlinear"]:
                  context["solvers"][i]["aderdg_kernel"]["nonlinear"]=token_s=="nonlinear"
                if token_s in ["Legendre","Lobatto"]:
                  context["solvers"][i]["aderdg_kernel"]["basis"]=token_s
              # terms
              result, n_point_sources = map_kernel_terms(aderdg_kernel_terms)
              context["solvers"][i]["aderdg_kernel"].update(result)
              # opts
              context["solvers"][i]["aderdg_kernel"].update(map_aderdg_kernel_opts(aderdg_kernel_opts))
              
          # fv
          if solver["solver_type"]=="Finite-Volumes":
              context["solvers"][i]["fv_kernel"]={}
              if "language" in context["solvers"][i]:
                context["solvers"][i]["fv_kernel"]["language"]=context["solvers"][i].pop("language")
              if "type" in context["solvers"][i]:
                fv_kernel_type  = context["solvers"][i].pop("type")
              if "terms" in context["solvers"][i]:
               fv_kernel_terms = context["solvers"][i].pop("terms")
              if "optimisation" in context["solvers"][i]:
                fv_kernel_opts  = context["solvers"][i].pop("optimisation")
          if solver["solver_type"]=="Limiting-ADER-DG":
              context["solvers"][i]["fv_kernel"]={}
              if "language" in context["solvers"][i]:
                context["solvers"][i]["fv_kernel"]["language"]=context["solvers"][i].pop("limiter-language")
              if "type" in context["solvers"][i]:
                fv_kernel_type  = context["solvers"][i].pop("limiter-type") 
              if "type" in context["solvers"][i]:
                fv_kernel_terms = context["solvers"][i].pop("limiter-terms")
              if "type" in context["solvers"][i]:
                fv_kernel_opts  = context["solvers"][i].pop("limiter-optimisation")
          # fv type
          for token in fv_kernel_type.split(","):
            token_s = token.strip() 
            if token_s in ["godunov","musclhancock"]:
              context["solvers"][i]["fv_kernel"]["scheme"]=token_s
          # fv terms
          result, n_point_sources = map_kernel_terms(fv_kernel_terms)
          context["solvers"][i]["fv_kernel"].update(result)
          # fv opts
          context["solvers"][i]["fv_kernel"].update(map_aderdg_kernel_opts(aderdg_kernel_opts))
    }
    
     
         
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
    
    return post_process_options(convert_to_dict(config,n_solvers,n_plotters))
