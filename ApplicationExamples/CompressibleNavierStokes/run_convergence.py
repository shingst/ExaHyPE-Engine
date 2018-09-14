#!/usr/bin/env python3
import argparse
import hashlib
import logging
import os
import subprocess
import json
import copy
import shutil

template = json.loads(r'''
{
  "project_name": "NavierStokes",
  "paths": {
    "peano_kernel_path": "./Peano/",
    "exahype_path": "./ExaHyPE",
    "output_directory": "./ApplicationExamples/CompressibleNavierStokes"
  },
  "computational_domain": {
    "dimension": 2,
    "end_time": 0.6,
    "offset": [
      0.0,
      0.0
    ],
    "width": [
      10.0,
      10.0
    ]
  },
  "shared_memory": {
    "cores": 16,
    "properties_file": "sharedmemory.properties",
    "autotuning_strategy": "dummy"
  },
  "optimisation": {
    "fuse_algorithmic_steps": true,
    "fuse_algorithmic_steps_factor": 0.99,
    "spawn_predictor_as_background_thread": true,
    "spawn_amr_background_threads": true,
    "disable_vertex_exchange_in_time_steps": true,
    "time_step_batch_factor": 1.0,
    "disable_metadata_exchange_in_batched_time_steps": true,
    "double_compression": 0.0,
    "spawn_double_compression_as_background_thread": true
  },
  "solvers": [
    {
      "type": "ADER-DG",
      "name": "NavierStokesSolverDG",
      "order": 1,
      "maximum_mesh_size": 0.37,
      "maximum_mesh_depth": 10,
      "time_stepping": "global",
      "aderdg_kernel": {
        "language": "C",
        "nonlinear": true,
        "terms": [
          "flux",
          "viscous_flux",
          "source"
        ],
        "space_time_predictor": {},
        "optimised_terms": [],
        "optimised_kernel_debugging": [],
        "implementation": "generic"
      },
      "point_sources": 0,
      "variables": [
        {
          "name": "rho",
          "multiplicity": 1
        },
        {
          "name": "j",
          "multiplicity": 2
        },
        {
          "name": "E",
          "multiplicity": 1
        }
      ],
      "parameters": {
        "viscosity": 0.1,
        "scenario": "convergence"
      },
      "plotters": [
        {
          "type": "vtu::Legendre::vertices::ascii",
          "name": "Plotter",
          "time": 0.5,
          "repeat": 0.5,
          "output": "./results/solution",
          "variables": 5
        }
      ]
    }
  ]
}

''')

def get_meshsize(factor):
    # We need to make sure we are slightly larger than actual mesh-size
    # otherwise wrong mesh size might be chosen.
    eps = 10e-3
    return 10/3**factor + eps

def get_exahype_root():
    return os.path.realpath(os.path.dirname(__file__) + "/../..") 

def get_application_path():
    # TODO
    return os.path.realpath(os.path.dirname(__file__))

def render_template(template, config, file_name):
    template['solvers'][0]['order'] = config['order']
    template['solvers'][0]['maximum_mesh_size'] = config['mesh_size']
    template['solvers'][0]['plotters'][0]['output'] = \
        bin_dir = get_application_path() + '/convergence/results/solution_order_{}_{}'.format(
            config['order'], config['factor'])
    rendered_template = json.dumps(template, indent=4)

    # Quick hack to ensure semi-unique name.
    template_hash = hashlib.sha512(rendered_template.encode('utf-8')).hexdigest()[:8]
    file_name = file_name.format(template_hash=template_hash)
    logging.info("Created file {file_name}".format(file_name=file_name))

    with open(file_name, 'w') as f:
        f.write(rendered_template)

    return file_name

def run(template, my_env, config):
    exahype_bin = get_application_path() +\
        "/convergence/bin/exahype_order_{}".format(config['order'])
    
    file_name = 'convergence/rendered_template_{template_hash}.exahype2'
    template_path = render_template(template, config, file_name)

    logging.info("Starting with {config}".format(config=config))
    subprocess.run([exahype_bin, template_path], env=my_env)
    logging.info("Ran file {file_name}.".format(file_name=template_path))
    

def build(template, my_env, config):
    # TODO: Make sure all directories exist.
    
    os.chdir(get_application_path()) # move to correct directory

    toolkit_bin = get_exahype_root() + '/Toolkit/toolkit.sh'

    file_name = 'convergence/rendered_template_{template_hash}.exahype2'
    template_path = render_template(template, config, file_name)
    
    logging.info("Started building with {}".format(config))

    # 1. Clean
    subprocess.run(['make', 'clean' ,'-j{}'.format(os.cpu_count())], env=my_env)
    logging.info("Cleaned.")

    # 2. Run toolkit
    subprocess.run([toolkit_bin, template_path], env=my_env)
    logging.info("Ran toolkit.")
    
    # 3. Compile
    subprocess.run(['make', '-j{}'.format(os.cpu_count())], env=my_env)
    logging.info("Compiled.")

    # 4. Move binary
    compiled_name = 'ExaHyPE-{}'.format(template['project_name'])
    new_name =  'exahype_order_{}'.format(config['order'])
    bin_dir = get_application_path() + "/convergence/bin/"
    new_path = os.path.realpath(bin_dir + new_name)

    shutil.copy2(compiled_name, new_path)
    logging.info("Copied file to {}".format(new_path))

    return new_path
    

def main():
    parser = argparse.ArgumentParser(description='Run simulations for various orders and mesh-sizes.')
    parser.add_argument('--build', '-b', action='store_true')
    parser.add_argument('--run', '-r', action='store_true')

    args = parser.parse_args()
    
    print(__file__)
    print(get_exahype_root())

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    my_env = os.environ.copy()
    my_env['COMPILER'] = 'GNU'
    my_env['TBB_INC'] = '/usr/include/tbb'
    my_env['TBB_SHLIB'] = '-L/lib64 -ltbb -lpthread'
    my_env['SHAREDMEM'] = 'TBB'
    my_env['MODE'] = 'Release'
    
    max_factor = 3

    order_grid = [1,2,3,4]

    # Build for various orders.
    if args.build:
        logging.info("Start compiling.")
        for order in order_grid:
            config = {'order': order,
                    'factor': max_factor,
                    'mesh_size': get_meshsize(max_factor)}
            build(template=template, my_env=my_env, config=config) 
    
    if args.run:
        logging.info("Start running.")
        for order in order_grid:
            for factor in range(1, max_factor + 1):
                config = {'order': order,
                        'factor': factor,
                        'mesh_size': get_meshsize(factor)}
                run(template=template, my_env=my_env, config=config)

if __name__ == '__main__':
    main()
