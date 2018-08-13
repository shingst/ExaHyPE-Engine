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
#
# @section DESCRIPTION
#
# Generates the application-specific Makefile
#

import sys, os, json
import shlex # > python 3.3
from string import Formatter

from .abstractModelBaseClass import AbstractModelBaseClass

escape = shlex.quote # always escape your strings, please
formatter_keys = lambda template: [i[1] for i in Formatter().parse(template)]

class KernelCallsModel(AbstractModelBaseClass):
    @staticmethod
    def obtainRuntimeProperties():
        data = {}
        data['binary'] = sys.executable

        if not data['binary']:
            # well, erhm.
            data['binary'] = "python3"

        data["env"] = {
            'PYTHONPATH': ":".join(sys.path)
        }

        # assuming that no chdir since startup
        data["working_directory"] = escape(os.getcwd())
        data["script"] = os.path.abspath(sys.argv[0])
        data["args"] = [ data['script'], "--format=any", "--validate-only"]

        return data
        #print("Stopping for debugging")
        #import ipdb; ipdb.set_trace()
    
    @classmethod
    def get_exahype_external_parser_command(cls):
        """
        Note: Considering KernelCallsImplementation.template, you'll notice that the filename of a specification file
        will be appended to this command.
        """
        runtimeProperties = cls.obtainRuntimeProperties()
        # the environment string is really "fingers crossed", without proper escaping.
        runtimeProperties["env_string"] = " ".join([ "%s='%s'" % kv for kv in runtimeProperties["env"].items() ])
        runtimeProperties["args_string"] = " ".join(runtimeProperties["args"])
        # stay readable in the Jinja. We could just do command.format(**runtimeProperties),
        # but instead we delegate the insertion until C++ sees it.

        # running it from the toolkit directory would break paths from the ExahyPE executable to it's specfile.
        # To avoid this, we could pass the specfile as stdin to the toolkit, thought. However, this is only the
        # last resort.

        # command = "cd {working_directory} && 
        command = "{env_string} {binary} {args_string}"
        command_strings = {k:v for k,v in runtimeProperties.items() if k in formatter_keys(command) }
        command_template = command.format(**{k:'" + ' + k + ' + "' for k in runtimeProperties.keys() if k in formatter_keys(command) })
        #import ipdb; ipdb.set_trace()
        return (command_template, command_strings)
    
    @staticmethod
    def specfile_as_hex(spec):
        """
        Given a native python nested dict/list object, dump it as string and then hex-encode that string
        character by character. This is safest way to include something in C++ without dealing with
        character sets or anything.
        """
        text = json.dumps(spec, sort_keys=True, indent=4)
        hex_tokens = [ "0x%02x"%ord(char) for char in text ] + ["0x00"] # null-terminated list of hex numbers
        return ", ".join(hex_tokens)

    def generateCode(self):
        return self.render("KernelCallsImplementation.template", "KernelCalls.cpp") #return path to generated file
