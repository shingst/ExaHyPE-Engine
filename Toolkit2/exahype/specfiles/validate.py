# Python3

"""
Validator for JSON specfiles ("specfiles 2.0"), using the a JSON-Schema library

Test it like this:

> testfile = json.load(open("../../examples/EulerFlow-example.exahype2", "r"))
> enriched_with_defaults = validate(testfile)

It will raise an Exception if the file is not valid.
"""

# python-included ("batteries")
import sys, json, pathlib

import specfile1_reader

# https://pypi.org/project/jsonschema/
# Current stable version is 2.6, version 3.0 brings Draft6Validator,
# but that's still alpha.

sys.path.append("../../")     # in case of local installation

from jsonschema import Draft4Validator, validators, validate

exahype_schema_filename = "../../exahype-specfile.schema.json"
schema = json.load(pathlib.Path(__file__).parent.joinpath(exahype_schema_filename).open("r"))

def extend_with_default(validator_class):
  """
  JSON-Schema Validator which sets the default values.
  This is a bit experimental but in genreal it is good practise to have config
  files which only describe the deviation from reasonable default values.
  
  Code comes from https://python-jsonschema.readthedocs.io/en/latest/faq/#why-doesn-t-my-schema-s-default-property-set-the-default-on-my-instance
  """  
  validate_properties = validator_class.VALIDATORS["properties"]

  def set_defaults(validator, properties, instance, schema):
    for property, subschema in properties.items():
      if "default" in subschema:
        instance.setdefault(property, subschema["default"])

    for error in validate_properties(
      validator, properties, instance, schema,
    ):
      yield error

  return validators.extend(
    validator_class, {"properties" : set_defaults},
  )

SimpleValidator = Draft4Validator # should try to use Draft6Validator if available...
ExtendingValidator = extend_with_default(Draft4Validator)

def get_validator(set_defaults=True):
  validator = ExtendingValidator if set_defaults else SimpleValidator
  return validator(schema)

def validate(json_filename_or_file, set_defaults=True):
  if not hasattr(json_filename_or_file, "read"): # is a file name then
    json_filename_or_file = open(json_filename_or_file, "r")
  input_structure = json.load(json_filename_or_file)
  get_validator(set_defaults).validate(input_structure)
  return input_structure
  
def validate_specfile1(specfile1_filename):
  input_structure = specfile1_reader.SpecFile1Reader().read(open(specfile1_filename, "r"))
  print("Reading legacy ExaHyPE specification file format... OK")
  print("Translating file to JSON format ... OK")
  print("Result:")
  print(json.dumps(input_structure,indent=2))
  get_validator(set_defaults=True).validate(input_structure)
  return input_structure
  
