#include <stdio.h>
#include <iostream>

#include "mexa.h"

#include "exahype/Parser.h"

#include <utility> // pair
#include <regex> 

int main(int argc, char** argv) {
	std::string progname = argv[0];
	if (argc < 2) {
		printf("Usage: %s <someSpecfile.exahype>\n", progname.c_str());
		return -1;
	}

	// cmdlineargs contains all argv expect the progname.
	std::vector<std::string> cmdlineargs(argv + 1, argv + argc);
	std::string firstarg = cmdlineargs[0];

	printf("Reading specfile %s...\n", firstarg.c_str());
	
	// run the program with EXAHYPE_VERBOSE_PARSER="TRUE"
	// to see how the parser parses the tokenStream.
	exahype::Parser parser;
	parser.readFile(firstarg);

	if (!parser.isValid()) {
		printf("Invalid specfile\n");
		return -2;
	}
	
	exahype::Parser::ParserView view(parser.getParserView(0));
	mexa::mexafile mf = mexa::fromSpecfile(view.getAllAsOrderedMap(), view.toString());
	
	std::cout << "Read in mexafile:\n" << mf.toString();
	
	std::cout << "\nTesting string access:";
	std::cout << "\ntest1 = " << mf("string/test1").get_string();
	std::cout << "\ntest2 = " << mf("string").get_string("test2");
	

	std::cout << "\nTesting queries in general:";
	std::cout << "\nQueried: " << mf.query("boundaries").toString();
	std::cout << "\nFiltered: " << mf.filter("boundaries").toString();
	std::cout << "\nSugar: " << mf("boundaries").toString();
	
	mexa::value v;
	std::cout << "\nTesting casting:";
	v = mf("casting/float").get_value();
	std::cout << "\nFrom float: "
		<< "double: " << v.as_double() << ", "
		<< "int: " << v.as_int() << ", "
		<< "bool: " << v.as_bool() << ", "
		<< "string: " << v.as_string();
	v = mf("casting/int").get_value();
	std::cout << "\nFrom int:   "
		<< "double: " << v.as_double() << ", "
		<< "int: " << v.as_int() << ", "
		<< "bool: " << v.as_bool() << ", "
		<< "string: " << v.as_string();
	v = mf("casting/bool").get_value();
	std::cout << "\nFrom bool:  "
		<< "double: " << v.as_double() << ", "
		<< "int: " << v.as_int() << ", "
		<< "bool: " << v.as_bool() << ", "
		<< "string: " << v.as_string();

	std::cout << "\nTesting intvec:\n";
	std::vector<int> ivec = mf("intvec/beta").vec(3).get_int();
	for(size_t i=0; i<ivec.size(); i++)
		std::cout << "ivec[" << i << "] = " << ivec[i] << "\n";
	
	std::cout << "\nTesting mixedvec as double:\n";
	std::vector<double> dvec = mf("mixedvec/vel").vec(3).as_double();
	for(size_t i=0; i<dvec.size(); i++)
		std::cout << "dvec[" << i << "] = " << dvec[i] << "\n";
}
