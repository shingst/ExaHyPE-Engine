#ifndef MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER
#define MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER

// C++ STL:
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <map>
#include <vector>
#include <cstring>
#include <string>
#include <cstdarg>

/**
 * Mexa is a C++11 standalone implementation for the simple-mexa format. This format
 * consists of line-by-line text files in the format
 *
 *   foo/baz = "bar"
 *   foo/bar = 0.9
 *   foo/foo = Yes
 *
 * i.e. it knows hierarchy, assignments and datatypes. In contrast to the fully
 * specified mexa format (which understands inclusion, extension, symbolic variables,
 * abstract trees/lists), this format is very simple.
 * For details about the file format, read the mexa README, for instance at
 * https://bitbucket.org/svek/mexa
 *
 * Written by SvenK in Nov 2017 for ExaHyPE.
 **/
namespace mexa {
	typedef std::vector<std::string> stringvec;

	/**
	 * A symbol represents an hierarchic path. Such a path can point either
	 * to a node in the file or to a leaf where actual information are stored.
	 * This class allows several arithmetic operations with these paths.
	 **/
	struct symbol {
		typedef stringvec path_t;
		path_t path;
		static path_t parseSymbol(std::string name);
		
		symbol() {} ///< the root symbol
		symbol(const std::string& name) : path(parseSymbol(name)) {}
		symbol(const char* const name) : symbol(std::string(name)) {} // for convenience
		symbol(const path_t& path) : path(path) {}
		
		/// Canonical representation of symbol
		std::string toString() const;
		
		/// Lists all parenting symbols, for instance for /a/b/c it is the list [/, /a, /a/b].
		std::vector<symbol> ancestors() const;
		
		/// Returns a new symbol which is prefixed by this
		symbol prefix_add(const symbol& other) const;
		
		/// Remove common prefix
		symbol prefix_remove(const symbol& other) const;
		
		/// For the usage as key in std::map
		bool operator <(const symbol& rhs) const;
		
		/// For comparison of equality (contains, ==)
		bool operator==(const symbol& rhs) const;
	};

	/// The sourcemap represents a line in a file for understandable error messages.
	struct sourcemap {
		std::string filename;
		int linenumber;
		sourcemap(std::string filename, int linenumber) : filename(filename), linenumber(linenumber) {}
		std::string toString() const;
	};

	/// The sourced POD represents a string (line of a file) with a sourcemap information.
	struct sourced {
		std::string line;
		sourcemap src;
		sourced() : src("unknown",-1) {} // for std::map :-/
		sourced(const std::string line, const sourcemap src) : line(line),src(src) {}
		sourced(const sourced&) = default;
		sourced& operator=(const sourced&) = default;
		std::string toString() const;
	};

	/**
	 * The mexafile represents a parameter file in the simple-mexa file format, i.e. it
	 * is only a map of symbols onto lines/strings (with source information for speaking
	 * error messages).
	 *
	 * This is a plain old datatype, it can be constructed empty. There are functions
	 * which return this structure, parsed from a file.
	 *
	 * The class supports various querying mechanisms inspired by XML-Xpath as well as
	 * actual getters for several C++ types. In case of errors, we rely on exceptions.
	 **/
	struct mexafile {
		std::map<symbol, sourced> assignments;

		/**
		 * Returns a string representation which actually is a valid simple-mexa
		 * file itself (if the input is correct).
		 **/
		std::string toString() const;

		/**
		 * Returns true if this file holds the questioned symbol.
		 * @arg doRaise if true, raise an exception if symbol not given
		 **/
		bool contains(symbol leaf, bool doRaise=false);

		///////////////////////////////
		// querying
		///////////////////////////////
		
		/**
		 * Filter the file for a given symbol. Returns a new mexafile.
		 **/
		mexafile query(const symbol root) const;
		
		/**
		* Query and give relative paths to root. Ie. if assignments="a/x=1,a/y=2,b/z=3",
		* then query(a) = "a/x=1,a/y=2" while query_root(a)="x=1,y=2".
		**/
		mexafile query_root(const symbol root) const;

		/**
		 * Returns a copy of this mexafile where the symbol `root' is added to
		 * every entry in the file.
		 **/
		mexafile prefix_add(const symbol root) const;

		/**
		 * Returns a copy of this mexafile where the symbol `root' is subtracted
		 * from every entry in the file.
		 **/
		mexafile prefix_remove(const symbol root) const;
		
		///////////////////////////////
		// getters
		///////////////////////////////

		/// An abstract getter for a templated type
		template<typename T>
		T get(symbol leaf, std::string type_as_str);
		
		/// Read a string at the questioned symbol
		std::string get_string(symbol leaf);
		
		/// Read a boolean value at the questioned symbol
		bool get_bool(symbol leaf);
		
		/// Read an integer
		int get_int(symbol leaf);
		
		/// Read a floating point value (double).
		double get_double(symbol leaf);
	};


	/**
	* A parser for a simple mexa file. You pass an input stream and we give you
	* a mexafile instance.
	* 
	* @param filename_or_desc Used in case of errors to trace back the errnous input.
	* 
	* (One could also come up with a parser for the full mexa file format.)
	**/
	mexafile SimpleMexa(std::istream& fh, const std::string filename_or_desc="unknown");

	/// A magic byte you can/should/might want to put at the beginning of a mexa file.
	constexpr const char* magicString = "##mexa";
	
	/// Reads the first line of a (potential) Mexa string and checks for the magic string.
	/// This line is then removed from the istream unless you rewind.
	bool MexaHasMagicString(std::istream& fh);

	/**
	 * Read a simple-mexa file from an embedded, i.e. encoded, place. Typical formats
	 * are base64, base58, base16 as well as `quoted printable' (pass "quotedprintable")
	 * where all textual parts are readable and only non-alphanumeric characters are
	 * escaped.
	 * 
	 * This is your entry point if you want to read mexa files out of ExaHyPE specification
	 * files (embedding).
	 **/
	mexafile SimpleMexa_fromEmbedded(
		const std::string& format, const std::string& content,
		const std::string filename_or_desc="unknown");
	
} // ns mexa

#endif /* MEXA_META_EXAHYPE_PARFILE_CPP_IMPL_HEADER */