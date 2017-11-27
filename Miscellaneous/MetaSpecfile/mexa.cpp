/**
 * Implementation of the simple-mexa format in C++. See mexa.h for further comments.
 * 
 * You can compile this file also in standalone-testmode, i.e. with a main() function.
 * Otherwise it will just provide the implementation of the mexa.h library interface.
 * 
 * Compile with 
 * 
 *    g++ -std=c++11 simple-mexa.cpp -Wall -DMEXA_TEST_MAIN
 * 
 * Prepare input file with
 * 
 *   ./mexa.py test.par --outformat simple-mexa > simple-mexa.tmp.par
 * 
 * SvenK, 2017-11-20
 **/

#include "mexa.h"

using namespace mexa;

///////////////////////////////////////////////////////////////////////////////
/////
///// String and vector helper routines
/////
///////////////////////////////////////////////////////////////////////////////


/// string to lowercase
void toLower(std::string& data) {
	std::transform(data.begin(), data.end(), data.begin(), ::tolower);
}

// a buffer-overflow-safe version of sprintf
// source: http://stackoverflow.com/a/69911
std::string vformat (const char *fmt, va_list ap) {
    // Allocate a buffer on the stack that's big enough for us almost
    // all the time.  Be prepared to allocate dynamically if it doesn't fit.
    size_t size = 1024;
    char stackbuf[1024];
    std::vector<char> dynamicbuf;
    char *buf = &stackbuf[0];
    va_list ap_copy;

    while (1) {
        // Try to vsnprintf into our buffer.
        va_copy(ap_copy, ap);
        int needed = vsnprintf (buf, size, fmt, ap);
        va_end(ap_copy);

        // NB. C99 (which modern Linux and OS X follow) says vsnprintf
        // failure returns the length it would have needed.  But older
        // glibc and current Windows return -1 for failure, i.e., not
        // telling us how much was needed.

        if (needed <= (int)size && needed >= 0) {
            // It fit fine so we're done.
            return std::string (buf, (size_t) needed);
        }

        // vsnprintf reported that it wanted to write more characters
        // than we allotted.  So try again using a dynamic buffer.  This
        // doesn't happen very often if we chose our initial size well.
        size = (needed > 0) ? (needed+1) : (size*2);
        dynamicbuf.resize (size);
        buf = &dynamicbuf[0];
    }
}

std::string sformat(const char *fmt, ...) {
	va_list ap;
	va_start (ap, fmt);
	std::string buf = vformat (fmt, ap);
	va_end (ap);
	return buf;
}

/// Split a string into a list of strings based on a single character
// Here is a version which can simply be extended to a list of delimiters: https://stackoverflow.com/a/36563096
stringvec split(const std::string &s, char delim) {
	stringvec elems;
	std::stringstream ss(s); std::string item;
	while (std::getline(ss, item, delim)) if (item.length() > 0) elems.push_back(item);  
	return elems;
}

/// Inplace strip whitespace at beginning and end of string
// and old version not using c++11 lambdas.
void strip(std::string& str) {
	char white_characters[] = " \t\n\r";
	for(char w : white_characters) {
		std::string::iterator end_pos = std::remove(str.begin(), str.end(), w);
		str.erase(end_pos, str.end());
	}
}

// instead several functions from https://stackoverflow.com/a/217605

/// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); })); }
/// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); }).base(), s.end()); }
/// trim from both ends (in place)
static inline void trim(std::string &s) { ltrim(s); rtrim(s); }
/// trim from start (copying)
static inline std::string ltrim_copy(std::string s) { ltrim(s); return s; }
/// trim from end (copying)
static inline std::string rtrim_copy(std::string s) { rtrim(s); return s; }
/// trim from both ends (copying)
static inline std::string trim_copy(std::string s) { trim(s); return s; }


/// Remove empty strings in list of strings
void remove_empty(stringvec& vec) {
	struct StringNotEmpty{bool operator()(const std::string& s) { return !s.empty(); }};
	vec.erase(std::find_if(vec.rbegin(), vec.rend(), StringNotEmpty()).base(), vec.end());
	vec.erase(vec.begin(), std::find_if(vec.begin(), vec.end(), StringNotEmpty()));
}

/// Implode a list of strings with delimiter, with an offset (begin).
std::string join(const stringvec& vec, const char* delim, const stringvec::const_iterator begin, const stringvec::const_iterator end) {
	std::stringstream res;
	std::copy(begin, end, std::ostream_iterator<std::string>(res, delim));
	std::string ret = res.str();
	if(std::strlen(delim)>0)
		ret.erase(ret.end() - std::strlen(delim), ret.end()); // remove last trailing delim...
	return ret;
}
/// Just implode a full vector
std::string join(const stringvec& vec, const char* delim) {
	return join(vec,delim,vec.begin(),vec.end());
}

/// Returns the slice of a vector
template<typename T>
std::vector<T> slice(const std::vector<T>& v, int start=0, int end=-1) {
	int oldlen = v.size();
	int newlen = (end == -1 || end >= oldlen) ? (oldlen-start) : (end-start);
	std::vector<T> nv(newlen);
	for (int i=0; i<newlen; i++) nv[i] = v[start+i];
	return nv;
}

/// Checks wether x is in vector v
template<typename T>
bool contains(const std::vector<T>& v, const T x) {
	return std::find(v.begin(), v.end(), x) != v.end();
}

/// Checks wether x is in map m
template<typename A, typename B>
bool contains(const std::map<A,B>& m, const A x) {
	return m.count(x) > 0; // this is inefficient. improve.
}

/// quickly look into a list of toString()'ables
template<typename T>
void listToString(const T& v) {
	int i=0;
	for(auto c:v){
		printf("%d. %s\n", i++, c.toString().c_str());
	}
}

template<typename T>
std::vector<T> concat(const std::vector<T>& a, const std::vector<T>& b) {
	std::vector<T> c = a;
	c.insert(c.end(), b.begin(), b.end());
	return c;
}
/// whether str starts with search
bool startswith(const std::string& str, const std::string& search) {
	return str.find(search) == 0;
}
std::string remove_common_prefix(std::string a, std::string b) {
	if(b.length() > a.length()) std::swap(a,b);
	return startswith(a,b) ? a.substr(b.length(),a.length()) : a;
}

void stripComment(std::string& line) {
	line.erase( std::find( line.begin(), line.end(), '#' ), line.end() );
}

// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

/// Parse stuff. Result value is true if output is usable and false if
/// it could not succeed.
template<typename T>
bool parse(const std::string& input, T& output) {
	std::istringstream in(input);
	//in.exceptions(std::fstream::failbit);
	in >> output;
	return !(in.fail() || in.bad());
}

/// An ASCII hex '0-9a-FA-F' to integer resolving. Returns <0 in case of error.
char hex2byte(char hex) {
	if('0' <= hex && hex <= '9')
		return hex - '0';
	if('a' <= hex && hex <= 'f')
		return hex - 'a' + 10;
	if('A' <= hex && hex <= 'F')
		return hex - 'A' + 10;
	else	return -1;
}

/**
 * Unescapes (parses) a quoted printable string (think of urlencode). We don't use the
 * MIME typical = escape character (as in =F0) but %F0 because this is less problematic
 * in ExaHyPE specfiles. We encode *all* whitespace, no other rules apply.
 * This implementation uses a super simple finite state machine.
 **/
std::string unescape_quotedprintable(const std::string& input, char escape='%') {
	// Written by SvenK at 2017-11-20
	std::string unescaped;
	enum { E0, E1, E2, ERR } state = E0; // FSM: escape sequence
	char o1,o2; int pos=0;
	for(char c : input) {
		switch(state) {
			case E0:
				if(c == escape) state = E1;
				else unescaped += c;
				break;
			case E1:
				o1 = hex2byte(c);
				state = (o1<0) ? ERR : E2;
				break;
			case E2:
				o2 = hex2byte(c);
				state = (o2<0) ? ERR : E0;
				unescaped += o1*16 + o2;
				break;
			case ERR:
				break;
		}
		if(state != ERR) pos++;
	}
	// error checking
	if(state != E0) {
		std::stringstream errmsg;
		errmsg << "Invalid Character at position " << pos << ", expected hexadecimal specifier %00-%FF.";
		errmsg << " Read so far: '"<<unescaped<<"'"; // debugging
		throw std::runtime_error(errmsg.str());
	}
	return unescaped;
}

///////////////////////////////////////////////////////////////////////////////
/////
///// Class methods
/////
///////////////////////////////////////////////////////////////////////////////



// represents lhs
// class  SYMBOL
	symbol::path_t symbol::parseSymbol(std::string name) {
		toLower(name);
		symbol::path_t path = split(name, '/');
		for(std::string &str : path) strip(str);
		remove_empty(path);
		return path;
	}
	
	/// Canonical representation of symbol
	std::string symbol::toString() const {
		return join(path, "/");
	}
	
	/// Lists all parenting symbols, for instance for /a/b/c it is the list [/, /a, /a/b].
	std::vector<symbol> symbol::ancestors() const {
		std::vector<symbol> ancestors(path.size());
		for(size_t i=0; i<path.size(); i++)
			ancestors[i] = symbol(slice(path,0,i));
		return ancestors;
	}
	
	/// Returns a new symbol which is prefixed by this
	symbol symbol::prefix_add(const symbol& other) const {
		return symbol(concat(path, other.path));
	}
	
	/// Remove common prefix
	symbol symbol::prefix_remove(const symbol& other) const {
		return symbol(remove_common_prefix(toString(), other.toString()));
	}
	
	/// For the usage as key in std::map
	bool symbol::operator <(const symbol& rhs) const {
		return toString() < rhs.toString();
	}
	/// For comparison of equality (contains, ==)
	bool symbol::operator==(const symbol& rhs) const {
		return toString() == rhs.toString();
	}
// end of class symbol

std::string sourcemap::toString() const { return filename + ":" + ::toString(linenumber); }
std::string sourced::toString() const { return std::string("sourced(")+line+","+src.toString()+")"; }


// CLASS mexafile
	mexafile mexafile::query(const symbol root) const {
		// todo
		// Search for all symbols which are *below* the root
		mexafile res;
		for(auto it : assignments) {
			// check case root="foo", include "foo/bar" and "foo/bar/baz"
			// check case root="foo", include "foo" itself.
			if(::contains(it.first.ancestors(), root) || root == it.first) {
				res.assignments[it.first] = it.second;
			}
		}
		return res;
	}
	
	/**
	 * Query and give relative paths to root. Ie. if assignments="a/x=1,a/y=2,b/z=3",
	 * then query(a) = "a/x=1,a/y=2" while query_root(a)="x=1,y=2".
	 **/
	mexafile mexafile::query_root(const symbol root) const {
		return query(root).prefix_remove(root);
	}
	
	mexafile mexafile::prefix_add(const symbol root) const {
		mexafile res;
		for(auto it : assignments) {
			res.assignments[ root.prefix_add(it.first) ] = it.second;
		}
		return res;
	}
	
	mexafile mexafile::prefix_remove(const symbol root) const {
		mexafile res;
		for(auto it: assignments) {
			res.assignments[ root.prefix_remove(it.first) ] = it.second;
		}
		return res;
	}

	std::string mexafile::toString() const {
		std::string ret;
		for(auto it : assignments)
			ret += it.first.toString() + "=" + it.second.toString() + "\n";
		return ret;
	}

	// contains helper
	bool mexafile::contains(symbol leaf, bool doRaise) {
		bool res = ::contains(assignments,leaf);
		if(!res && doRaise) {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' not in assignment list. ";
			if(assignments.empty()) errmsg << "The assignment list is empty";
			else errmsg << "The assignment list is given by: " << toString();
			throw std::runtime_error((errmsg).str());
		}
		return res;
	}
	
	//// getters

	template<typename T>
	T mexafile::get(symbol leaf, std::string type_as_str) {
		T value;
		contains(leaf, true);
		sourced& rhs = assignments[leaf];
		if(!parse(rhs.line, value)) {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' with value '"<< rhs.line << "' cannot be casted as '" << type_as_str << "'. It was given on " <<  rhs.src.toString();
			throw std::runtime_error(errmsg.str());
		}
		return value;
	}
	
	std::string mexafile::get_string(symbol leaf) {
		contains(leaf, true);
		std::string& line = assignments[leaf].line;
		// check for string enclosement characters
		ltrim(line);
		const int quotation_length = 1; // length of the string enclosing quotation: one character
		bool start_quotation_given = (line.find_first_of("'\"") == 0);
		unsigned int end_quotation_position = line.find_first_of(line[0], quotation_length);
		bool end_quotation_given = (end_quotation_position != std::string::npos);
		std::string string_content = line.substr(quotation_length, end_quotation_position-quotation_length);
		std::string rest_of_line = line.substr(end_quotation_position+quotation_length);
		stripComment(rest_of_line);
		
		if(!start_quotation_given) {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must start with \"double\" or 'single' quotation marks. It was given on " <<  assignments[leaf].src.toString();
			throw std::runtime_error(errmsg.str());
		}
		if(!end_quotation_given) {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must end with in \"double\" or 'single' quotation marks. It was given on " <<  assignments[leaf].src.toString();
			throw std::runtime_error(errmsg.str());
		}
		if(!rest_of_line.empty()) {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must be enclosed with in \"double\" or 'single' quotation marks. Only comments can be given afterwards on the same line. It was given on " <<  assignments[leaf].src.toString();
			throw std::runtime_error(errmsg.str());
		}
		return string_content;
	}
	
	bool mexafile::get_bool(symbol leaf) {
		contains(leaf, true);
		std::string value = assignments[leaf].line;
		stripComment(value); toLower(value); strip(value);
		bool isTrue = (value == "yes") || (value == "true") || (value == "on");
		bool isFalse = (value == "no") || (value == "false") || (value == "off");
		if( isTrue && !isFalse) return true;
		if(!isTrue &&  isFalse) return false;
		else {
			std::stringstream errmsg;
			errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as bool, allowed values are only true/yes/on and false/no/off. It was given on " <<  assignments[leaf].src.toString();
			throw std::runtime_error(errmsg.str());
		}
	}
	
	int mexafile::get_int(symbol leaf) { return get<int>(leaf, "int"); }
	
	double mexafile::get_double(symbol leaf) { return get<double>(leaf, "double"); }
// end of class mexafile


///////////////////////////////////////////////////////////////////////////////
/////
///// Functions
/////
///////////////////////////////////////////////////////////////////////////////


mexafile mexa::SimpleMexa(std::istream& fh, const std::string filename_or_desc) {
	mexafile mf;
	std::string line;
	int linecounter = 1; // count for humans
	while(std::getline(fh, line)) {
		sourcemap src(filename_or_desc, linecounter++);
		strip(line);
		if(line.empty() || line[0] == '#') continue;
		stringvec parts = split(line, '=');
		symbol lhs(parts.at(0));
		// We do this crazy joining in case of "=" signs are in the RHS
		std::string rhs = join(parts,"=",parts.begin()+1,parts.end());
		strip(rhs);
		mf.assignments[lhs] = sourced(rhs,src);
	}
	return mf;
}

bool mexa::MexaHasMagicString(std::istream& fh) {
	std::string firstline;
	if(!std::getline(fh, firstline) || firstline.empty())
		return false; // could not even read first line.
	toLower(firstline);
	return startswith(firstline, magicString);
}


mexafile mexa::SimpleMexa_fromEmbedded(const std::string& format, const std::string& content, const std::string filename_or_desc) {
	if(format == "quotedprintable") {
		std::stringstream is(unescape_quotedprintable(content));
		if(!MexaHasMagicString(is)) {
			throw std::runtime_error("Could not detect Mexa file.");
		}
		return SimpleMexa(is, format+":"+filename_or_desc);
	} else {
		throw std::runtime_error("Format not supported so far, only quotedprintable is supported.");
	}
}

///////////////////////////////////////////////////////////////////////////////
/////
///// Main for independent testing
/////
///////////////////////////////////////////////////////////////////////////////

#ifdef MEXA_TEST_MAIN
int main() {
	std::string filename = "simple-mexa.tmp.par";
	std::ifstream infile(filename);
	mexafile mf = SimpleMexa(infile, filename);
	
	std::string query = "exahype/solvers/solver/constants";
	printf("Querying %s:\n", query.c_str());
	printf(mf.query(query).prefix_add("hannes/moritz/").prefix_remove("hannes/").toString().c_str());
	// .query(query)
	// .prefix_add("hannes/moritz")

	// try to get values:
	//printf("back = '%s'\n", mf.query_root("exahype/solvers/solver/constants/boundaries").get_string("back").c_str());
	//printf("back = '%d'\n", mf.query_root("exahype/solvers/solver/constants/boundaries").get_int("back"));
	
	//symbol a("foo/bar/baz");
	//printf("Ancestors of %s:\n", a.toString().c_str());
	//listToString(a.ancestors());
	
	//int z = parse<int>("schlecht");
	//printf("17 = %d\n", z);
	
	// Works
	mexafile mf2 = SimpleMexa_fromEmbedded("quotedprintable","%23%23MEXA%2Dsimple%20configuration%20file%0Avariables%20%3D%20%22rho%3A1%2Cvel%3A3%2CE%3A1%2CB%3A3%2Cpsi%3A1%2Clapse%3A1%2Cshift%3A3%2Cgij%3A6%2Ccoordinates%3A3%2Ccheck%3A1%22%0Atime%20%3D%200%2E0%0Arepeat%20%3D%200%2E001%0Aoutput%20%3D%20%22%2E%2Foutput%2Fglobal%2Dintegrals%22%0Aprimitives%20%3D%20True%0Aconserved%20%3D%20True%0Aerrors%20%3D%20True%0A","test");
	printf(mf2.toString().c_str());
}
#endif 