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

///////////////////////////////////////////////////////////////////////////////
/////
///// String and vector helper routines
/////
///////////////////////////////////////////////////////////////////////////////

namespace mexa {
namespace tools {
	
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
	if(vec.empty()) return "";
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


// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
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

} // ns tools
} // ns mexa

using namespace mexa;
using namespace tools; // mexa::tools

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

// My poor man's value variant type
bool mexa::value::canCastTo(Type type) const {
	bool targetIsNumerical = (type == mexa::value::Type::INT) || (type == mexa::value::Type::DOUBLE) || (type == mexa::value::Type::BOOL);
	bool sourceIsNumerical = (a == mexa::value::Type::INT) || (a == mexa::value::Type::DOUBLE) || (a == mexa::value::Type::BOOL);
	return targetIsNumerical && sourceIsNumerical;
}

template<typename T>
T mexa::value::cast_as_numerical(mexa::value::Type type) const {
	switch(a) {
		case mexa::value::Type::INT: return get_int();
		case mexa::value::Type::DOUBLE: return get_double();
		case mexa::value::Type::BOOL: return get_bool();
		default: break;
	}
	// value cast error:
	std::stringstream errmsg;
	errmsg << "Value '"<<toString()<< "' not castable to " << mexa::value::type2str(type);
	throw std::runtime_error(errmsg.str());
}

int mexa::value::as_int() const { return cast_as_numerical<int>(mexa::value::Type::INT); }
double mexa::value::as_double() const { return cast_as_numerical<double>(mexa::value::Type::DOUBLE); }
bool mexa::value::as_bool() const { return cast_as_numerical<bool>(mexa::value::Type::BOOL); }

std::string mexa::value::as_string() const {
	std::stringstream ret;
	switch(a) {
		case mexa::value::Type::INT: ret << get_int(); break;
		case mexa::value::Type::DOUBLE: ret << get_double(); break;
		case mexa::value::Type::BOOL: ret << get_bool(); break;
		case mexa::value::Type::STRING: ret << get_string(); break;
		case mexa::value::Type::UNDEF: ret << "---"; break;
		default: ret << "error";
	}
	return ret.str();
}

std::string mexa::value::type2str(mexa::value::Type type) {
	switch(type) {
		case mexa::value::Type::INT: return "int";
		case mexa::value::Type::DOUBLE: return "double";
		case mexa::value::Type::BOOL: return "bool";
		case mexa::value::Type::STRING: return "string";
		case mexa::value::Type::UNDEF: return "undefined";
		default: return "error";
	}
}

void mexa::value::assertActive(mexa::value::Type requested) const {
	if(!isActive(requested)) {
		std::stringstream errmsg;
		errmsg << "Value '"<<toString()<< "' does not hold a " << type2str(requested);
		throw std::runtime_error(errmsg.str());
	}
}

std::string mexa::value::toString() const {
	std::stringstream ret;
	ret << type2str(a) << "(" << as_string() << ")";
	return ret.str();
}

// VECTOR VALUES

vector_value::vector_value(mexafile mq, size_t required_length) : mq(mq), required_length(required_length) {
	if(mq.isEmpty() && required_length != no_required_length) {
		std::stringstream errmsg;
		errmsg << "The given assignment list is empty:" << mq.toString();
		throw std::runtime_error(errmsg.str());
	}
}

template<typename T>
std::vector<T> vector_value::get( T (value::*getter)() const , value::Type type, bool doCast) const {
	std::vector<T> ret;
	int itemcount = 1; // just for error output
	for(auto it : mq.assignments()) {
		bool isExact = it.val.isActive(type);
		bool canCast = it.val.canCastTo(type);
		if((doCast && !canCast) || (!doCast && !isExact)) {
			std::stringstream errmsg;
			errmsg << "While reading in the " << value::type2str(type) << "-vector, the " << itemcount << ". symbol " << it.key.toString() << " with value '" << it.val.toString() << "' cannot be casted as " << value::type2str(type) << ". It was given on " << it.src.toString();
			throw std::runtime_error(errmsg.str());
		}
		T unpacked_value = (it.val.*getter)(); // should call e.g. value::get_int for T=int.
		ret.push_back(unpacked_value);
		itemcount++;
	}
	// assert_vec_size(node, mq, type_as_str, require_length, ret);
	if(required_length != no_required_length) { // asked for a require_length
                if(ret.size() != required_length) {
			std::stringstream errmsg;
			errmsg << "Having read in a " << value::type2str(type) << "-vector with size "<<ret.size()<<", but a vector of size "<< required_length
			<< " is required. The vector values are given by: [";
			for(auto j : ret) errmsg << j << ",";
			errmsg << "] and where read by " << mq.toString();
			throw std::runtime_error(errmsg.str());
	        }
        }

	return ret;
}

std::vector<int> vector_value::get_int() const { return get(&value::get_int, value::Type::INT, false); }
std::vector<double> vector_value::get_double() const { return get(&value::get_double, value::Type::DOUBLE, false); }
std::vector<bool> vector_value::get_bool() const { return get(&value::get_bool, value::Type::BOOL, false); }
std::vector<std::string> vector_value::get_string() const { return get(&value::get_string, value::Type::STRING, false); }
std::vector<int> vector_value::as_int() const { return get(&value::as_int, value::Type::INT, true); }
std::vector<double> vector_value::as_double() const { return get(&value::as_double, value::Type::DOUBLE, true); }
std::vector<bool> vector_value::as_bool() const { return get(&value::as_bool, value::Type::BOOL, true); }
std::vector<std::string> vector_value::as_string() const { return get(&value::as_string, value::Type::STRING, true); }


// The multiline string is probably something we could reintroduce if required.
/*	
	// specialization
	std::string mexafile::get_multiline_string(symbol node) const {
		std::string ret;
		mexa::mexafile mq = query(node);
		for(auto it : mq.assignments()) {
			symbol& leaf(it.first);
			ret += get_string(leaf) + "\n";
		}
		return ret;
	}*/

// CLASS assignment:
assignment assignment::make(symbol key, value val, sourcemap src) {
	assignment ret;
	ret.key = key;
	ret.val = val;
	ret.src = src;
	return ret;
}

std::string assignment::toString() const {
	std::stringstream ret;
	ret << "eq(",
	ret << key.toString();
	ret << ", ",
	ret << val.toString();
	ret << ", source=";
	ret << src.toString();
	ret << ")";
	return ret.str();
}

// CLASS mexafile
void mexafile::add(const symbol key, const value val, const sourcemap src) {
	assignments().push_back(assignment::make(key,val,src));
}

value mexafile::get() const {
	if(assignments().empty()) {
		std::stringstream errmsg;
		errmsg << "Assignment list is empty, cannot get an element: " << toString();
		throw std::runtime_error((errmsg).str());
	} else if(assignments().size() >1) {
		std::stringstream errmsg;
		errmsg << "Assignment list contains more then one element, cannot get a single element: " << toString();
		throw std::runtime_error((errmsg).str());
	} else if(!(assignments()[0].key == symbol())) {
		std::stringstream errmsg;
		errmsg << "Assignment list does not contain a root. You probably request a leaf where there is only a node. The assignment list is given by: " << toString();
		throw std::runtime_error((errmsg).str());
	}
	
	return assignments()[0].val;
}

vector_value mexafile::vec(size_t required_length) const {
	return vector_value(*this, required_length);
}

value mexafile::get(const symbol key) const { return query(key).get(); }
vector_value mexafile::vec(const symbol node, size_t required_length) const { return query(node).vec(); }

mexafile mexafile::filter(const symbol root) const {
	// todo
	// Search for all symbols which are *below* the root
	mexafile res;
	res.source = source + ".filter(\"" + root.toString() + "\")";
	for(auto it : assignments()) {
		// check case root="foo", include "foo/bar" and "foo/bar/baz"
		// check case root="foo", include "foo" itself.
		if(::contains(it.key.ancestors(), root) || root == it.key) {
			res.assignments().push_back(it);
		}
	}
	return res;
}

const mexafile& mexafile::requireNonEmpty() const {
	if(isEmpty()) {
		std::stringstream errmsg;
		errmsg << "The mexafile is empty: "<< toString();
		throw std::runtime_error((errmsg).str());
	}
	return *this;
}

/**
 * Query and give relative paths to root. Ie. if assignments="a/x=1,a/y=2,b/z=3",
 * then query(a) = "a/x=1,a/y=2" while query_root(a)="x=1,y=2".
 **/
mexafile mexafile::query(const symbol root) const {
	mexafile res = filter(root).prefix_remove(root);
	res.source = source + ".query(\"" + root.toString() + "\")";
	return res;
}

mexafile mexafile::prefix_add(const symbol root) const {
	mexafile res;
	res.source = source + ".prefix_add(\"" + root.toString() + "\")";
	for(auto it : assignments())
		res.add(root.prefix_add(it.key), it.val, it.src);
	return res;
}

mexafile mexafile::prefix_remove(const symbol root) const {
	mexafile res;
	res.source = source + ".prefix_remove(\"" + root.toString() + "\")";
	for(auto it: assignments())
		res.add(root.prefix_remove(it.key), it.val, it.src);
	return res;
}

std::string mexafile::toString() const {
	std::string ret = source;
	if(!isEmpty()) {
		ret += " = (\n";
		for(auto it : assignments())
			ret += "\t" + it.toString() + ",\n";
		ret +=")\n";
	} else {
		ret += " = (empty)";
	}
	return ret;
}

// contains helper
bool mexafile::contains(symbol leaf, bool doRaise) const {
	bool res = false;
	for(auto it : assignments())
		if(leaf == it.key)
			res = true;
	if(!res && doRaise) {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' not in assignment list. ";
		if(assignments().empty()) errmsg << "The assignment list is empty";
		else errmsg << "The assignment list is given by: " << toString();
		throw std::runtime_error((errmsg).str());
	}
	return res;
}
	
//// PARSERS

void parser::stripComment(std::string& line) {
	line.erase( std::find( line.begin(), line.end(), '#' ), line.end() );
}

template<typename T> bool parser::tryCast(T& output) const {
	std::istringstream in(input);
	//in.exceptions(std::fstream::failbit);
	in >> output;
	bool successfulCast = !(in.fail() || in.bad());
	// Make sure nothing follows the casted thing.
	std::string remnant;
	in >> remnant;
	strip(remnant); stripComment(remnant); strip(remnant);
	return successfulCast && remnant.empty();
}

template<typename T> T parser::get_as(std::string type_as_str) const {
	T value;
	if(!tryCast(value)) {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' with value '"<< input << "' cannot be casted as '" << type_as_str << "'. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
	return value;
}

bool parser::is_a_bool() const {
	std::string value = input;
	stripComment(value); toLower(value); strip(value);
	bool isTrue = (value == "yes") || (value == "true") || (value == "on");
	bool isFalse = (value == "no") || (value == "false") || (value == "off");
	bool meansTrue  = ( isTrue && !isFalse);
	bool meansFalse = (!isTrue &&  isFalse);
	return (meansTrue || meansFalse);
}

bool parser::get_as_bool() const {
	std::string value = input;
	stripComment(value); toLower(value); strip(value);
	bool isTrue = (value == "yes") || (value == "true") || (value == "on");
	bool isFalse = (value == "no") || (value == "false") || (value == "off");
	if( isTrue && !isFalse) return true;
	if(!isTrue &&  isFalse) return false;
	else {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as bool, allowed values are only true/yes/on and false/no/off. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
}

int parser::get_as_int() const { return get_as<int>("int"); }

bool parser::is_a_int() const {
	int val; return tryCast(val);
}

double parser::get_as_double() const { return get_as<double>("double"); }

bool parser::is_a_double() const {
	double val; return tryCast(val);
}
std::string parser::get_as_quoted_string() const {
	std::string line = input;
	// check for string enclosement characters
	ltrim(line);
	const int quotation_length = 1; // length of the string enclosing quotation: one character
	bool start_quotation_given = (line.find_first_of("'\"") == 0);
	size_t end_quotation_position = line.find_first_of(line[0], quotation_length);
	bool end_quotation_given = (end_quotation_position != std::string::npos);
	std::string string_content = line.substr(quotation_length, end_quotation_position-quotation_length);
	std::string rest_of_line = line.substr(end_quotation_position+quotation_length);
	stripComment(rest_of_line);
	
	if(!start_quotation_given) {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must start with \"double\" or 'single' quotation marks. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
	if(!end_quotation_given) {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must end with in \"double\" or 'single' quotation marks. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
	if(!rest_of_line.empty()) {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' cannot be casted as string. Strings must be enclosed with in \"double\" or 'single' quotation marks. Only comments can be given afterwards on the same line. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
	return string_content;
}

bool parser::is_a_quoted_string() const {
	std::string line = input;
	// check for string enclosement characters
	ltrim(line);
	const int quotation_length = 1; // length of the string enclosing quotation: one character
	bool start_quotation_given = (line.find_first_of("'\"") == 0);
	size_t end_quotation_position = line.find_first_of(line[0], quotation_length);
	bool end_quotation_given = (end_quotation_position != std::string::npos);
	std::string string_content = line.substr(quotation_length, end_quotation_position-quotation_length);
	std::string rest_of_line = line.substr(end_quotation_position+quotation_length);
	stripComment(rest_of_line);
	
	return start_quotation_given && end_quotation_given && rest_of_line.empty();
}

value parser::getValue() const {
	if(is_a_int())
		return value(get_as_int());
	if(is_a_double())
		return value(get_as_double());
	if(is_a_bool())
		return value(get_as_bool());
	if(is_a_quoted_string())
		return value(get_as_quoted_string());
	else {
		std::stringstream errmsg;
		errmsg << "Leaf '" << leaf.toString() << "' with value '" << input << "' cannot be understood neither as int, double, bool nor as a String. Strings must be enclosed with in \"double\" or 'single' quotation marks. Only comments can be given afterwards on the same line. It was given on " <<  src.toString();
		throw std::runtime_error(errmsg.str());
	}
}

value parser::getSloppyValue() const {
	if(is_a_int())
		return value(get_as_int());
	if(is_a_double())
		return value(get_as_double());
	if(is_a_bool())
		return value(get_as_bool());
	// in any other case, treat this as a string.
	return value(input);
}



///////////////////////////////////////////////////////////////////////////////
/////
///// Functions
/////
///////////////////////////////////////////////////////////////////////////////

mexafile mexa::fromFile(std::istream& fh, const std::string filename_or_desc) {
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
		mf.add(lhs, parser(lhs,rhs,src).getValue(), src);
	}
	mf.source = "mexa::fromFile(" + filename_or_desc + ")";
	return mf;
}

bool mexa::hasMagicString(std::istream& fh, const std::string magicString) {
	std::string firstline;
	if(!std::getline(fh, firstline) || firstline.empty())
		return false; // could not even read first line.
	toLower(firstline);
	return startswith(firstline, magicString);
}

mexafile mexa::fromEmbedded(const std::string& format, const std::string& content, const std::string filename_or_desc) {
	if(format == "quotedprintable") {
		std::stringstream is(unescape_quotedprintable(content));
		if(!hasMagicString(is)) {
			throw std::runtime_error("Could not detect Mexa file.");
		}
		return fromFile(is, format+":"+filename_or_desc);
	} else {
		throw std::runtime_error("Format not supported so far, only quotedprintable is supported.");
	}
}

mexafile mexa::fromSpecfile(const std::vector<std::pair<std::string, std::string>>& list, const std::string source_description) {
	mexafile mf;
	int linecounter = 1; // counting like humans 
	for(auto it : list) {
		sourcemap src(source_description, linecounter++);
		symbol lhs(it.first);
		std::string rhs(it.second);
		mf.add(lhs, parser(lhs, rhs, src).getSloppyValue(), src);
	}
	mf.source = "mexa::fromSpecfile(" + source_description + ")";
	return mf;
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
	mexafile mf = mexa::fromFile(infile, filename);
	
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
	mexafile mf2 = mexa::fromEmbedded("quotedprintable","%23%23MEXA%2Dsimple%20configuration%20file%0Avariables%20%3D%20%22rho%3A1%2Cvel%3A3%2CE%3A1%2CB%3A3%2Cpsi%3A1%2Clapse%3A1%2Cshift%3A3%2Cgij%3A6%2Ccoordinates%3A3%2Ccheck%3A1%22%0Atime%20%3D%200%2E0%0Arepeat%20%3D%200%2E001%0Aoutput%20%3D%20%22%2E%2Foutput%2Fglobal%2Dintegrals%22%0Aprimitives%20%3D%20True%0Aconserved%20%3D%20True%0Aerrors%20%3D%20True%0A","test");
	printf(mf2.toString().c_str());
}
#endif 
