#ifndef __GRMHD_RUNTIME_PARAMETERS__
#define __GRMHD_RUNTIME_PARAMETERS__

/*
 * RuntimeParameters was an approach to parse strings such as
 * 
 *   {foo:bar,baz:bla}
 *
 * i.e. select statements or untokenized strings. We do not encounter such
 * strings any more and the successor of this technique is the mexa.h
 * interface for hierarchical data.
 *
 */

namespace RuntimeParameters {

// abstract class
struct MapView {
	virtual bool is_string(const std::string& key) const = 0;
	virtual std::string get_string(const std::string& key) const = 0;
	
	// todo: Add also other types
};

/// Translate to ExaHyPE's language
template<class ParserViewOrSo>
struct ParserView : public MapView {
	ParserViewOrSo& p;
	ParserView(ParserViewOrSo& _p) : p(_p) {}
	bool is_string(const std::string& key) const override {
		return p.isValueValidString(key);
	}
	std::string get_string(const std::string& key) const override {
		return p.getValueAsString(key);
	}
};

/** Helper class to parse key value strings */
class StringMapView : public MapView {
   std::string base;
   public:
    StringMapView(std::string _base) : base(_base) {}
    std::string get_string(const std::string& key) const override {
	std::size_t startIndex = base.find(key);
	startIndex = base.find(":", startIndex);
	std::size_t endIndex = base.find_first_of("}, \n\r", startIndex + 1);
	return base.substr(startIndex + 1, endIndex - startIndex - 1);
    }
    bool is_string(const std::string& key) const override {
	std::size_t startIndex = base.find(key);
	return (startIndex != std::string::npos);
    }
};

// we will add here shortly:
// a) a correct specfile parserview
// b) a mexa map view

} // ns RuntimeParameters

#endif /* __GRMHD_RUNTIME_PARAMETERS__ */
