#include "fontface.h"

gen::FontFace::FontFace() {
	
}

gen::FontFace::FontFace(std::string datafilename) : _defaultFont("Arial") {
	std::ifstream file(datafilename);
	std::string jsonstr((std::istreambuf_iterator<char>(file)),
	                     std::istreambuf_iterator<char>());
	jsoncons::json json = jsoncons::json::parse(jsonstr);

	auto it = json.find(_defaultFont);
	if (it != json.members().end()) {
	    _fontFace = _defaultFont;
	} else {
		for (const auto& member : json.members()) {
		    _fontFace = member.name();
		    _defaultFont = member.name();
		    break;
		}
	}

	for (const auto& member : json[_fontFace].members()) {
	    _fontSize = member.name();
	    _defaultFontSize = member.name();
	    break;
	}

	_jsonData = json;
}

std::string gen::FontFace::getFontFace() {
	return _fontFace;
}

bool gen::FontFace::setFontFace(std::string name) {
	auto it = _jsonData.find(name);
	if (it == _jsonData.members().end()) {
	    return false;
	}

	_fontFace = name;
	it = _jsonData[_fontFace].find(_fontSize);
	if (it == _jsonData[_fontFace].members().end()) {
	    for (const auto& member : _jsonData[_fontFace].members()) {
		    _fontSize = member.name();
		    break;
		}
	}

	return true;
}

bool gen::FontFace::setFontFace(std::string name, int size) {
	auto it = _jsonData.find(name);
	if (it == _jsonData.members().end()) {
	    return false;
	}

	std::string sizestr = _intToString(size);
	it = _jsonData[_fontFace].find(sizestr);
	if (it == _jsonData[_fontFace].members().end()) {
	    return false;
	}

	_fontFace = name;
	_fontSize = sizestr;

	return true;
}

bool gen::FontFace::setFontFace() {
	_fontFace = _defaultFont;
	_fontSize = _defaultFontSize;

	return true;
}

std::vector<std::string> gen::FontFace::getFontFaces() {
	std::vector<std::string> fonts;
	for (const auto& member : _jsonData.members()) {
	    fonts.push_back(member.name());
	}

	return fonts;
}

std::vector<int> gen::FontFace::getFontSizes() {
	std::vector<int> sizes;
	for (const auto& member : _jsonData[_fontFace].members()) {
	    sizes.push_back(_stringToInt(member.name()));
	}

	return sizes;
}

std::vector<int> gen::FontFace::getFontSizes(std::string font) {
	auto it = _jsonData.find(font);
	if (it == _jsonData.members().end()) {
	    return std::vector<int>();
	}

	std::vector<int> sizes;
	for (const auto& member : _jsonData[font].members()) {
	    sizes.push_back(_stringToInt(member.name()));
	}

	return sizes;
}

int gen::FontFace::getFontSize() {
	return _stringToInt(_fontSize);
}

bool gen::FontFace::setFontSize(int size) {
	std::string sizestr = _intToString(size);
	auto it = _jsonData[_fontFace].find(sizestr);
	if (it == _jsonData[_fontFace].members().end()) {
	    return false;
	}
	_fontSize = sizestr;

	return true;
}

bool gen::FontFace::setFontSize() {
	for (const auto& member : _jsonData[_fontFace].members()) {
	    _fontSize = member.name();
	    return true;
	}

	return false;
}

gen::TextExtents gen::FontFace::getTextExtents(std::string str) {
	TextExtents extents;
	if (str.size() == 0) {
		return extents;
	}

	if (str.size() == 1) {
		return _getCharExtents(str[0]);
	}

	TextExtents temp = _getCharExtents(str[0]);
	extents.offx = temp.offx;

	double ymin = 0.0;
	double ymax = 0.0;
	double dx = 0.0;
	for (unsigned int i = 0; i < str.size(); i++) {
		temp = _getCharExtents(str[i]);
	    ymin = fmin(ymin, temp.offy);
	    ymax = fmax(ymax, temp.offy + temp.height);
	    dx += temp.dx;
	}
	extents.offy = ymin;
	
	temp = _getCharExtents(str[str.size() - 1]);
	extents.width = dx + extents.offx - (temp.dx - temp.width);
	extents.height = ymax - ymin;
	extents.dx = dx;
	extents.dy = 0.0;

	return extents;
}

gen::TextExtents gen::FontFace::_getCharExtents(char c) {
	std::string key(1, c);
	jsoncons::json chardata = _jsonData[_fontFace][_fontSize][key];
	double data[6];

	int i = 0;
	for (const auto& member : chardata.elements()) {
	    data[i] = member.as<double>();
	    i++;

	    if (i >= 6) {
	    	break;
	    }
	}

	TextExtents extents;
	extents.offx = data[0];
	extents.offy = data[1];
	extents.width = data[2];
	extents.height = data[3];
	extents.dx = data[4];
	extents.dy = data[5];

	return extents;
}

std::string gen::FontFace::_intToString(int number) {
    std::ostringstream ss;
    ss << number;
    return ss.str();
}

int gen::FontFace::_stringToInt(std::string number) {
	int s = 0;
	std::istringstream(number) >> s;
	return s;
}
