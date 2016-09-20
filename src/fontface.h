#ifndef FONTFACE_H
#define FONTFACE_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "jsoncons/json.hpp"

namespace gen {

struct TextExtents {
    double offx = 0.0;
    double offy = 0.0;
    double width = 0.0;
    double height = 0.0;
    double dx = 0.0;
    double dy = 0.0;
};

class FontFace {

public:
    FontFace();
    FontFace(std::string datafilename);

    std::string getFontFace();
    bool setFontFace(std::string name);
    bool setFontFace(std::string name, int size);
    bool setFontFace();
    std::vector<std::string> getFontFaces();
    std::vector<int> getFontSizes();
    std::vector<int> getFontSizes(std::string font);
    int getFontSize();
    bool setFontSize(int size);
    bool setFontSize();
    TextExtents getTextExtents(std::string str);
    std::vector<TextExtents> getCharacterExtents(std::string str);

private:
    std::string _intToString(int number);
    int _stringToInt(std::string number);
    TextExtents _getCharExtents(char c);

    jsoncons::json _jsonData;

    std::string _defaultFont;
    std::string _defaultFontSize;
    std::string _fontFace;
    std::string _fontSize;
};

}

#endif